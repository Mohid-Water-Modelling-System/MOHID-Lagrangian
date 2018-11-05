    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the simulation class and its methods. This is the only
    !> class that is exposed to an external program, as it encapsulates every other
    !> class and method.
    !------------------------------------------------------------------------------
    module simulation_mod

    use simulation_about_mod
    use simulation_initialize_mod
    use boundingbox_mod
    use blocks_mod
    use emitter_mod
    use sources_mod
    use tracers_mod
    use background_mod
    use simulation_output_streamer_mod
    use common_modules

    implicit none
    private

    type :: simulation_class   !< Parameters class
        type(timer_class) :: timerInit          !< timer for the initialization routines
        type(timer_class) :: timerTotalRun      !< timer for the total wall time spent on the simulation
        type(timer_class) :: timerOutput        !< timer for the output writting routines
    contains
    procedure, public  :: initialize => initSimulation
    procedure, public  :: run
    procedure, public  :: finalize   => closeSimulation
    procedure, private :: decompose  => DecomposeDomain
    procedure, private :: ToggleSources
    procedure, private :: BlocksEmitt
    procedure, private :: BlocksDistribute
    procedure, private :: BlocksConsolidateArrays
    procedure, private :: BlocksTracersToAoT
    procedure, private :: BlocksRunSolver
    procedure, private :: BlocksAoTtoTracers
    procedure, private :: BlocksCleanAoT
    procedure, private :: setInitialState
    procedure, private :: getTracerTotals
    procedure, private :: printTracerTotals
    procedure, private :: setTracerMemory
    end type simulation_class

    !Exposed public class
    public :: simulation_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation run method. Runs the initialized case main time cycle.
    !---------------------------------------------------------------------------
    subroutine run(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext, aux

    outext = '====================================================================='
    call Log%put(outext,.false.)
    outext = '->Simulation starting'
    call Log%put(outext)
    outext = '====================================================================='
    call Log%put(outext,.false.)

    !main time cycle
    do while (Globals%SimTime .lt. Globals%Parameters%TimeMax)
        call Globals%Sim%increment_numdt()
        call self%timerTotalRun%Tic()
        !activate suitable Sources
        call self%ToggleSources()
        !emitt Tracers from active Sources
        call self%BlocksEmitt()
        !Distribute Tracers and Sources by Blocks
        call self%BlocksDistribute()
        !Optimize Block Tracer lists
        call self%BlocksConsolidateArrays()
        !Build AoT
        call self%BlocksTracersToAoT()
        !load hydrodynamic fields from files (curents, wind, waves, ...)
        !Update all tracers with base behavior (AoT) - Integration step
        call self%BlocksRunSolver()
        !AoT to Tracers
        call self%BlocksAoTtoTracers()
        !Update Tracers with type-specific behavior
        !Write results if time to do so
        call self%timerOutput%Tic()
        call OutputStreamer%WriteStepSerial(DBlock)
        call self%timerOutput%Toc()
        !Print some stats from the time step
        call self%printTracerTotals()
        !Clean AoT
        call self%BlocksCleanAoT()
        !update Simulation time and counters
        Globals%SimTime = Globals%SimTime + Globals%SimDefs%dt
        !print*, 'Global time is ', Globals%SimTime
        !print*, 'Can we continue?'
        !read (*,*)
        call self%timerTotalRun%Toc()
    enddo
    call self%setTracerMemory()
    call SimMemory%detailedprint()

    call self%timerTotalRun%print()
    call self%timerOutput%print()

    end subroutine run

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation initialization method. Effectively builds and populates the
    !> simulation objects that will be used latter on.
    !> @param[in] self, casefilename, outpath
    !---------------------------------------------------------------------------
    subroutine initSimulation(self, casefilename, outpath)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string), intent(in) :: casefilename         !< case file name
    type(string), intent(in) :: outpath              !< Output path
    type(string) :: outext, aux
    !type(generic_field_class) :: testField
    !type(background_class) :: testBackground
    
    aux = 'Simulation::initialization'
    call self%timerInit%initialize(aux)
    call self%timerInit%Tic()

    aux = 'Simulation::Run'
    call self%timerTotalRun%initialize(aux)
    aux = 'Simulation::Output'
    call self%timerOutput%initialize(aux)

    ! Initialize logger
    call Log%initialize(outpath)
    !Print licences and build info
    call PrintLicPreamble
    !initialize parallel services manager
    call OMPManager%initialize()
    !initializing memory log
    call SimMemory%initialize()
    !setting every global variable and input parameter to their default
    call Globals%initialize(outpath = outpath)
    !initializing geometry class
    call Geometry%initialize()
    !Check if case file has .xml extension
    if (casefilename%extension() == '.xml') then
        ! Initialization routines to build the simulation from the input case file
        call InitFromXml(casefilename)
    else
        outext='[initSimulation]: only .xml input files are supported at the time. Stopping'
        call Log%put(outext)
        stop
    endif
    !Case was read and now we can build/initialize our simulation objects that are case-dependent
    !initilize simulation bounding box
    call BBox%initialize()
    !decomposing the domain and initializing the Simulation Blocks
    call self%decompose()
    !Distributing Sources
    call self%setInitialState()
    !printing memory occupation at the time
    call SimMemory%detailedprint()
    !Initializing output file streamer
    call OutputStreamer%initialize()
    !Writing the domain to file
    call OutputStreamer%WriteDomain(Globals%Names%casename, BBox, Geometry%getnumPoints(BBox), DBlock)
    
    call self%timerInit%Toc()
    call self%timerInit%print()

    !call testField%test()
    !call testBackground%test()

    end subroutine initSimulation

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to activate and deactivate Sources based on the
    !> Global%SimTime
    !---------------------------------------------------------------------------
    subroutine ToggleSources(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(DBlock)
        call DBlock(i)%ToogleBlockSources()
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine ToggleSources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to emitt tracers at current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksEmitt(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    do i=1, size(DBlock)
        call DBlock(i)%CallEmitter()
    enddo
    end subroutine BlocksEmitt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to distribute Tracers at
    !> current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksDistribute(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    do i=1, size(DBlock)
        call DBlock(i)%DistributeTracers()
    enddo
    !need to distribute Sources also! TODO
    end subroutine BlocksDistribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to consolidate the Tracer array at
    !> current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksConsolidateArrays(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(DBlock)
        call DBlock(i)%ConsolidateArrays()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine BlocksConsolidateArrays

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to build their Array of
    !> Tracers (AoT) from the Tracer list at current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksTracersToAoT(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(DBlock)
        call DBlock(i)%TracersToAoT()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine BlocksTracersToAoT

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to run the Solver at
    !> current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksRunSolver(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(DBlock)
        call DBlock(i)%RunSolver()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine BlocksRunSolver

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to print their Array of
    !> Tracers (AoT) back to the Tracer objects on the list at current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksAoTtoTracers(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(DBlock)
        call DBlock(i)%AoTtoTracers()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine BlocksAoTtoTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to clean their Array of
    !> Tracers (AoT) at current SimTime
    !---------------------------------------------------------------------------
    subroutine BlocksCleanAoT(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i
    do i=1, size(DBlock)
        call DBlock(i)%CleanAoT()
    enddo
    end subroutine BlocksCleanAoT

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to distribute the Sources to the Blocks, allocate the
    !> respective Tracers and redistribute if needed
    !---------------------------------------------------------------------------
    subroutine setInitialState(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext
    integer :: i, blk, ntrc
    !iterate every Source to distribute
    ntrc = 0
    do i=1, size(tempSources%src)
        blk = getBlockIndex(Geometry%getCenter(tempSources%src(i)%par%geometry))
        call DBlock(blk)%putSource(tempSources%src(i))
        ntrc = ntrc + tempSources%src(i)%stencil%total_np
    end do
    call tempSources%finalize() !destroying the temporary Sources now they are shipped to the Blocks
    outext='-->Sources allocated to their current Blocks'
    call Log%put(outext,.false.)
    outext = ntrc
    outext='-->'//outext//' Tracers on the emission stack'
    call Log%put(outext,.false.)
    call self%setTracerMemory(ntrc)
    end subroutine setInitialState

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to count Tracer numbers
    !---------------------------------------------------------------------------
    integer function getTracerTotals(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: i, total
    total = 0
    do i=1, size(DBlock)
        total = total + DBlock(i)%numAllocTracers()
    enddo
    getTracerTotals = total
    end function getTracerTotals

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to count Tracer numbers
    !---------------------------------------------------------------------------
    subroutine printTracerTotals(self)
    implicit none
    class(simulation_class), intent(in) :: self
    type(string) :: outext, temp
    temp = self%getTracerTotals()
    outext='-->'//temp //' Tracers allocated'
    call Log%put(outext,.false.)
    end subroutine printTracerTotals

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to account for Tracer memory consumption
    !---------------------------------------------------------------------------
    subroutine setTracerMemory(self, ntrc)
    implicit none
    class(simulation_class), intent(in) :: self
    integer, optional, intent(in) :: ntrc
    integer :: sizem, i
    type(tracer_class) :: dummyTracer
    sizem = 0
    do i=1, size(DBlock)
        sizem = sizem + sizeof(DBlock(i)%LTracer) !this accounts for the array structure
        sizem = sizem + sizeof(dummyTracer)*DBlock(i)%LTracer%getSize() !this accounts for the contents
    enddo
    call SimMemory%setracer(sizem)
    if(present(ntrc)) then
        call SimMemory%setNtrc(ntrc)
        call SimMemory%setsizeTrc(sizeof(dummyTracer))
    end if
    end subroutine setTracerMemory

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to do domain decomposition and define the Blocks
    !---------------------------------------------------------------------------
    subroutine DecomposeDomain(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext
    if (Globals%SimDefs%autoblocksize) then
        call allocBlocks(Globals%SimDefs%numblocks)
    else
        outext='[DecomposeDomain]: Only automatic Block sizing at the moment, stoping'
        call Log%put(outext)
        stop
    end if
    ! Initializing the Blocks
    call setBlocks(Globals%SimDefs%autoblocksize,Globals%SimDefs%numblocks,Globals%SimDefs%numblocksx,Globals%SimDefs%numblocksy)
    end subroutine DecomposeDomain

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation finishing method. Closes output files and writes the final messages
    !---------------------------------------------------------------------------
    subroutine closeSimulation(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext
    outext='Simulation ended, freeing resources. See you next time'
    call Log%put(outext)
    call Log%finalize()
    end subroutine closeSimulation


    end module simulation_mod
