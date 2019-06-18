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

    use simulationAbout_mod
    use simulationInitialize_mod
    use boundingbox_mod
    use blocks_mod
    use emitter_mod
    use sources_mod
    use tracers_mod
    use background_mod
    use simulationOutputStreamer_mod
    use simulationInputStreamer_mod
    use common_modules

    use hdf5Writter_mod

    implicit none
    private

    type :: simulation_class   !< Parameters class
        !Timers
        type(timer_class) :: timerInit          !< timer for the initialization routines
        type(timer_class) :: timerTotalRun      !< timer for the total wall time spent on the simulation
        type(timer_class) :: timerInput        !< timer for the output writting routines
        type(timer_class) :: timerOutput        !< timer for the output writting routines
        type(timer_class) :: timerPrep          !< timer for the toggling, emission and consolidation phase of every time-step
        type(timer_class) :: timerAoTOps        !< timer for the AoT operations (encoding and decoding)
        type(timer_class) :: timerSolver        !< timer for the solver runs
        !Output objects
        type(output_streamer_class) :: OutputStreamer !< Writter that streams simulation data out
        !Input objects
        type(input_streamer_class)  :: InputStreamer  !< Reader that streams data into the simulation
    contains
    procedure, public  :: initialize => initSimulation
    procedure, public  :: run
    procedure, public  :: finalize   => closeSimulation
    procedure, private :: decompose  => DecomposeDomain
    procedure, private :: ToggleSources
    procedure, private :: BlocksEmitt
    procedure, private :: BlocksDistribute
    procedure, private :: BlocksConsolidate
    procedure, private :: BlocksTracersToAoT
    procedure, private :: BlocksRunSolver
    procedure, private :: BlocksAoTtoTracers
    procedure, private :: BlocksCleanAoT
    procedure, private :: InputData
    procedure, private :: OutputStepData
    procedure, private :: setInitialState
    procedure, private :: updateSimDateTime
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
    class(simulation_class), intent(inout) :: self

    call self%OutputStreamer%writeOutputHeader()

    !main time cycle
    do while (Globals%SimTime%CurrTime .lt. Globals%Parameters%TimeMax)
        call Globals%Sim%increment_numdt()
        call self%timerTotalRun%Tic()
        !activate suitable Sources
        call self%ToggleSources()
        !emitt Tracers from active Sources
        call self%BlocksEmitt()
        !Optimize Block Tracer lists and memory use
        call self%BlocksConsolidate()
        !Distribute Tracers and Sources by Blocks
        call self%BlocksDistribute()
        !Build AoT
        call self%BlocksTracersToAoT()
        !load hydrodynamic fields from files (curents, wind, waves, ...)
        call self%InputData()
        !Update all tracers with base behavior (AoT) - Integration step
        call self%BlocksRunSolver()
        !AoT to Tracers
        call self%BlocksAoTtoTracers()
        !Update Tracers with type-specific behavior

        !Write results if time to do so
        call self%OutputStepData()
        !Clean AoT
        call self%BlocksCleanAoT()
        !update Simulation time
        call self%updateSimDateTime()
        !print*, 'Global time is ', Globals%SimTime%CurrTime
        !print*, 'Can we continue?'
        !read (*,*)
        call self%timerTotalRun%Toc()
    enddo
    
    call self%setTracerMemory()
    call SimMemory%detailedprint()
    call self%timerTotalRun%print()
    call self%timerPrep%print()
    call self%timerAoTOps%print()
    call self%timerSolver%print()
    call self%timerInput%print()
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
    class(simulation_class), intent(inout) :: self
    type(string), intent(in) :: casefilename         !< case file name
    type(string), intent(in) :: outpath              !< Output path
    type(string) :: outext, aux

    aux = 'Simulation::initialization'
    call self%timerInit%initialize(aux)
    call self%timerInit%Tic()

    aux = 'Simulation::Total'
    call self%timerTotalRun%initialize(aux)
    aux = 'Simulation::Input'
    call self%timerInput%initialize(aux)
    aux = 'Simulation::Output'
    call self%timerOutput%initialize(aux)
    aux = 'Simulation::Preparation'
    call self%timerPrep%initialize(aux)
    aux = 'Simulation::AoT encoding/decoding'
    call self%timerAoTOps%initialize(aux)
    aux = 'Simulation::Solver'
    call self%timerSolver%initialize(aux)

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
    !Initializing input file streamer
    call self%InputStreamer%initialize()
    call self%InputStreamer%print()
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
    call self%OutputStreamer%initialize()
    !Writing the domain to file
    call self%OutputStreamer%WriteDomain(Globals%Names%casename, BBox, Geometry%getnumPoints(BBox), sBlock)

    call self%timerInit%Toc()
    call self%timerInit%print()

    end subroutine initSimulation

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to activate and deactivate Sources based on the
    !> Globals%SimTime%CurrTime
    !---------------------------------------------------------------------------
    subroutine ToggleSources(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerPrep%Tic()
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(sBlock)
        call sBlock(i)%ToogleBlockSources()
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    call self%timerPrep%Toc()
    end subroutine ToggleSources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to emitt tracers at current Time
    !---------------------------------------------------------------------------
    subroutine BlocksEmitt(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerPrep%Tic()
    do i=1, size(sBlock)
        call sBlock(i)%CallEmitter()
    enddo
    call self%timerPrep%Toc()
    end subroutine BlocksEmitt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to distribute Tracers at
    !> current Time
    !---------------------------------------------------------------------------
    subroutine BlocksDistribute(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerPrep%Tic()
    do i=1, size(sBlock)
        call sBlock(i)%DistributeTracers()
    enddo
    call self%timerPrep%Toc()
    !need to distribute Sources also! TODO
    end subroutine BlocksDistribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to consolidate the Tracer array at
    !> current Time and optimize background memory useage
    !---------------------------------------------------------------------------
    subroutine BlocksConsolidate(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerPrep%Tic()
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(sBlock)
        call sBlock(i)%ConsolidateArrays()
        call sBlock(i)%ShedMemory()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call self%timerPrep%Toc()
    end subroutine BlocksConsolidate

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to build their Array of
    !> Tracers (AoT) from the Tracer list at current Time
    !---------------------------------------------------------------------------
    subroutine BlocksTracersToAoT(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerAoTOps%Tic()
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(sBlock)
        call sBlock(i)%TracersToAoT()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call self%timerAoTOps%Toc()
    end subroutine BlocksTracersToAoT

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to run the Solver at
    !> current Time
    !---------------------------------------------------------------------------
    subroutine BlocksRunSolver(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    if (Globals%Sim%getnumdt() /= 1 ) then
        call self%timerSolver%Tic()
        !$OMP PARALLEL PRIVATE(i)
        !$OMP DO
        do i=1, size(sBlock)
            call sBlock(i)%RunSolver()
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        call self%timerSolver%Toc()
    end if
    end subroutine BlocksRunSolver

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to print their Array of
    !> Tracers (AoT) back to the Tracer objects on the list at current Time
    !---------------------------------------------------------------------------
    subroutine BlocksAoTtoTracers(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerAoTOps%Toc()
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(sBlock)
        call sBlock(i)%AoTtoTracers()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call self%timerAoTOps%Toc()
    end subroutine BlocksAoTtoTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the Blocks to clean their Array of
    !> Tracers (AoT) at current Time
    !---------------------------------------------------------------------------
    subroutine BlocksCleanAoT(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    call self%timerAoTOps%Tic()
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1, size(sBlock)
        call sBlock(i)%CleanAoT()
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call self%timerAoTOps%Toc()
    end subroutine BlocksCleanAoT
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to update the date and time variables
    !---------------------------------------------------------------------------
    subroutine updateSimDateTime(self)
    class(simulation_class), intent(inout) :: self
    integer :: i
    if (Globals%Sim%getnumdt() /= 1 ) call Globals%SimTime%setCurrDateTime(Globals%SimDefs%dt)
    end subroutine updateSimDateTime
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the input streamer reading methods at
    !> current Time
    !---------------------------------------------------------------------------
    subroutine InputData(self)
    class(simulation_class), intent(inout) :: self
    call self%timerInput%Tic()
    call self%InputStreamer%loadDataFromStack(BBox, sBlock)
    call self%timerInput%Toc()
    end subroutine InputData

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to call the output streamer writting methods at
    !> current Time
    !---------------------------------------------------------------------------
    subroutine OutputStepData(self)
    class(simulation_class), intent(inout) :: self
    call self%timerOutput%Tic()
    call self%OutputStreamer%WriteStep(sBlock, self%getTracerTotals(), self%timerTotalRun)
    call self%timerOutput%Toc()
    end subroutine OutputStepData

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to distribute the Sources to the Blocks, allocate the
    !> respective Tracers and redistribute if needed
    !---------------------------------------------------------------------------
    subroutine setInitialState(self)
    class(simulation_class), intent(inout) :: self
    type(string) :: outext
    integer :: i, blk, ntrc
    !iterate every Source to distribute
    ntrc = 0
    do i=1, size(tempSources%src)
        blk = getBlockIndex(Geometry%getCenter(tempSources%src(i)%par%geometry))
        call sBlock(blk)%putSource(tempSources%src(i))
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
    class(simulation_class), intent(in) :: self
    integer :: i, total
    total = 0
    do i=1, size(sBlock)
        total = total + sBlock(i)%numAllocTracers()
    enddo
    getTracerTotals = total
    end function getTracerTotals

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to count Tracer numbers
    !---------------------------------------------------------------------------
    subroutine printTracerTotals(self)
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
    class(simulation_class), intent(in) :: self
    integer, optional, intent(in) :: ntrc
    integer :: sizem, i
    type(tracer_class) :: dummyTracer
    sizem = 0
    do i=1, size(sBlock)
        sizem = sizem + sizeof(sBlock(i)%LTracer) !this accounts for the array structure
        sizem = sizem + sizeof(dummyTracer)*sBlock(i)%LTracer%getSize() !this accounts for the contents
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
    class(simulation_class), intent(inout) :: self
    type(string) :: outext
    outext='Simulation ended, freeing resources. See you next time'
    call Log%put(outext)
    call Log%finalize()
    call self%OutputStreamer%finalize()
    end subroutine closeSimulation


    end module simulation_mod