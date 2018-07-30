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
    !> Module to hold the simulation class and its methods
    !------------------------------------------------------------------------------
    module simulation_mod

    use common_modules
    use initialize_mod
    use boundingbox_mod
    use emitter_mod
    use sources_mod
    use tracers_mod
    use blocks_mod
    use about_mod
    use simulation_output_streamer_mod

    !use simulation_objects_mod

    implicit none
    private

    type :: simulation_class   !< Parameters class
        private
        integer :: nbx, nby               !< number of blocks in 2D
    contains
    procedure, public  :: initialize => initSimulation
    procedure, public  :: run
    procedure, public  :: finalize   => closeSimulation
    procedure, private :: decompose  => DecomposeDomain
    procedure, private :: ToggleSources
    procedure, private :: BlocksEmitt
    procedure, private :: setInitialState
    procedure, private :: getTracerTotals
    procedure, private :: printTracerTotals
    procedure, private :: setTracerMemory
    procedure, private :: saveStateToFile
    end type

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
    type(string) :: outext

    !main time cycle
    do while (Globals%SimTime .lt. Globals%Parameters%TimeMax)
        !activate suitable Sources
        call self%ToggleSources()
        !emitt Tracers from active Sources
        call self%BlocksEmitt()
        !Distribute Tracers by Blocks
        !Optimize Block Tracer arrays (sort and consider building AoT)
        !load hydrodynamic fields from files
        !interpolate fields to tracer coordinates
        !Update all tracers with base behavior
        !Update Tracers with type-specific behavior
        !Write results if time to do so

        call self%printTracerTotals()
        !update Simulation time and counters
        Globals%SimTime = Globals%SimTime + Globals%SimDefs%dt
        call Globals%Sim%increment_numdt()
        !print*, 'Global time is ', Globals%SimTime
        !print*, 'Can we continue?'
        !read (*,*)
    enddo

    end subroutine run

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation initialization method. Effectively builds and populates the
    !> simulation objects that will be used latter on.
    !> @param[in] casefilename, outpath
    !---------------------------------------------------------------------------
    subroutine initSimulation(self, casefilename, outpath)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string), intent(in) :: casefilename         !< case file name
    type(string), intent(in) :: outpath              !< Output path
    type(string) :: outext
    type(vector) :: tempvec

    ! Initialize logger
    call Log%initialize(outpath)
    !Print licences and build info
    call PrintLicPreamble

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

    !Distributing Sources and trigerring Tracer allocation and distribution
    call self%setInitialState()

    !printing memory occupation at the time
    call SimMemory%detailedprint()
    
    !Initializing output file streamer
    call OutputStreamer%initialize()
    
    !Writing the domain to file
    call OutputStreamer%WriteDomain(Globals%Names%casename, BBox, Geometry%getnumPoints(BBox), DBlock)

    end subroutine initSimulation

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to write the current state to file
    !---------------------------------------------------------------------------
    subroutine saveStateToFile(self)
        implicit none
        class(simulation_class), intent(in) :: self
               
    end subroutine saveStateToFile

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
        do i=1, size(DBlock)
            call DBlock(i)%ToogleBlockSources()
        enddo        
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
    !> Simulation method to distribute the Sources to the Blocks, allocate the 
    !> respective Tracers and redistribute if needed
    !---------------------------------------------------------------------------
    subroutine setInitialState(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext, temp(2)
    integer :: i, ix, iy, blk
    real(prec) :: dx, dy
    type(vector) :: coords
    
    !this is easy because all the blocks are the same
    dx = DBlock(1)%extents%size%x
    dy = DBlock(1)%extents%size%y
    !iterate every Source to distribute
    do i=1, size(tempSources%src)
        coords = Geometry%getCenter(tempSources%src(i)%par%geometry)
        !call Geometry%getCenter(tempSources%src(i)%par%geometry, coords)
        !finding the 2D coordinates of the corresponding Block
        ix = min(int((coords%x + BBox%offset%x)/dx) + 1, self%nbx)
        iy = min(int((coords%y + BBox%offset%y)/dy) + 1, self%nby)
        !print*, 'Source pt position'
        !print*, tempSources%src(i)%now%pos
        !print*, 'Source center position'
        !print*, coords
        !print*, 'Source grid position'
        !print*, ix, iy
        !Converting to the 1D index - Notice how the blocks were built in [Blocks::setBlocks]
        blk = 2*ix + iy -2
        !print*, blk
        if (blk > size(DBlock)) then
            outext='[DistributeSources]: problem in getting correct Block index, stoping'
            call Log%put(outext)
            stop
        end if
        call DBlock(blk)%putSource(tempSources%src(i))
    end do
    
    call tempSources%finalize() !destroying the temporary Sources now they are shipped to the Blocks
    outext='-->Sources allocated to their current Blocks'
    call Log%put(outext,.false.)
    
    call self%printTracerTotals()
    call self%setTracerMemory()
    
    end subroutine setInitialState
    
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to count Tracer numbers
    !---------------------------------------------------------------------------
    subroutine getTracerTotals(self, alloc, active)
    implicit none
    class(simulation_class), intent(in) :: self
    integer, intent(out) :: alloc, active
    integer :: i
    alloc = 0
    active = 0
    do i=1, size(DBlock)
        alloc = alloc + DBlock(i)%numAllocTracers()
        active = active + DBlock(i)%numActiveTracers()
    enddo        
    end subroutine getTracerTotals
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to count Tracer numbers
    !---------------------------------------------------------------------------
    subroutine printTracerTotals(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: alloc, active
    type(string) :: outext, temp(2)
    call self%getTracerTotals(alloc, active)
    temp(1) = alloc
    temp(2) = active
    outext='-->'//temp(1) //' Tracers allocated, '//temp(2) //' Tracers active'    
    call Log%put(outext,.false.)    
    end subroutine printTracerTotals
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Simulation method to account for Tracer memory consumption
    !---------------------------------------------------------------------------
    subroutine setTracerMemory(self)
    implicit none
    class(simulation_class), intent(in) :: self
    integer :: alloc, active
    integer :: sizem, i
    sizem = 0
    do i=1, size(DBlock)
        sizem = sizem + DBlock(i)%Tracer%getMemSize() !this accounts for the array structure
        sizem = sizem + sizeof(dummyTracer)*DBlock(i)%Tracer%getLength() !this accounts for the contents
    enddo  
    call SimMemory%addtracer(sizem)
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
    call setBlocks(Globals%SimDefs%autoblocksize,Globals%SimDefs%numblocks,self%nbx,self%nby)
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
