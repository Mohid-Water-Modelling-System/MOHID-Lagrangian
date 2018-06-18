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

    use commom_modules
    use initialize_mod
    use boundingbox_mod
    use emitter_mod
    use sources_mod
    use blocks_mod
    use about_mod

    !use simulation_objects_mod

    implicit none
    private

    type :: simulation_class   !< Parameters class

    contains
    procedure :: initialize => initSimulation
    procedure :: finalize   => closeSimulation
    procedure :: decompose  => DecomposeDomain
    procedure :: DistributeSources
    procedure :: run
    end type

    !Simulation variables
    public :: simulation_class
    !Public access vars

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !
    !> @brief
    !> Simulation run method. Runs the initialized case main time cycle.
    !---------------------------------------------------------------------------
    subroutine run(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext

    !main time cycle
    do while (Globals%SimTime .LT. Globals%Parameters%TimeMax)

        !Do your Lagrangian things here :D

        Globals%SimTime = Globals%SimTime + Globals%SimDefs%dt
    enddo

    end subroutine run

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Simulation initialization method. Effectively builds and populates the
    !> simulation objects that will be used latter on.
    !
    !> @param[in] casefilename, outpath
    !---------------------------------------------------------------------------
    subroutine initSimulation(self, casefilename, outpath)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string), intent(in) :: casefilename         !< case file name
    type(string), intent(in) :: outpath              !< Output path
    type(string) :: outext

    ! Initialize logger
    call Log%initialize(outpath)
    !Print licences and build info
    call PrintLicPreamble

    !setting every global variable and input parameter to their default
    call Globals%initialize()
    !initializing memory log
    call SimMemory%initialize()
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
    !call BBox%print()

    call self%decompose()
    
    call self%DistributeSources()

    !printing memory occupation at the time
    call SimMemory%detailedprint()

    end subroutine initSimulation

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Simulation to distribute the Sources to the blocks
    !---------------------------------------------------------------------------
    subroutine DistributeSources(self)
    implicit none
    class(simulation_class), intent(inout) :: self
    type(string) :: outext
    integer :: i, ix, iy
    real(prec) :: dx, dy
    
    !this is easy because all the blocks are the same
    dx = DBlock(1)%extents%size%x
    dy = DBlock(1)%extents%size%y
    !Now we find the 2D coordinates
    do i=1, size(tempSources%src)
        ix = int(abs(tempSources%src(i)%now%pos%x/dx)) + 1
        iy = int(abs(tempSources%src(i)%now%pos%y/dy)) + 1
        print*, 'Source position'
        print*, tempSources%src(i)%now%pos
        print*, 'Source grid position'
        print*, ix, iy
    end do
    
    end subroutine DistributeSources


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
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
    ! Initializing the blocks
    call setBlocks(Globals%SimDefs%autoblocksize,Globals%SimDefs%numblocks)

    end subroutine DecomposeDomain

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
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
