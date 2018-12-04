    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : blocks
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a block class and related methods. A block is a fundamental
    !> type of the model. It contains a sub-domain of the simulation bounding box,
    !> holding all entities inside that sub-domain. It maps to a domain decomposition
    !> parallelization strategy, if needed.
    !------------------------------------------------------------------------------

    module blocks_mod

    use common_modules
    use simulation_globals_mod
    use boundingbox_mod
    use tracer_list_mod
    use sources_list_mod
    use sources_mod
    use tracers_mod
    use emitter_mod
    use AoT_mod
    use solver_mod
    use background_mod

    use simulation_testmaker_mod


    implicit none
    private

    type :: block_class
        integer :: id
        type(box) :: extents                  !< shape::box that defines the extents of this block
        type(sourceList_class) :: LSource     !< List of Sources currently on this block
        type(emitter_class)    :: Emitter     !< Block Emitter
        type(tracerList_class) :: LTracer     !< List of Tracers currently on this block
        type(aot_class)        :: AoT         !< Block Array of Tracers for actual numerical work
        type(solver_class)     :: Solver      !< Block Solver
        type(background_class), allocatable, dimension(:) :: Background !< Solution Backgrounds for the Block
    contains
    private
    procedure, public :: initialize => initBlock
    procedure, public :: putSource
    procedure, public :: CallEmitter
    procedure, public :: DistributeTracers
    procedure, public :: ToogleBlockSources
    procedure, public :: ConsolidateArrays
    procedure, public :: TracersToAoT
    procedure, public :: RunSolver
    procedure, public :: AoTtoTracers
    procedure, public :: CleanAoT
    procedure, public :: numAllocTracers
    procedure, public :: print => printBlock
    procedure, public :: detailedprint => printdetailBlock
    end type block_class

    !Simulation variables
    type(block_class), allocatable, dimension(:) :: DBlock

    !Public access vars
    public :: DBlock, block_class
    !Public access procedures
    public :: allocBlocks, setBlocks, getBlockIndex

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that returns the total allocated Tracers in the Block
    !---------------------------------------------------------------------------
    function numAllocTracers(self)
    implicit none
    class(block_class), intent(in) :: self
    integer :: numAllocTracers
    numAllocTracers = self%LTracer%getSize()
    end function numAllocTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to allocate and initialize Blocks and their Emitters
    !> @param[in] self, id, templatebox
    !---------------------------------------------------------------------------
    subroutine initBlock(self, id, templatebox)
    implicit none
    class(block_class), intent(inout) :: self
    integer, intent(in) :: id
    type(box), intent(in) :: templatebox
    integer :: sizem, i
    self%id = id
    !setting the block sub-domain
    self%extents = templatebox
    !initializing the block emitter
    call self%Emitter%initialize()
    !initializing the block solver
    i = Globals%Parameters%Integrator
    call self%Solver%initialize(i, Globals%Parameters%IntegratorNames(i))
    sizem = sizeof(self)
    call SimMemory%addblock(sizem)

    allocate(self%Background(1))
    call TestMaker%initialize(2, self%extents, self%Background(1))
    !call self%print()
    !call self%Background(1)%print()

    end subroutine initBlock

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to place a Source on the Block sourceList_class object. Adds the
    !> Source info to the Block Emitter
    !> @param[in] self, sourcetoadd
    !---------------------------------------------------------------------------
    subroutine putSource(self, sourcetoadd)
    implicit none
    class(block_class), intent(inout) :: self
    class(source_class), intent(inout) :: sourcetoadd !< Source object to store
    call self%LSource%add(sourcetoadd)
    end subroutine putSource

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to activate and deactivate the sources on this block, based on
    !> Globa%SimTime
    !---------------------------------------------------------------------------
    subroutine ToogleBlockSources(self)
    implicit none
    class(block_class), intent(inout) :: self
    integer :: i
    class(*), pointer :: aSource
    type(string) :: outext

    call self%LSource%reset()                   ! reset list iterator
    do while(self%LSource%moreValues())         ! loop while there are values
        aSource => self%LSource%currentValue()  ! get current value
        select type(aSource)
        class is (source_class)
            if (Globals%SimTime <= aSource%par%stoptime) then       !SimTime smaller than Source end time
                if (Globals%SimTime >= aSource%par%startime) then   !SimTime larger than source start time
                    aSource%now%active = .true.
                end if
            else            !SimTime larger than Source end time
                aSource%now%active = .false.
            end if
            class default
            outext = '[Block::ToogleBlockSources] Unexepected type of content, not a Source'
            call Log%put(outext)
            stop
        end select
        call self%LSource%next()            ! increment the list iterator
    end do
    call self%LSource%reset()               ! reset list iterator

    end subroutine ToogleBlockSources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to emitt Tracers from currently active Sources on the Block
    !---------------------------------------------------------------------------
    subroutine CallEmitter(self)
    implicit none
    class(block_class), intent(inout) :: self
    call self%Emitter%emitt(self%LSource, self%LTracer)
    end subroutine CallEmitter

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to distribute the Tracers to their correct Blocks
    !---------------------------------------------------------------------------
    subroutine DistributeTracers(self)
    implicit none
    class(block_class), intent(inout) :: self
    integer :: i, blk
    class(*), pointer :: aTracer
    type(string) :: outext
    logical :: notremoved

    call self%LTracer%reset()                   ! reset list iterator
    do while(self%LTracer%moreValues())         ! loop while there are values
        notremoved = .true.
        aTracer => self%LTracer%currentValue()  ! get current value
        select type(aTracer)
        class is (tracer_class)
            if (aTracer%now%active) then
                blk = getBlockIndex(aTracer%now%pos)
                if (blk /= self%id) then        !tracer is on a different block than the current one
                    !PARALLEL this is a CRITICAL section, need to ensure correct tracer index attribution
                    print*, 'Trc ', aTracer%par%id, 'changing block from ', self%id, 'to ', blk
                    call sendTracer(blk,aTracer)
                    call self%LTracer%removeCurrent() !this also advances the iterator to the next position
                    !call self%LTracer%lowerNumActive()
                    notremoved = .false.
                end if
            end if
            class default
            outext = '[Block::DistributeTracers]: Unexepected type of content, not a Tracer'
            call Log%put(outext)
            stop
        end select
        if (notremoved) call self%LTracer%next()    ! increment the list iterator
    end do
    call self%LTracer%reset()                   ! reset list iterator

    end subroutine DistributeTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to clean the Tracer list from inactive Tracers. This includes
    !> tracers that escaped from the Simulation Bounding Box
    !---------------------------------------------------------------------------
    subroutine ConsolidateArrays(self)
    implicit none
    class(block_class), intent(inout) :: self
    class(*), pointer :: aTracer
    type(string) :: outext
    logical :: notremoved

    call self%LTracer%reset()                   ! reset list iterator
    do while(self%LTracer%moreValues())         ! loop while there are values
        notremoved = .true.
        aTracer => self%LTracer%currentValue()  ! get current value
        select type(aTracer)
        class is (tracer_class)
            if (aTracer%now%active) aTracer%now%active = TrcInBox(aTracer%now%pos, BBox) !check that the Tracer is inside the Simulation domain
            if (aTracer%par%id == MV) print*, 'found a problem '
            if (aTracer%now%active .eqv. .false.) then
                call self%LTracer%removeCurrent() !this advances the iterator to the next position
                !call self%LTracer%lowerNumActive()
                notremoved = .false.
            end if
            class default
            outext = '[Block::ConsolidateArrays]: Unexepected type of content, not a Tracer'
            call Log%put(outext)
            stop
        end select
        if (notremoved) call self%LTracer%next()    ! increment the list iterator
    end do
    call self%LTracer%reset()                       ! reset list iterator

    end subroutine ConsolidateArrays

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to build the AoT object at this timestep for actual numerical work
    !---------------------------------------------------------------------------
    subroutine TracersToAoT(self)
    implicit none
    class(block_class), intent(inout) :: self
    print*, '----printing List'
    call self%LTracer%print()
    self%AoT = AoT(self%LTracer)
    if (self%LTracer%getSize() > 0) then
        print*, 'From Block ', self%id
        call self%AoT%print()
    end if
    end subroutine TracersToAoT

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to run the solver on the data on this Block for the current
    !> timestep. Time for some actual numerical work!
    !---------------------------------------------------------------------------
    subroutine RunSolver(self)
    implicit none
    class(block_class), intent(inout) :: self
    if (size(self%AoT%id) > 0) then             !There are Tracers in this Block
        if (allocated(self%Background)) then    !There are Backgrounds in this Block
            call self%Solver%runStep(self%AoT, self%Background, Globals%SimTime, Globals%SimDefs%dt)
        end if
    end if
    end subroutine RunSolver

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to write the data in the AoT back to the Tracer objects in the list
    !---------------------------------------------------------------------------
    subroutine AoTtoTracers(self)
    implicit none
    class(block_class), intent(inout) :: self
    !call self%AoT%print()
    call self%AoT%toTracers()
    end subroutine AoTtoTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to clean out the AoT object
    !---------------------------------------------------------------------------
    subroutine CleanAoT(self)
    implicit none
    class(block_class), intent(inout) :: self
    call self%AoT%Clean()
    end subroutine CleanAoT

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to send a Tracer from the current Block to another Block. Checks
    !> if Block index exists, if not, Tracer is not added to any Block Tracer 
    !> list
    !---------------------------------------------------------------------------
    subroutine sendTracer(blk,trc)
    implicit none
    integer, intent(in) :: blk
    class(tracer_class), intent(inout) :: trc
    !PARALLEL this is a CRITICAL section, need to ensure correct tracer
    !index attribution at the new block
    print*, 'trc to send'
    call trc%print()
    print*, 'block list to add to'
    call DBlock(blk)%LTracer%print()
    call DBlock(blk)%LTracer%add(trc)
    ! if (blk <= size(DBlock)) then
    !     if (blk > 0) call DBlock(blk)%LTracer%add(trc)
    ! end if
    print*, 'block list to added to'
    call DBlock(blk)%LTracer%print()

    end subroutine sendTracer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the index of a Block for a given set of coordinates.
    !> @param[in] pt
    !---------------------------------------------------------------------------
    integer function getBlockIndex(pt)
    implicit none
    type(vector), intent(in) :: pt
    integer :: ix, iy
    ix = min(int((pt%x + BBox%offset%x)/Globals%SimDefs%blocksize%x) + 1, Globals%SimDefs%numblocksx)
    iy = min(int((pt%y + BBox%offset%y)/Globals%SimDefs%blocksize%y) + 1, Globals%SimDefs%numblocksy)
    getBlockIndex = 2*ix + iy -2
    end function getBlockIndex

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns true if the point is inside the requested box.
    !> @param[in] pt, testbox
    !---------------------------------------------------------------------------
    logical function TrcInBox(trc, testbox)
    implicit none
    type(vector), intent(in) :: trc
    type(boundingbox_class), intent(inout) :: testbox
    TrcInBox = .false.
    if (trc%x >= testbox%pt%x) then
        if (trc%x <= testbox%pt%x + testbox%size%x) then
            if (trc%y >= testbox%pt%y) then
                if (trc%y <= testbox%pt%y + testbox%size%y) then
                    if (trc%z >= testbox%pt%z) then
                        if (trc%z <= testbox%pt%z + testbox%size%z) then
                            TrcInBox = .true.
                        end if
                    end if
                end if
            end if
        end if
    end if
    end function TrcInBox

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print basic info about the block
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printBlock(self)
    implicit none
    class(block_class), intent(inout) :: self
    type(string) :: outext, temp_str
    temp_str = self%id
    outext='-->Block '//temp_str//' is a'
    call Log%put(outext,.false.)
    call Geometry%print(self%extents)
    temp_str = self%LSource%getSize()
    outext='      and has '//temp_str//' Sources'
    call Log%put(outext,.false.)
    end subroutine printBlock

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print detailed info about the block
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printdetailBlock(self)
    implicit none
    class(block_class), intent(inout) :: self
    type(string) :: outext, temp_str
    integer :: i
    temp_str = self%id
    outext='-->Block '//temp_str//' is a'
    call Log%put(outext,.false.)
    call Geometry%print(self%extents)
    temp_str = self%LSource%getSize()
    outext='      and has '//temp_str//' Sources'
    call Log%put(outext,.false.)
    call self%LSource%print()
    end subroutine printdetailBlock

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> routine to set the simulation blocks extents and call the block initializer
    !> @param[in] auto, nblk, nxi, nyi
    !---------------------------------------------------------------------------
    subroutine setBlocks(auto, nblk, nxi, nyi)
    implicit none
    logical, intent(in) ::  auto
    integer, intent(in) ::  nblk
    integer, intent(out) :: nxi, nyi
    type(string) :: outext, temp(2)
    integer :: i, j, b
    real(prec) :: ar
    type(box) :: tempbox

    if (auto) then
        ar = BBox%size%x/BBox%size%y
        ar = Utils%get_closest_twopow(ar) !aspect ratio of our bounding box
        nyi = sqrt(nblk/ar)
        if (nyi == 0) then
            temp(1) = ar
            outext='[setBlocks]: block auto sizing failed. Bouding box aspect ratio = '//temp(1)//'. Stoping'
            call Log%put(outext)
            stop
        endif
        nxi = (nblk/nyi)

        b=1
        do i=1, nxi
            do j=1, nyi
                tempbox%pt = BBox%pt + BBox%size%x*(i-1)/nxi*ex + BBox%size%y*(j-1)/nyi*ey
                tempbox%size = BBox%size%x/nxi*ex + BBox%size%y/nyi*ey + BBox%size%z*ez
                call DBlock(b)%initialize(b, tempbox)
                b=b+1
            end do
        end do
        temp(1) = nxi
        temp(2) = nyi
        outext='-->Automatic domain decomposition sucessful. Domain is '//temp(1)// ' X ' //temp(2)//' Blocks'
        call Log%put(outext,.false.)
    end if
    Globals%SimDefs%blocksize = DBlock(1)%extents%size
    !do i=1, size(DBlock)
    !    call DBlock(i)%print()
    !enddo

    return
    end subroutine setBlocks

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> routine to allocate the simulation blocks
    !> @param[in] nblk
    !---------------------------------------------------------------------------
    subroutine allocBlocks(nblk)
    implicit none
    integer, intent(in) ::  nblk
    type(string) :: outext, temp
    integer err
    allocate(DBlock(nblk), stat=err)
    if(err/=0)then
        outext='[allocBlobks]: Cannot allocate Blocks, stoping'
        call Log%put(outext)
        stop
    else
        temp = nblk
        outext = 'Allocated '// temp // ' Blocks.'
        call Log%put(outext)
    endif
    end subroutine allocBlocks

    end module blocks_mod
