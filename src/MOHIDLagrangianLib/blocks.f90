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
    use simulationGlobals_mod
    use boundingbox_mod
    use tracerList_mod
    use sourcesList_mod
    use sources_mod
    use tracers_mod
    use emitter_mod
    use solver_mod
    use background_mod
    use stateVector_mod

    !use simulationTestMaker_mod

    implicit none
    private

    type :: block_class
        integer :: id
        type(box) :: extents                  !< shape::box that defines the extents of this block
        type(sourceList_class) :: LSource     !< List of Sources currently on this block
        type(emitter_class)    :: Emitter     !< Block Emitter
        type(tracerList_class) :: LTracer     !< List of Tracers currently on this block
        type(solver_class)     :: Solver      !< Block Solver
        type(stateVector_class), allocatable, dimension(:) :: BlockState     !< State vector of the Tracers in the Block, per type
        type(background_class), allocatable, dimension(:) :: Background !< Solution Backgrounds for the Block
        integer, allocatable, dimension(:,:) :: trcType
    contains
    private
    procedure, public :: initialize => initBlock
    procedure, public :: putSource
    procedure, public :: CallEmitter
    procedure, public :: DistributeTracers
    procedure, public :: ToogleBlockSources
    procedure, public :: ConsolidateArrays
    procedure, public :: ShedMemory
    procedure, public :: RunSolver
    procedure, public :: numAllocTracers
    procedure, public :: TracersToSV
    procedure, public :: SVtoTracers    
    procedure, public :: CleanSV
    procedure, public :: print => printBlock
    procedure, public :: detailedprint => printdetailBlock
    end type block_class

    !Simulation variables
    type(block_class), allocatable, dimension(:) :: sBlock

    !Public access vars
    public :: sBlock, block_class
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

    !Tests
    !allocate(self%Background(1))
    !call TestMaker%initialize(3, self%extents, self%Background(1))

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
    !> Globals%SimTime%CurrTime
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
            aSource%now%active = aSource%par%activeTime(Globals%Sim%getnumdt())
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
                    call sendTracer(blk,aTracer)
                    call self%LTracer%removeCurrent() !this also advances the iterator to the next position
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
    !> tracers that escaped from the Simulation Bounding Box. Also builds the
    !> counter of types and tracers per type.
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
            if (aTracer%now%active) aTracer%now%active = TrcInBBox(aTracer%now%pos, BBox) !check that the Tracer is inside the Simulation domain
            if (aTracer%now%active .eqv. .false.) then
                call self%LTracer%removeCurrent() !this advances the iterator to the next position
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
    !> Method that calls background functions to clean data no longer needed
    !---------------------------------------------------------------------------
    subroutine ShedMemory(self)
    class(block_class), intent(inout) :: self
    integer :: i
    if (allocated(self%Background)) then
        do i=1, size(self%Background)
            call self%Background(i)%ShedMemory()
        end do
    end if
    end subroutine ShedMemory
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to build the State Vector at this timestep for actual numerical work
    !---------------------------------------------------------------------------
    subroutine TracersToSV(self)
    class(block_class), intent(inout) :: self    
    integer :: i, tType, idx
    class(*), pointer :: aTracer
    logical :: builtState, notFound
    type(string) :: outext
    integer, allocatable, dimension(:,:) :: tempTrcType
    
    if (allocated(self%trcType)) deallocate(self%trcType)
    call self%LTracer%reset()                   ! reset list iterator
    do while(self%LTracer%moreValues())         ! loop while there are values
        notFound = .true.
        aTracer => self%LTracer%currentValue()  ! get current value
        select type(aTracer)
        class is (tracer_class)
                tType = aTracer%par%ttype
                if (allocated(self%trcType)) then
                    do i =1, size(self%trcType,1)
                        if (self%trcType(i,1) == tType) then
                            self%trcType(i,2) = self%trcType(i,2) + 1
                            notFound = .false.
                            exit
                        end if
                    end do
                end if
                if (notFound) then
                    if (allocated(self%trcType)) then
                        allocate(tempTrcType(size(self%trcType,1),size(self%trcType,2)))
                        tempTrcType = self%trcType
                        deallocate(self%trcType)
                        allocate(self%trcType(size(tempTrcType,1)+1,size(tempTrcType,2)))
                        self%trcType(:size(tempTrcType,1),:) = tempTrcType
                        deallocate(tempTrcType)
                        idx = size(self%trcType,1)
                    else
                        allocate(self%trcType(1,3))
                        idx = 1
                    end if
                    self%trcType(idx,1) = tType
                    self%trcType(idx,2) = 1
                    self%trcType(idx,3) = aTracer%getNumVars()
                end if
            class default
            outext = '[Block::TracersToSV]: Unexepected type of content, not a Tracer'
            call Log%put(outext)
            stop
        end select
        call self%LTracer%next()    ! increment the list iterator
    end do
    call self%LTracer%reset()                       ! reset list iterator
    
    if (allocated(self%trcType)) then
        allocate(self%BlockState(size(self%trcType,1)))
        do i=1, size(self%BlockState)
            self%BlockState(i)%idx = 1
            self%BlockState(i)%ttype = self%trcType(i,1)
            allocate(self%BlockState(i)%state(self%trcType(i,2),self%trcType(i,3)))
            allocate(self%BlockState(i)%trc(self%trcType(i,2)))
            allocate(self%BlockState(i)%active(self%trcType(i,2)))
            allocate(self%BlockState(i)%source(self%trcType(i,2)))
            allocate(self%BlockState(i)%id(self%trcType(i,2)))
            allocate(self%BlockState(i)%landMask(self%trcType(i,2)))
            self%BlockState(i)%landMask = Globals%Mask%waterVal
            allocate(self%BlockState(i)%landIntMask(self%trcType(i,2)))            
            self%BlockState(i)%landIntMask = Globals%Mask%waterVal
            allocate(self%BlockState(i)%resolution(self%trcType(i,2)))
        end do
        call self%LTracer%reset()                   ! reset list iterator
        do while(self%LTracer%moreValues())         ! loop while there are values
            builtState = .false.
            aTracer => self%LTracer%currentValue()  ! get current value
            select type(aTracer)
            class is (tracer_class)
                tType = aTracer%par%ttype
                do i = 1, size(self%BlockState)
                    if (tType == self%BlockState(i)%ttype) then
                        self%BlockState(i)%state(self%BlockState(i)%idx,:) = aTracer%getStateArray()
                        self%BlockState(i)%source(self%BlockState(i)%idx) = aTracer%par%idsource
                        self%BlockState(i)%id(self%BlockState(i)%idx) = aTracer%par%id
                        self%BlockState(i)%active(self%BlockState(i)%idx) = aTracer%now%active
                        self%BlockState(i)%trc(self%BlockState(i)%idx)%ptr => aTracer
                        self%BlockState(i)%idx = self%BlockState(i)%idx + 1
                        if(.not.allocated(self%BlockState(i)%varName)) then
                            allocate(self%BlockState(i)%varName(aTracer%getNumVars()))
                            self%BlockState(i)%varName = aTracer%varName
                        end if
                        builtstate = .true.
                        exit
                    end if
                end do
                aTracer%now%active = .false. !only gets flaged as active if the state vector says so when copying back to the Tracers
                if (.not.builtState) then
                    outext = '[Block::TracersToSV]: Tracer did not find correspoding State Vector, stoping'
                    call Log%put(outext)
                    stop
                end if
                class default
                outext = '[Block::TracersToSV]: Unexepected type of content, not a Tracer'
                call Log%put(outext)
                stop
            end select
            call self%LTracer%next()    ! increment the list iterator
        end do
        call self%LTracer%reset()                       ! reset list iterator
    end if
    end subroutine TracersToSV

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to run the solver on the data on this Block for the current
    !> timestep. Time for some actual numerical work!
    !---------------------------------------------------------------------------
    subroutine RunSolver(self)
    implicit none
    class(block_class), intent(inout) :: self
    if (allocated(self%BlockState)) then             !There are Tracers in this Block
        if (allocated(self%Background)) then    !There are Backgrounds in this Block
            !print*, 'From Block ', self%id
            call self%Solver%runStep(self%BlockState, self%Background, Globals%SimTime%CurrTime, Globals%SimDefs%dt)
        end if
    end if
    end subroutine RunSolver

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to write the data in the SV back to the Tracer objects in the list
    !---------------------------------------------------------------------------
    subroutine SVtoTracers(self)
    class(block_class), intent(inout) :: self
    integer :: i
    if (allocated(self%BlockState)) then
        do i=1, size(self%BlockState)
            call self%BlockState(i)%toTracers()
        end do
    end if
    end subroutine SVtoTracers
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to clean out the State Vector object
    !---------------------------------------------------------------------------
    subroutine CleanSV(self)
    class(block_class), intent(inout) :: self
    integer :: i
    if (allocated(self%BlockState)) then
        do i=1, size(self%BlockState)
            call self%BlockState(i)%finalize()
        end do
        deallocate(self%BlockState)
    end if
    end subroutine CleanSV

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to send a Tracer from the current Block to another Block. Checks
    !> if Block index exists, if not, Tracer is not added to any Block Tracer 
    !> list
    !> @param[in] blk, trc
    !---------------------------------------------------------------------------
    subroutine sendTracer(blk, trc)
    implicit none
    integer, intent(in) :: blk
    class(tracer_class), intent(inout) :: trc
    !PARALLEL this is a CRITICAL section, need to ensure correct tracer
    !index attribution at the new block
    call sBlock(blk)%LTracer%add(trc)
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
    getBlockIndex = Globals%SimDefs%numblocksy*(ix-1) + iy
    end function getBlockIndex

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns true if the point is inside the requested bounding box.
    !> @param[in] trc, testbox
    !---------------------------------------------------------------------------
    logical function TrcInBBox(trc, testbox)
    implicit none
    type(vector), intent(in) :: trc
    type(boundingbox_class), intent(inout) :: testbox
    TrcInBBox = .false.
    if (trc%x >= testbox%pt%x) then
        if (trc%x <= testbox%pt%x + testbox%size%x) then
            if (trc%y >= testbox%pt%y) then
                if (trc%y <= testbox%pt%y + testbox%size%y) then
                    if (trc%z >= testbox%pt%z) then
                        if (trc%z <= testbox%pt%z + testbox%size%z) then
                            TrcInBBox = .true.
                        end if
                    end if
                end if
            end if
        end if
    end if
    end function TrcInBBox

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
        if (nblk == 1) nyi = 1
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
                call sBlock(b)%initialize(b, tempbox)
                b=b+1
            end do
        end do
        temp(1) = nxi
        temp(2) = nyi
        outext='-->Automatic domain decomposition sucessful. Domain is '//temp(1)// ' X ' //temp(2)//' Blocks'
        call Log%put(outext,.false.)
    end if
    Globals%SimDefs%blocksize = sBlock(1)%extents%size
    !do i=1, size(sBlock)
    !    call sBlock(i)%print()
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
    allocate(sBlock(nblk), stat=err)
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
