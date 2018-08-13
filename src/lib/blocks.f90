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
    use tracer_array_mod
    use sources_array_mod
    use sources_mod
    use tracers_mod
    use emitter_mod
    use AoT_mod

    implicit none
    private

    type block_class
        integer :: id
        type(box) :: extents            !< shape::box that defines the extents of this block
        type(sourcearray_class) :: Source     !< List of Sources currently on this block
        type(tracerarray_class) :: Tracer     !< List of Tracers currently on this block
        type(aot_class) :: AoT
        type(emitter_class) :: Emitter  !< Block Emitter
        real(prec) :: resize_factor = 1.20 !< factor to resize the Tracer array once needed
    contains
    private
    procedure, public :: initialize => initBlock
    procedure, public :: putSource
    procedure, public :: CallEmitter
    procedure, public :: DistributeTracers
    procedure, public :: numActiveTracers
    procedure, public :: numAllocTracers
    procedure, public :: ToogleBlockSources
    procedure, public :: ConsolidateArrays
    procedure, public :: TracersToAoT
    procedure, public :: CleanAoT
    procedure, private :: removeTracer
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
    numAllocTracers = self%Tracer%getLength()
    end function numAllocTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that returns the total active Tracers in the Block
    !---------------------------------------------------------------------------
    function numActiveTracers(self)
    implicit none
    class(block_class), intent(in) :: self
    integer :: numActiveTracers
    numActiveTracers = self%Tracer%numActive
    end function numActiveTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to allocate and initialize blocks and their emitters
    !> @param[in] self, templatebox
    !---------------------------------------------------------------------------
    subroutine initBlock(self, id, templatebox)
    implicit none
    class(block_class), intent(inout) :: self
    integer, intent(in) :: id
    type(box), intent(in) :: templatebox
    integer :: sizem
    self%id = id
    !setting the block sub-domain
    self%extents%pt = templatebox%pt
    self%extents%size = templatebox%size
    !initializing the block emitter
    call self%Emitter%initialize(self%id)
    !initializing the Sources and Tracers arrays
    call self%Source%init(1)   !Starting the Sources array with one position
    self%Source%usedLength = 0 !But there are no stored Sources
    call self%Tracer%init(1, initvalue = dummyTracer)   !Starting the Tracers array with one position
    self%Tracer%lastActive = 0
    self%Tracer%numActive = 0 !But there are no stored Tracers
    !logging the ocupied space by the block
    sizem = sizeof(self)
    call SimMemory%addblock(sizem)
    end subroutine initBlock

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to place a Source on the Block sourcearray_class object.
    !> Checks for space and allocates more if needed.
    !> The array gets incremented by one unit at a time
    !> Allocates space in the Blocks Tracer array with a dummy Tracer
    !> @param[in] self, sourcetoput
    !---------------------------------------------------------------------------
    subroutine putSource(self, sourcetoput)
    implicit none
    class(block_class), intent(inout) :: self
    class(source_class), intent(inout) :: sourcetoput !< Source object to store
    !Check if the array is at capacity and needs to be resized
    if (self%Source%usedLength == self%Source%getLength()) then
        call self%Source%resize(self%Source%getLength()+1) !incrementing one entry
    end if
    self%Source%usedLength = self%Source%usedLength + 1
    call self%Source%put(self%Source%usedLength, sourcetoput)
    !adding this Source to the Block Emitter pool
    call self%Emitter%addSource(sourcetoput)
    !Resizing the Tracer array for the maximum possible emmited Tracers by the Sources in this Block (+1)
    call self%Tracer%resize(self%Tracer%getLength() + sourcetoput%stencil%total_np, initvalue = dummyTracer)
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

    do i=1, self%Source%usedLength
        aSource => self%Source%get(i)
        select type(aSource)
        class is (source_class)
            if (Globals%SimTime <= aSource%par%stoptime) then !SimTime smaller than Source end time
                if (Globals%SimTime >= aSource%par%startime) then !SimTime larger than source start time
                    aSource%now%active = .true.
                end if
            else !SimTime larger than Source end time
                aSource%now%active = .false.
            end if
            class default
            outext = '[Block::ToogleBlockSources] Unexepected type of content, not a Source'
            call Log%put(outext)
            stop
        end select
    end do

    end subroutine ToogleBlockSources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to activate and deactivate the sources on this block, based on
    !> Global%SimTime
    !---------------------------------------------------------------------------
    subroutine CallEmitter(self)
    implicit none
    class(block_class), intent(inout) :: self
    integer :: i
    class(*), pointer :: aSource
    type(string) :: outext

    do i=1, self%Source%usedLength
        aSource => self%Source%get(i)
        select type(aSource)
        class is (source_class)
            if (aSource%now%active) then
                aSource%now%emission_stride = aSource%now%emission_stride - 1 !decreasing the stride at this dt
                if (aSource%now%emission_stride == 0) then !reached the bottom of the stride stack, time to emitt
                    !print*, 'emitting from Block ', self%id, ' Source ', aSource%par%id
                    call self%Emitter%emitt(aSource, self%Tracer)
                    aSource%now%emission_stride = aSource%par%emitting_rate !reseting the stride after the Source emitts
                end if
            end if
            class default
            outext = '[Block::CallEmitter]: Unexepected type of content, not a Source'
            call Log%put(outext)
            stop
        end select
    end do
    end subroutine CallEmitter

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to distribute the Tracers to their correct Blocks
    !---------------------------------------------------------------------------
    subroutine DistributeTracers(self)
    implicit none
    class(block_class), intent(inout) :: self
    integer :: i, ix, iy, blk
    class(*), pointer :: aTracer
    type(string) :: outext

    do i=1, self%Tracer%lastActive
        aTracer => self%Tracer%get(i)
        select type(aTracer)
        class is (tracer_class)
            if (aTracer%now%active) then
                blk = getBlockIndex(aTracer%now%pos)
                if (blk /= self%id) then !tracer is on a different block than the current one
                    !PARALLEL this is a CRITICAL section, need to ensure correct tracer index attribution
                    call sendTracer(blk,aTracer)
                    call self%removeTracer(i)
                end if
            end if
            class default
            outext = '[Block::DistributeTracers]: Unexepected type of content, not a Tracer'
            call Log%put(outext)
            stop
        end select
    end do
    end subroutine DistributeTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to remove a Tracer from the Block. Optionally doesn't decrease the
    !> active count from the Block
    !> param[in] trc,fromCopy
    !---------------------------------------------------------------------------
    subroutine removeTracer(self,trc,fromCopy)
    implicit none
    class(block_class), intent(inout) :: self
    integer, intent(in) :: trc
    logical, optional, intent(in) :: fromCopy
    call self%Tracer%put(trc,dummyTracer)
    self%Tracer%numActive = self%Tracer%numActive - 1
    if (present(fromCopy)) then
        if (fromCopy) then
            self%Tracer%numActive = self%Tracer%numActive + 1
        end if
    end if
    if (trc == self%Tracer%lastActive) then
        self%Tracer%lastActive = self%Tracer%findLastActive()
    end if
    end subroutine removeTracer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to optimize the Tracer array from the Block - checks for memory
    !> use, sorts by active, ...
    !---------------------------------------------------------------------------
    subroutine ConsolidateArrays(self)
    implicit none
    class(block_class), intent(inout) :: self
    integer :: i
    class(*), pointer :: aTracer, bTracer
    type(string) :: outext
    !sorts the array by active tracers    
    do i=1, self%Tracer%lastActive
        aTracer => self%Tracer%get(i)
        select type(aTracer)
        class is (tracer_class)
            if (aTracer%now%active .eqv. .false.) then
                !bring the last active tracer to this position
                bTracer => self%Tracer%get(self%Tracer%lastActive)
                call self%Tracer%put(i,bTracer)                
                call self%removeTracer(self%Tracer%lastActive,.true.)                
            end if
            class default
            outext = '[Block::ConsolidateArrays]: Unexepected type of content, not a Tracer'
            call Log%put(outext)
            stop
        end select
    end do
    !resizes the array if it is too big - this needs more thinking. probably need to
    !know the trend of the numActive on this block to do this correctly - TODO
    !if (1.0*self%Tracer%getLength()/self%Tracer%lastActive .gt. self%resize_factor) then
    !    print*, 'array is too big, could be trimmed a bit...'
    !end if
    end subroutine ConsolidateArrays

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to build the AoT object at this timestep for code efficiency
    !---------------------------------------------------------------------------
    subroutine TracersToAoT(self)
    implicit none
    class(block_class), intent(inout) :: self
    self%AoT = MakeAoT(self%Tracer)
    !if (self%Tracer%numActive > 0) then
    !    print*, 'From Block ', self%id
    !    call self%AoT%print()
    !end if
    end subroutine TracersToAoT
    
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
    !> Method to send a Tracer from the current Block to another Block
    !---------------------------------------------------------------------------
    subroutine sendTracer(blk,trc)
    implicit none
    integer, intent(in) :: blk
    class(*), intent(in) :: trc
    integer :: idx
    type(string) :: outext
    !PARALLEL this is a CRITICAL section, need to ensure correct tracer
    !index attribution at the new block
    if (DBlock(blk)%Tracer%lastActive == DBlock(blk)%Tracer%getLength()) then !OPTIMIZATION - this could be numActive if the array was sorted and optimized - TODO
        call DBlock(blk)%Tracer%resize(max(int(DBlock(blk)%Tracer%getLength()*DBlock(blk)%resize_factor),10))
    end if
    idx = DBlock(blk)%Tracer%lastActive + 1
    call DBlock(blk)%Tracer%put(idx,trc)
    DBlock(blk)%Tracer%lastActive = idx
    DBlock(blk)%Tracer%numActive = DBlock(blk)%Tracer%numActive + 1
    end subroutine sendTracer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the index of the Block for a given set of coordinates.
    !> @param[in] pt
    !---------------------------------------------------------------------------
    integer function getBlockIndex(pt)
    implicit none
    type(vector), intent(in) :: pt
    integer :: ix, iy, temp
    type(string) :: outext
    ix = min(int((pt%x + BBox%offset%x)/Globals%SimDefs%blocksize%x) + 1, Globals%SimDefs%numblocksx)
    iy = min(int((pt%y + BBox%offset%y)/Globals%SimDefs%blocksize%y) + 1, Globals%SimDefs%numblocksy)
    temp = 2*ix + iy -2
    if (temp > Globals%SimDefs%numblocks) then
        outext='[Blocks::getBlockIndex]: problem in getting correct Block index, stoping'
        call Log%put(outext)
        stop
    end if
    getBlockIndex = temp
    end function getBlockIndex

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
    temp_str = self%Source%usedLength
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
    temp_str = self%Source%usedLength
    outext='      and has '//temp_str//' Sources'
    call Log%put(outext,.false.)
    call self%Source%printArray()
    end subroutine printdetailBlock


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> routine to set the simulation blocks extents and call the block initializer
    !> @param[in] self
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
        ar = get_closest_twopow(ar) !aspect ratio of our bounding box
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
                tempbox%pt = BBox%pt + BBox%size%x*(i-1)/nxi*ex + BBox%size%y*(j-1)/nyi*ey - BBox%pt%z*ez
                tempbox%size = BBox%size%x/nxi*ex + BBox%size%y/nyi*ey
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
