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

    use commom_modules
    use simulation_globals_mod
    use boundingbox_mod
    use tracer_array_mod
    use sources_array_mod
    use emitter_mod

    implicit none
    private

    type block_class
        integer :: id
        type(box) :: extents            !< shape::box that defines the extents of this block
        type(SourceArray) :: Source     !< List of Sources currently on this block
        type(TracerArray) :: Tracer     !< List of Tracers currently on this block
        type(emitter_class) :: Emitter  !< Block Emitter
    contains
    private
    procedure, public :: initialize => initBlock
    procedure, public :: print => printBlock
    end type block_class

    !Simulation variables
    type(block_class), allocatable, dimension(:) :: DBlock

    !Public access vars
    public :: DBlock, block_class
    !Public access procedures
    public :: allocBlocks, setBlocks

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to allocate and initialize blocks and their emitters
    !
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
    call self%Emitter%initialize()
    !initializing the Sources and Tracers arrays

    
    !logging the ocupied space by the block
    sizem = sizeof(self)
    call SimMemory%addblock(sizem)
    end subroutine initBlock

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Method to populate the Blocks with their initial Sources and Tracers
    !> This copies the Sources from their temporary global position and then
    !> allocates the foreseable tracers in their arrays.
    !
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine populate(self)
    implicit none
    class(block_class), intent(inout) :: self
    type(string) :: outext
    integer :: sizem, err




    end subroutine populate
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Method to print basic info about the block
    !
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
    end subroutine printBlock


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> routine to set the simulation blocks extents and call the block initializer
    !
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
    !do i=1, size(DBlock)
    !    call DBlock(i)%print()
    !enddo

    return
  end subroutine setBlocks

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> routine to allocate the simulation blocks
    !
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
