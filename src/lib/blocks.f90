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
    !> parallelization strategy, in any memory architecture.
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
        type(box) :: extents !< shape::box that defines the extents of this block

        type(SourceArray) :: BlockSource
        type(TracerArray) :: BlockTracer

        type(emitter_class) :: BlockEmitter
    contains
    private
    procedure, public :: initialize => initBlocks
    end type block_class

    !Simulation variables
    type(block_class), allocatable, dimension(:) :: DBlock

    !Public access vars
    public :: DBlock
    !Public access procedures
    public :: allocBlocks

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to allocate and initialize all the simulation blocks
    !
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine initBlocks(self)
    implicit none
    class(block_class), intent(inout) :: self
    type(string) :: outext
    integer :: sizem, err

    call set_blocks_extents(Globals%SimDefs%autoblocksize,Globals%SimDefs%numblocks)

    sizem = sizeof(self)*Globals%SimDefs%numblocks
    call SimMemory%addblock(sizem)

    return
    end subroutine initBlocks
    
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
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

    

    return
    end subroutine populate
    
    

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> routine to set the simulation blocks extents
    !
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine set_blocks_extents(auto, nblk)
    implicit none
    logical, intent(in) ::  auto
    integer, intent(in) ::  nblk
    type(string) :: outext, temp(2)
    integer :: i, j, b, nxi, nyi
    real(prec) :: ar

    if (auto) then
        ar = BBox%size%x/BBox%size%y
        ar = get_closest_twopow(ar) !aspect ratio of our bounding box        
        nyi = sqrt(nblk/ar)
        if (nyi == 0) then
            temp(1) = ar
            outext='[set_blocks_extents]: block auto sizing failed. Bouding box aspect ratio = '//temp(1)//'. Stoping'
            call Log%put(outext)
            stop
        endif
        nxi = (nblk/nyi)
        
        b=1
        do i=1, nxi
            do j=1, nyi
                DBlock(b)%extents%pt = BBox%pt + BBox%size%x*(i-1)/nxi*ex + BBox%size%y*(j-1)/nyi*ey - BBox%pt%z*ez
                DBlock(b)%extents%size = BBox%size%x/nxi*ex + BBox%size%y/nyi*ey
                b=b+1
            end do
        end do
        temp(1) = nxi
        temp(2) = nyi
        outext='-->Automatic domain decomposition finished. Domain is '//temp(1)// ' X ' //temp(2)//' Blocks'
        call Log%put(outext,.false.)
    end if

    !sizem = sizeof(Tracer)
    !call SimMemory%addtracer(sizem)

    return
    end subroutine set_blocks_extents

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
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
