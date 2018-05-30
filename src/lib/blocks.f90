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
    procedure :: initialize => initBlocks

  end type block_class

  !Simulation variables
  type(block_class), allocatable, dimension(:) :: Block

  !Public access vars
  public :: Block
  !Public access procedures

contains

  !---------------------------------------------------------------------------
  !> @Ricardo Birjukovs Canelas - MARETEC
  ! Routine Author Name and Affiliation.
  !
  !> @brief
  !> method to allocate and initialize the simulation blocks
  !
  !> @param[in] self
  !---------------------------------------------------------------------------
  subroutine initBlocks(self)
    implicit none
    class(block_class), intent(inout) :: self
    type(string) :: outext
    integer :: sizem, err

    if (Globals%SimDefs%autoblocksize) then
      call allocBlocks(Globals%SimDefs%numblocks)
    else
      outext='[initblocks]: Only automatic block sizing at the moment, stoping'
      call Log%put(outext)
      stop
    end if
    call set_blocks_extents(Globals%SimDefs%autoblocksize,Globals%SimDefs%numblocks)

    !sizem = sizeof(Tracer)
    !call SimMemory%addtracer(sizem)

    return
  end subroutine initBlocks

  !---------------------------------------------------------------------------
  !> @Ricardo Birjukovs Canelas - MARETEC
  ! Routine Author Name and Affiliation.
  !
  !> @brief
  !> routine to set the simulation blocks extents
  !
  !> @param[in] self
  !---------------------------------------------------------------------------
  subroutine allocBlocks(nblk)
    implicit none
    integer, intent(in) ::  nblk
    type(string) :: outext
    integer err
    allocate(Block(nblk), stat=err)
    if(err/=0)then
      outext='[allocBlobks]: Cannot allocate Blocks, stoping'
      call Log%put(outext)
      stop
    endif
  end subroutine allocBlocks

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
    type(string) :: outext
    integer :: i, j, b, nxi, nyi
    real(prec) :: ar, nx, ny

    if (auto) then
      ar = BBox%size%x/BBox%size%y !aspect ratio of our bounding box
      ny = floor(sqrt(nblk/ar))
      nxi = ceiling(nblk/ny)
      nyi = nblk/nxi

      if (nxi*nyi == nblk) then
        b=1
        do i=1, nxi
          do j=1, nyi
            Block(b)%extents%pt = BBox%pt + BBox%size%x*(i-1)*ex + BBox%size%y*(j-1)*ey
            b=b+1
          end do
        end do
      else
        outext='[set_blocks_extents]: block auto sizing failed. Check BB aspect ratio. Stoping'
        call Log%put(outext)
        stop
      end if
    end if

    !sizem = sizeof(Tracer)
    !call SimMemory%addtracer(sizem)

    return
  end subroutine set_blocks_extents

end module blocks_mod
