    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : boundingbox_mod
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.3
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a simulation Bounding Box.
    !------------------------------------------------------------------------------

    module boundingbox_mod

    use common_modules

    implicit none
    private

    type, extends(box) :: boundingbox_class
        type(vector) :: offset
        type(vector) :: center
    contains
    procedure :: initialize => initboundingbox
    procedure :: print      => printboundingbox
    end type boundingbox_class

    type(boundingbox_class), public :: BBox

    public :: boundingbox_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to initialize the simulation Bounding Box
    !---------------------------------------------------------------------------
    subroutine initboundingbox(self)
    implicit none
    class(boundingbox_class), intent(inout) :: self
    self%pt = Globals%SimDefs%Pointmin
    self%size = Globals%SimDefs%Pointmax - Globals%SimDefs%Pointmin
    self%offset = -self%pt !distance to the origin - local reference
    self%center = (Globals%SimDefs%Pointmin + Globals%SimDefs%Pointmax)/2.0
    call Globals%SimDefs%setCenter(self%center)
    end subroutine initboundingbox

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the simulation Bounding Box
    !---------------------------------------------------------------------------
    subroutine printboundingbox(self)
    implicit none
    class(boundingbox_class), intent(inout) :: self
    type(string) :: outext
    type(string) :: temp_str(3)
    outext = '-->Main bounding box is '//new_line('a')
    temp_str(1)=self%pt%x
    temp_str(2)=self%pt%y
    temp_str(3)=self%pt%z
    outext = outext//'       Point = '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')
    temp_str(1)=self%size%x
    temp_str(2)=self%size%y
    temp_str(3)=self%size%z
    outext = outext//'       Size = '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)
    call Log%put(outext,.false.)
    end subroutine printboundingbox

    end module boundingbox_mod
