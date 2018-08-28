    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : background
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> 
    !------------------------------------------------------------------------------

module background_mod
    
    use common_modules
    use field_types_mod
    
    implicit none
    private
            
    type :: background_class      !< a background solution class
    integer :: id
    type(box) :: extents                  !< shape::box that defines the extents of this background solution
    type(scalar2d_field_class), allocatable, dimension(:) :: scalar2d
    type(scalar3d_field_class), allocatable, dimension(:) :: scalar3d
    type(scalar4d_field_class), allocatable, dimension(:) :: scalar4d
    type(vectorial2d_field_class), allocatable, dimension(:) :: vectorial2d
    type(vectorial3d_field_class), allocatable, dimension(:) :: vectorial3d
    type(vectorial4d_field_class), allocatable, dimension(:) :: vectorial4d
    contains
    end type background_class

    
    !Public access vars
    public :: background_class

    contains
    
    
end module background_mod