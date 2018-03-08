!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : source
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module that defines geometry classes and related methods.
!------------------------------------------------------------------------------
    
module geometry
    
    use vecfor
    use tracer_precision
    
    implicit none
    private
    
    type point              !>Type - point class
        type(vector) :: coord   !> Coordinates of the point
    end type
    
    type line               !>Type - line class
        type(vector) :: first   !> Coordinates of the initial point
        type(vector) :: last    !> Coordinates of the end point
    end type
    
    type sphere             !>Type - sphere class
        type(vector) :: coord   !> Coordinates of the point
        real(prec) :: radius    !> Sphere radius
    end type
    
    type box                !>Type - point class
        type(vector) :: coord   !> Coordinates of the lower left corner point
        type(vector) :: size    !> Box size
    end type
    
    
    
    
    
    
end module geometry 
