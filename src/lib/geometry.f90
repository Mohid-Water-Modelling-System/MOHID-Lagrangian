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
    
    type shape                                  !>Type - extendable shape class
        type(vector) :: pt      !> Coordinates of a point
    end type
    
    type, extends(shape) :: point               !>Type - point class
    end type
    
    type, extends(shape) :: line               !>Type - line class
        type(vector) :: last    !> Coordinates of the end point
    end type
    
    type, extends(shape) :: sphere             !>Type - sphere class
        real(prec) :: radius    !> Sphere radius
    end type
    
    type, extends(shape) :: box                !>Type - point class
        !> Coordinates of the lower left corner point are defined by shape class
        type(vector) :: size    !> Box size
    end type
    
    
    
    
    
    
end module geometry 
