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
    use simulation_precision
    
    implicit none
    
    !>Type - extendable shape class
    type shape
        !> Coordinates of a point
        type(vector) :: pt      
    end type
    
    !>Type - point class
    type, extends(shape) :: point               
    end type
    
    !>Type - line class
    type, extends(shape) :: line
        !> Coordinates of the end point
        type(vector) :: last    
    end type
    
    !>Type - sphere class
    type, extends(shape) :: sphere  
        !> Sphere radius
        real(prec) :: radius    
    end type
    
    !>Type - point class
    type, extends(shape) :: box                
        ! Coordinates of the lower left corner point are defined by shape class
        !> Box size
        type(vector) :: size    
    end type
    
    
    
    
    
    
end module geometry 
