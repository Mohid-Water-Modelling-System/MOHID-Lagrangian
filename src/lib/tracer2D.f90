!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : tracer2D
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module that defines a pure Lagrangian 2D tracer class and related methods, as a subset of the tracer3D module.
!
! REVISION HISTORY:
! 26 02 2018 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
    
module tracer2D 
    ! Module to wrap a 2D (profile) version of the tracer model
    ! Makes calls to main tracer code by reshaping profile into
    ! 3D array with y-dimension thickness of 1. 

    use tracer_precision
    use tracer3D
    
    implicit none 

contains

   !---------------------------------------------------------------------------  
   !> @Ricardo Birjukovs Canelas - MARETEC
   !> Routine Author Name and Affiliation.
   !
   ! DESCRIPTION: 
   !> Brief description of routine. 
   !> @brief
   !> 2D Tracer inititialization routine - Generates a tracer collection and initializes their variables
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[out] trc      
   !> @param[in] filename
   !---------------------------------------------------------------------------
    subroutine tracer2D_init(trc,filename,time,x,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:)
        logical,    intent(IN) :: is_sigma 
        real(prec_time) :: time 

        real(prec) :: y(5) 

        ! Define the ghost y-dimension
        y(1:5) = [0.0,0.25,0.50,0.75,1.0] 

        ! Call 3D tracer_init
        call tracer_init(trc,filename,time,x,y,is_sigma)

        return 

    end subroutine tracer2D_init

end module tracer2D

