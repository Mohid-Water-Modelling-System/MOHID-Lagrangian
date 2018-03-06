!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Lagrangian Tracer
! MODULE        : tracer3D
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Feb 2018
! REVISION      : Canelas 0.1
!> @author
!> Ricardo Birjukovs Canelas
!
! DESCRIPTION: 
!> Module that defines a pure Lagrangian tracer class and related methods.
!
! REVISION HISTORY:
! 26 02 2018 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
    
module tracer3D 

    use tracer_precision
    use tracer_interp
    
    use commom_modules

    implicit none
    private

    type tracer_par_trans_class             !>Type - transient parameters of a pure Lagrangian tracer object        
        character(len=512) :: par_trans_file 
        logical            :: use_par_trans
    end type 

    type tracer_par_class                   !>Type - parameters of a pure Lagrangian tracer object
        integer :: id                           !> unique tracer identification (integer)
        integer :: group                        !> Group to which the tracer belongs (usually by source)
        real(prec) :: vel_max                   !> Maximum velocity of tracer to track (m/s)
        logical    :: noise                     !  Add noise to location - TODO
        character(len=56) :: interp_method      !> interpolation method this tracer calls
        ! Transient parameters         
        type(tracer_par_trans_class) :: tpar    !> access to the transient parameters is done through this
    end type 

    type tracer_state_class             !>Type - state variables of a pure Lagrangian tracer object
        real(prec_time) :: time, age        ! time variables
        real(prec_time) :: dt
        logical :: active                   !> active switch
        type(vector) :: pos                 !> Position of the tracer (m)
        type(vector) :: vel                 !> Velocity of the tracer (m s-1)
        type(vector) :: acc                 !> Acceleration of the tracer (m s-2)
        real(prec) :: depth                 !> Depth of the tracer (m)
        real(prec) :: T                     !> Temperature of the tracer (Celcius)
    end type 

    type tracer_stats_class             !>Type - statistical variables of a pure Lagrangian tracer object        
        ! All stats variables at writing precision (prec_wrt)
        ! Avegarge variable is computed by Accumulated_var / ns
        type(vector) :: acc_pos             !> Accumulated position of the tracer (m)
        type(vector) :: acc_vel             !> Accumulated velocity of the tracer (m s-1)
        real(prec_wrt) :: acc_depth         !> Accumulated depth of the tracer (m)
        real(prec_wrt) :: acc_T             !> Accumulated temperature of the tracer (Celcius)
        integer :: ns                       !> Number of sampling steps
    end type    

    type tracer_class                   !>Type - a pure Lagrangian tracer class
        type(tracer_par_class)   :: par     !>To access parameters
        type(tracer_state_class) :: now     !>To access state variables
        type(tracer_stats_class) :: stats   !>To access statistics
    end type 

    ! For other tracer modules, altough we want them to inherit from tracer_class directly
    public :: tracer_par_class
    public :: tracer_state_class
    public :: tracer_stats_class

    ! General public 
    public :: tracer_class 
    public :: tracer_init

    contains 
    
   !---------------------------------------------------------------------------  
   !> @Ricardo Birjukovs Canelas - MARETEC
   ! Routine Author Name and Affiliation.
   ! 
   !> @brief
   !> Tracer inititialization routine - Generates a tracer and initializes its variables
   ! 
   !> @param[out] trc      
   !> @param[in] filename
   !---------------------------------------------------------------------------
    subroutine tracer_init(trc,id,time,x,y,z)

        implicit none

        type(tracer_class),   intent(OUT) :: trc
        integer, intent(IN) :: id
        real(prec), intent(IN) :: x, y, z  
        real(prec_time), intent(IN) :: time 

        ! Local variables 
        integer :: i         
        
        ! initialize parameters
        trc%par%id = id
        !vel_max, noise, interp_method
                
        ! Initialize position
        trc%now%pos = x*ex + y*ey + z*ez    !ex ... are versors along the x... axis
        !for depth, T we need info from the hydrodynamic solution and the mesh
        
        ! Initialize statistical accumulator variables
        
        !Initialize transient parameters
        
        return 

    end subroutine tracer_init   

end module tracer3D 


