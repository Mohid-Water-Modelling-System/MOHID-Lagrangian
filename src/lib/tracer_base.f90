    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_base
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a pure Lagrangian tracer class and related methods.
    !------------------------------------------------------------------------------

    module tracer_base_mod

    use tracer_interp_mod
    use commom_modules

    implicit none
    private

    type tracer_par_class               !<Type - parameters of a pure Lagrangian tracer object
        integer :: id                       !< unique tracer identification
        integer :: idsource                 !< Source to which the tracer belongs
        real(prec) :: velmax                !< Maximum velocity of tracer to track (m/s)
        logical    :: noise                 !  Add noise to location
        type(string) :: interp_method       !< interpolation method this tracer calls
    end type

    type tracer_state_class             !<Type - state variables of a pure Lagrangian tracer object
        real(prec_time) :: age              ! time variables
        logical :: active                   !< active switch
        type(vector) :: pos                 !< Position of the tracer (m)
        type(vector) :: vel                 !< Velocity of the tracer (m s-1)
        type(vector) :: acc                 !< Acceleration of the tracer (m s-2)
        real(prec) :: depth                 !< Depth of the tracer (m)
        !real(prec) :: T                     !< Temperature of the tracer (Celcius)
    end type

    type tracer_stats_class             !<Type - statistical variables of a pure Lagrangian tracer object
        ! All stats variables at writing precision (prec_wrt)
        ! Avegarge variable is computed by Accumulated_var / ns
        type(vector) :: acc_pos             !< Accumulated position of the tracer (m)
        type(vector) :: acc_vel             !< Accumulated velocity of the tracer (m s-1)
        real(prec_wrt) :: acc_depth         !< Accumulated depth of the tracer (m)
        !real(prec_wrt) :: acc_T             !< Accumulated temperature of the tracer (Celcius)
        integer :: ns                       !< Number of sampling steps
    end type

    type tracer_class                   !<Type - The pure Lagrangian tracer class
        type(tracer_par_class)   :: par     !<To access parameters
        type(tracer_state_class) :: now     !<To access state variables
        type(tracer_stats_class) :: stats   !<To access statistics
    contains
    procedure :: initialize
    end type

    !Simulation variables
    type(tracer_class) :: dummyTracer !< Just a template to allocate the generic arrays to this size

    !Public access vars
    public :: tracer_class, dummyTracer 

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Tracer initialization method
    !
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine initialize(trc,id,id_source,time,pt)
    implicit none
    class(tracer_class) :: trc
    integer, intent(in) :: id
    integer, intent(in) :: id_source
    real(prec_time), intent(in) :: time
    type(vector), intent(in) :: pt

    ! initialize parameters
    trc%par%id = id
    trc%par%idsource = id_source
    trc%par%velmax = 15.0 !(m/s, just a placeholder)
    ! interp_method - TODO
    ! initialize tracer state
    trc%now%age=0.0
    trc%now%active = .false.
    trc%now%pos = pt
    trc%now%vel = 0.0
    trc%now%acc = 0.0
    trc%now%depth = 0.0
    ! Initialize statistical accumulator variables
    trc%stats%acc_pos = 0.0
    trc%stats%acc_vel = 0.0
    trc%stats%acc_depth = 0.0
    trc%stats%ns = 0

    end subroutine initialize

  end module tracer_base_mod
