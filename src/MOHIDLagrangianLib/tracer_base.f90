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

    use common_modules
    use sources_mod

    implicit none
    private

    type :: tracer_par_class               !<Type - parameters of a pure Lagrangian tracer object
        integer :: id = MV                       !< unique tracer identification
        integer :: idsource = MV                 !< Source to which the tracer belongs
        real(prec) :: velmax = MV                !< Maximum velocity of tracer to track (m/s)
    end type tracer_par_class

    type :: tracer_state_class             !<Type - state variables of a pure Lagrangian tracer object
        real(prec) :: age = MV             ! time variables
        logical :: active = .false.             !< active switch
        type(vector) :: pos                     !< Position of the tracer (m)
        type(vector) :: vel                     !< Velocity of the tracer (m s-1)
        type(vector) :: acc                     !< Acceleration of the tracer (m s-2)
        real(prec) :: level = MV                !< Depth of the tracer (m)
    end type tracer_state_class

    type :: tracer_stats_class             !<Type - statistical variables of a pure Lagrangian tracer object
        type(vector) :: acc_pos                 !< Accumulated position of the tracer (m)
        type(vector) :: acc_vel                 !< Accumulated velocity of the tracer (m s-1)
        real(prec_wrt) :: acc_level = MV        !< Accumulated depth of the tracer (m)
        integer :: ns = MV                      !< Number of sampling steps
    end type tracer_stats_class

    type :: tracer_class                   !<Type - The pure Lagrangian tracer class
        type(tracer_par_class)   :: par         !<To access parameters
        type(tracer_state_class) :: now         !<To access state variables
        type(tracer_stats_class) :: stats       !<To access statistics
    contains
    procedure :: print => printTracer
    end type tracer_class
        
    !Public access vars
    public :: tracer_class

    !Public access routines
    public :: Tracer

    interface Tracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print basic info about the Tracer
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printTracer(self)
    implicit none
    class(tracer_class), intent(inout) :: self
    type(string) :: outext, t(6)
    if (self%now%active .eqv. .false.) then
        outext = '-->Tracer is inactive'
        call Log%put(outext,.false.)
    else
        t(1) = self%par%id
        t(2) = self%now%pos%x
        t(3) = self%now%pos%y
        t(4) = self%now%pos%z
        outext = 'Tracer['//t(1)//']::xyz('//t(2)//','//t(3)//','//t(4)//')'
        call Log%put(outext,.false.)
    end if
    end subroutine printTracer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Base Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    implicit none
    type(tracer_class) :: constructor
    integer, intent(in) :: id
    class(source_class), intent(in) :: src
    real(prec), intent(in) :: time
    integer, intent(in) :: p

    ! initialize parameters
    constructor%par%id = id
    constructor%par%idsource = src%par%id
    constructor%par%velmax = 15.0 !(m/s, just a placeholder)
    ! initialize tracer state
    constructor%now%age=0.0
    constructor%now%active = .true.
    !print*, 'Source at'
    !print*, src%now%pos
    !print*, 'New tracer at'
    !print*, src%stencil%ptlist(p) + src%now%pos
    constructor%now%pos = src%stencil%ptlist(p) + src%now%pos
    constructor%now%vel = 0.0
    constructor%now%acc = 0.0
    constructor%now%level = 0.0
    ! Initialize statistical accumulator variables
    constructor%stats%acc_pos = 0.0
    constructor%stats%acc_vel = 0.0
    constructor%stats%acc_level = 0.0
    constructor%stats%ns = 0

    end function constructor

    end module tracer_base_mod
