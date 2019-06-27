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

    module tracerBase_mod

    use common_modules
    use sources_mod

    implicit none
    private

    type :: tracer_par_class               !<Type - parameters of a pure Lagrangian tracer object
        integer :: id = MV                       !< unique tracer identification
        integer :: idsource = MV                 !< Source to which the tracer belongs
        integer :: ttype
    end type tracer_par_class

    type :: tracer_state_class             !<Type - state variables of a pure Lagrangian tracer object
        real(prec) :: age = MV             ! time variables
        logical :: active = .false.             !< active switch
        type(vector) :: pos                     !< Position of the tracer (m)
        type(vector) :: vel                     !< Velocity of the tracer (m s-1)        
    end type tracer_state_class

    type :: tracer_class                   !<Type - The pure Lagrangian tracer class
        type(tracer_par_class)   :: par         !<To access parameters
        type(tracer_state_class) :: now         !<To access state variables
    contains
    procedure :: getNumVars
    procedure :: getStateArray
    procedure :: setStateArray
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
    !> Method that returns the number of variables used by this tracer
    !---------------------------------------------------------------------------
    integer function getNumVars(self)
    class(tracer_class), intent(in) :: self
    getNumVars = 6
    end function getNumVars
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the state array of this tracer
    !---------------------------------------------------------------------------
    function getStateArray(self)
    class(tracer_class), intent(in) :: self
    real(prec), allocatable, dimension(:) :: getStateArray
    allocate(getStateArray(self%getNumVars()))
    getStateArray(1) = self%now%pos%x
    getStateArray(2) = self%now%pos%y
    getStateArray(3) = self%now%pos%z
    getStateArray(4) = self%now%vel%x
    getStateArray(5) = self%now%vel%y
    getStateArray(6) = self%now%vel%z
    end function getStateArray
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the state array of this tracer
    !---------------------------------------------------------------------------
    subroutine setStateArray(self, stateArray)
    class(tracer_class), intent(inout) :: self
    real(prec), dimension(:), intent(in) :: stateArray
    !if(size(stateArray)<self%getNumVars())
    self%now%pos%x = StateArray(1)
    self%now%pos%y = StateArray(2)
    self%now%pos%z = StateArray(3)
    self%now%vel%x = StateArray(4)
    self%now%vel%y = StateArray(5)
    self%now%vel%z = StateArray(6)
    end subroutine setStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print basic info about the Tracer
    !---------------------------------------------------------------------------
    subroutine printTracer(self)
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
    type(tracer_class) :: constructor
    integer, intent(in) :: id
    class(source_class), intent(in) :: src
    real(prec), intent(in) :: time
    integer, intent(in) :: p
    ! initialize parameters
    constructor%par%id = id
    constructor%par%idsource = src%par%id
    constructor%par%ttype = Globals%Types%base
    ! initialize tracer state
    constructor%now%age=0.0
    constructor%now%active = .true.
    constructor%now%pos = src%stencil%ptlist(p) + src%now%pos
    constructor%now%vel = 0.0
    end function constructor

    end module tracerBase_mod
