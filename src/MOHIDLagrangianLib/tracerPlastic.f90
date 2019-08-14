    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_plastic
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a Lagrangian tracer class for plastic modelling and related methods.
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods
    !------------------------------------------------------------------------------

    module tracerPlastic_mod

    use tracerBase_mod
    use common_modules
    use sources_mod

    implicit none
    private

    type :: plastic_par_class               !<Type - parameters of a Lagrangian tracer object representing a plastic material
        real(prec) :: degradation_rate              !< degradation rate of the material
        logical    :: particulate                   !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec) :: size                          !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
    end type plastic_par_class

    type :: plastic_state_class             !<Type - State variables of a tracer object representing a plastic material
        real(prec) :: density                       !< density of the material
        real(prec) :: radius                        !< Tracer radius (m)
        real(prec) :: condition                     !< Material condition (1-0)
        real(prec) :: concentration                 !< Particle concentration
    end type plastic_state_class

    type, extends(tracer_class) :: plastic_class    !<Type - The plastic material Lagrangian tracer class
        type(plastic_par_class)   :: mpar     !<To access material parameters
        type(plastic_state_class) :: mnow     !<To access material state variables
    contains
    procedure :: getNumVars
    procedure :: getStateArray
    procedure :: setStateArray
    end type plastic_class

    !Public access vars
    public :: plastic_class

    !Public access routines
    public :: plasticTracer

    interface plasticTracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the number of variables used by this tracer
    !---------------------------------------------------------------------------
    integer function getNumVars(self)
    class(plastic_class), intent(in) :: self
    getNumVars = 15
    end function getNumVars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the state array of this tracer
    !---------------------------------------------------------------------------
    function getStateArray(self)
    class(plastic_class), intent(in) :: self
    real(prec), allocatable, dimension(:) :: getStateArray
    allocate(getStateArray(self%getNumVars()))
    getStateArray(1) = self%now%pos%x
    getStateArray(2) = self%now%pos%y
    getStateArray(3) = self%now%pos%z
    getStateArray(4) = self%now%vel%x
    getStateArray(5) = self%now%vel%y
    getStateArray(6) = self%now%vel%z
    getStateArray(7) = self%now%diffusionVel%x
    getStateArray(8) = self%now%diffusionVel%y
    getStateArray(9) = self%now%diffusionVel%z
    getStateArray(10) = self%now%usedMixingLenght
    getStateArray(11) = self%now%age
    getStateArray(12) = self%mnow%density
    getStateArray(13) = self%mnow%radius
    getStateArray(14) = self%mnow%condition
    getStateArray(15) = self%mnow%concentration
    end function getStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the state array of this tracer
    !---------------------------------------------------------------------------
    subroutine setStateArray(self, stateArray)
    class(plastic_class), intent(inout) :: self
    real(prec), dimension(:), intent(in) :: stateArray
    !if(size(stateArray)<self%getNumVars())
    self%now%pos%x = StateArray(1)
    self%now%pos%y = StateArray(2)
    self%now%pos%z = StateArray(3)
    self%now%vel%x = StateArray(4)
    self%now%vel%y = StateArray(5)
    self%now%vel%z = StateArray(6)
    self%now%diffusionVel%z = StateArray(7)
    self%now%diffusionVel%z = StateArray(8)
    self%now%diffusionVel%z = StateArray(9)
    self%now%usedMixingLenght = StateArray(10)
    self%now%age   = StateArray(11)
    self%mnow%density = StateArray(12)
    self%mnow%radius = StateArray(13)
    self%mnow%condition = StateArray(14)
    self%mnow%concentration = StateArray(15)
    end subroutine setStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Plastic Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    type(plastic_class) :: constructor
    integer, intent(in) :: id
    class(source_class), intent(in) :: src
    real(prec), intent(in) :: time
    integer, intent(in) :: p
    integer :: idx
    type(string) :: tag

    !use the base class constructor to build the base of our new derived type
    constructor%tracer_class = Tracer(id, src, time, p, constructor%getNumVars())
    !VERY NICE IFORT BUG (I think) - only some of the variables get used using the base constructor...
    constructor%par%id = id !forcing
    constructor%par%idsource = src%par%id !forcing
    !now initialize the specific components of this derived type
    constructor%par%ttype = Globals%Types%plastic
    !material parameters
    !constructor%mpar%degradation_rate = src%prop%degrd_rate
    constructor%mpar%particulate = src%prop%particulate
    !constructor%mpar%size = src%prop%radius
    !material state
    constructor%mnow%density = src%prop%density
    tag = 'condition'
    idx = Utils%find_str(src%prop%propName, tag, .true.)
    constructor%mnow%condition = src%prop%propValue(idx)
    constructor%mnow%radius = src%prop%radius
    !constructor%mnow%concentration = MV
    if (constructor%mpar%particulate) then
        !constructor%mpar%size = src%prop%pt_radius !correcting size to now mean particle size, not tracer size
        !constructor%mnow%concentration = src%prop%ini_concentration
    end if
    !filling the rest of the varName list
    constructor%varName(12) = Globals%Var%density
    constructor%varName(13) = 'radius'
    constructor%varName(14) = 'condition'
    constructor%varName(15) = 'concentration'
    end function constructor

    end module tracerPlastic_mod
