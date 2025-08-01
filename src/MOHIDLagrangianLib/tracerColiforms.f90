    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_coliform
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : Colab +Atlantic
    ! DATE          : Dec 2020
    ! REVISION      : Sobrinho 0.1
    !> @author
    !> Joao Sobrinho
    !
    ! DESCRIPTION:
    !> Module that defines a Lagrangian tracer class for fecal coliforms modelling and related methods.
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods
    !------------------------------------------------------------------------------

    module tracerColiform_mod

    use tracerBase_mod
    use common_modules
    use sources_mod

    implicit none
    private

    type :: coliform_par_class               !<Type - parameters of a Lagrangian tracer object representing a coliform cell
        integer    :: particulate                   !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec) :: size                          !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
    end type coliform_par_class

    type :: coliform_state_class             !<Type - State variables of a tracer object representing a coliform material
        real(prec) :: density                       !< density of the material
        real(prec) :: radius                        !< Tracer radius (m)
        real(prec) :: volume                        !< Tracer volume (m3)
        real(prec) :: initial_volume                !< Tracer initial volume (m3)
        real(prec) :: area                          !< Tracer area (m2)
        real(prec) :: condition                     !< Material condition (1-0)
        real(prec) :: T90                           !< T90 mortality of the coliforms (seconds)
        real(prec) :: sw_extinction_coef            !< Short wave radiation extinction coef (1/m)
        real(prec) :: sw_percentage                 !< Short wave radiation percentage of incoming radiation (0-1)
        integer    :: T90_variable                  !< Variable T90 decay 
        integer    :: T90_method                    !< Fecal decay according to 1: Canteras et al. (1995). 2: Chapra (1997)
        real(prec) :: concentration                 !< Particle concentration
        real(prec) :: temperature                   !< temperature of the tracer
        real(prec) :: salinity                      !< salinity of the tracer
        !logical :: beachPeriod                      !< consecutive period of time (in seconds) that the tracer has been beached
        !integer :: beachAreaId                      !< beaching area Id where the tracer last beached
        !integer :: beachedWaterLevel                !< Water level at the time the tracer was beachded
    end type coliform_state_class

    type, extends(tracer_class) :: coliform_class    !<Type - The coliform material Lagrangian tracer class
        type(coliform_par_class)   :: mpar     !<To access material parameters
        type(coliform_state_class) :: mnow     !<To access material state variables
    contains
    procedure :: getNumVars
    procedure :: getStateArray
    procedure :: setStateArray
    end type coliform_class

    !Public access vars
    public :: coliform_class

    !Public access routines
    public :: coliformTracer

    interface coliformTracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the number of variables used by this tracer
    !---------------------------------------------------------------------------
    integer function getNumVars(self)
    class(coliform_class), intent(in) :: self
    getNumVars = 32
    end function getNumVars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the state array of this tracer
    !---------------------------------------------------------------------------
    function getStateArray(self)
    class(coliform_class), intent(in) :: self
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
    getStateArray(12) = self%mpar%particulate
    getStateArray(13) = self%now%bathymetry
    getStateArray(14) = self%now%dwz
    getStateArray(15) = self%now%dist2bottom
    getStateArray(16) = self%now%beachPeriod
    getStateArray(17) = self%now%beachAreaId
    getStateArray(18) = self%now%beachedWaterLevel
    getStateArray(19) = self%mnow%density
    getStateArray(20) = self%mnow%radius
    getStateArray(21) = self%mnow%volume
    getStateArray(22) = self%mnow%area
    getStateArray(23) = self%mnow%condition
    getStateArray(24) = self%mnow%T90
    getStateArray(25) = self%mnow%T90_variable
    getStateArray(26) = self%mnow%T90_method
    getStateArray(27) = self%mnow%sw_percentage
    getStateArray(28) = self%mnow%sw_extinction_coef
    getStateArray(29) = self%mnow%concentration
    getStateArray(30) = self%mnow%initial_volume
    getStateArray(31) = self%mnow%temperature
    getStateArray(32) = self%mnow%salinity
    end function getStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the state array of this tracer
    !---------------------------------------------------------------------------
    subroutine setStateArray(self, stateArray)
    class(coliform_class), intent(inout) :: self
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
    self%mpar%particulate = StateArray(12)
    self%now%bathymetry   = StateArray(13)
    self%now%dwz          = StateArray(14)
    self%now%dist2bottom = StateArray(15)
    self%now%beachPeriod = StateArray(16)
    self%now%beachAreaId = StateArray(17)
    self%now%beachedWaterLevel = StateArray(18)
    self%mnow%density = StateArray(19)
    self%mnow%radius = StateArray(20)
    self%mnow%volume = StateArray(21)
    self%mnow%area = StateArray(22)
    self%mnow%condition = StateArray(23)
    self%mnow%T90        = StateArray(24)
    self%mnow%T90_variable  = StateArray(25)
    self%mnow%T90_method    = StateArray(26)
    self%mnow%sw_percentage = StateArray(27)
    self%mnow%sw_extinction_coef = StateArray(28)
    self%mnow%concentration = StateArray(29)
    self%mnow%initial_volume = StateArray(30)
    self%mnow%temperature = StateArray(31)
    self%mnow%salinity = StateArray(32)
    end subroutine setStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> coliform Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    type(coliform_class) :: constructor
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
    constructor%par%ttype = Globals%Types%coliform
    constructor%mpar%particulate = src%prop%particulate
    constructor%mpar%size = src%prop%radius
    
    !material state
    constructor%mnow%density = src%prop%density
    constructor%mnow%radius = src%prop%radius
    constructor%mnow%volume = src%prop%volume
    constructor%mnow%area = src%prop%area
    !default values
    constructor%mnow%T90 = 1/7200
    constructor%mnow%condition = 1.0
    constructor%mnow%T90_variable = 0
    constructor%mnow%concentration = 1000000
    
    constructor%mnow%temperature = 15.0
    constructor%mnow%salinity = 36.0
    
    !try to find value from material types files
    tag = 'condition'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%condition = src%prop%propValue(idx)
    end if
    tag = 'T90'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%T90 = 1 / src%prop%propValue(idx)
    end if
    
    tag = 'T90_variable'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%T90_variable = src%prop%propValue(idx)
    end if
    
    tag = 'T90_method'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%T90_method = src%prop%propValue(idx)
    end if
    
    tag = 'sw_percentage'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%sw_percentage= src%prop%propValue(idx)
    end if
    
    tag = 'sw_extinction_coef'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%sw_extinction_coef = src%prop%propValue(idx)
    end if
    
    tag = 'initial_volume'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%initial_volume = src%prop%propValue(idx)
    end if
    
    tag = 'concentration'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%concentration = src%prop%propValue(idx)
    end if
    
    if (constructor%mpar%particulate==1) then
        !constructor%mpar%size = src%prop%pt_radius !correcting size to now mean particle size, not tracer size
        !constructor%mnow%concentration = src%prop%ini_concentration
    end if
    
    !filling the rest of the varName list
    constructor%varName(19) = Globals%Var%density
    constructor%varName(20) = 'radius'
    constructor%varName(21) = 'volume'
    constructor%varName(22) = 'area'
    constructor%varName(23) = 'condition'
    constructor%varName(24) = 'T90'
    constructor%varName(25) = 'T90_variable'
    constructor%varName(26) = 'T90_method'
    constructor%varName(27) = 'sw_percentage'
    constructor%varName(28) = 'sw_extinction_coef'
    constructor%varName(29) = 'concentration'
    constructor%varName(30) = 'initial_volume'
    !constructor%varName(27) = 'particulate'
    constructor%varName(31) = 'temp'
    constructor%varName(32) = 'salt'
    
    end function constructor

    end module tracerColiform_mod
