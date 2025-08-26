    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_detritus
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a Lagrangian tracer class for detritus modelling and related methods.
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods
    !------------------------------------------------------------------------------

    module tracerdetritus_mod

    use tracerBase_mod
    use common_modules
    use sources_mod

    implicit none
    private

    type :: detritus_par_class               !<Type - parameters of a Lagrangian tracer object representing a detritus material
        integer    :: particulate                   !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec) :: size                          !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
    end type detritus_par_class

    type :: detritus_state_class             !<Type - State variables of a tracer object representing a detritus material
        real(prec) :: density                       !< density of the material
        real(prec) :: radius                        !< Tracer radius (m)
        real(prec) :: volume                        !< Tracer volume (m3)
        real(prec) :: area                          !< Tracer area (m2)
        real(prec) :: condition                     !< Material condition (1-0)
        real(prec) :: initial_volume                !< initial volume of tracer
        real(prec) :: temperature                   !< temperature of the tracer
        real(prec) :: salinity                      !< salinity of the tracer
        !logical    :: beachPeriod                   !< consecutive period of time (in seconds) that the tracer has been beached
        !integer    :: beachAreaId                   !< beaching area Id where the tracer last beached
        !integer    :: beachedWaterLevel             !< Water level at the time the tracer was beachded
    end type detritus_state_class

    type, extends(tracer_class) :: detritus_class    !<Type - The detritus material Lagrangian tracer class
        type(detritus_par_class)   :: mpar     !<To access material parameters
        type(detritus_state_class) :: mnow     !<To access material state variables
    contains
    procedure :: getNumVars
    procedure :: getStateArray
    procedure :: setStateArray
    end type detritus_class

    !Public access vars
    public :: detritus_class

    !Public access routines
    public :: detritusTracer

    interface detritusTracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the number of variables used by this tracer
    !---------------------------------------------------------------------------
    integer function getNumVars(self)
    class(detritus_class), intent(in) :: self
    getNumVars = 28
    end function getNumVars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the state array of this tracer
    !---------------------------------------------------------------------------
    function getStateArray(self)
    class(detritus_class), intent(in) :: self
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
    getStateArray(11) = self%now%VelStandardDeviation
    getStateArray(12) = self%now%TPathHor
    getStateArray(13) = self%now%age
    getStateArray(14) = self%mpar%particulate
    getStateArray(15) = self%now%bathymetry
    getStateArray(16) = self%now%dwz
    getStateArray(17) = self%now%dist2bottom
    getStateArray(18) = self%now%beachPeriod
    getStateArray(19) = self%now%beachAreaId
    getStateArray(20) = self%now%beachedWaterLevel
    getStateArray(21) = self%mnow%density
    getStateArray(22) = self%mnow%radius
    getStateArray(23) = self%mnow%volume
    getStateArray(24) = self%mnow%area
    getStateArray(25) = self%mnow%condition
    getStateArray(26) = self%mnow%initial_volume
    getStateArray(27) = self%mnow%temperature
    getStateArray(28) = self%mnow%salinity
    end function getStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the state array of this tracer
    !---------------------------------------------------------------------------
    subroutine setStateArray(self, stateArray)
    class(detritus_class), intent(inout) :: self
    real(prec), dimension(:), intent(in) :: stateArray
    !if(size(stateArray)<self%getNumVars())
    self%now%pos%x = StateArray(1)
    self%now%pos%y = StateArray(2)
    self%now%pos%z = StateArray(3)
    self%now%vel%x = StateArray(4)
    self%now%vel%y = StateArray(5)
    self%now%vel%z = StateArray(6)
    self%now%diffusionVel%x = StateArray(7)
    self%now%diffusionVel%y = StateArray(8)
    self%now%diffusionVel%z = StateArray(9)
    self%now%usedMixingLenght = StateArray(10)
    self%now%VelStandardDeviation = StateArray(11)
    self%now%TPathHor = StateArray(12)
    self%now%age          = StateArray(13)
    self%mpar%particulate = StateArray(14)
    self%now%bathymetry   = StateArray(15)
    self%now%dwz          = StateArray(16)
    self%now%dist2bottom = StateArray(17)
    self%now%beachPeriod = StateArray(18)
    self%now%beachAreaId = StateArray(19)
    self%now%beachedWaterLevel = StateArray(20)
    self%mnow%density = StateArray(21)
    self%mnow%radius = StateArray(22)
    self%mnow%volume = StateArray(23)
    self%mnow%area = StateArray(24)
    self%mnow%condition = StateArray(25)
    self%mnow%initial_volume = StateArray(26)
    self%mnow%temperature = StateArray(27)
    self%mnow%salinity = StateArray(28)
    end subroutine setStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> detritus Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    type(detritus_class) :: constructor
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
    constructor%par%ttype = Globals%Types%detritus
    constructor%mpar%particulate = src%prop%particulate
    constructor%mpar%size = src%prop%radius
    !material state
    constructor%mnow%density = src%prop%density
    constructor%mnow%radius = src%prop%radius
    constructor%mnow%volume = src%prop%volume
    constructor%mnow%area = src%prop%area
    constructor%mnow%initial_volume = src%prop%volume
    !default values
    constructor%mnow%condition = 1.0
    !constructor%mnow%degradation_rate = 1/(100*365*24*3600)
    constructor%mnow%temperature = 15.0
    constructor%mnow%salinity = 36.0
    
    !try to find value from material types files
    tag = 'condition'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%condition = src%prop%propValue(idx)
    end if

    if (constructor%mpar%particulate==1) then
        !constructor%mpar%size = src%prop%pt_radius !correcting size to now mean particle size, not tracer size
        !constructor%mnow%concentration = src%prop%ini_concentration
    end if
    
    !filling the rest of the varName list
    constructor%varName(21) = Globals%Var%density
    constructor%varName(22) = 'radius'
    constructor%varName(23) = 'volume'
    constructor%varName(24) = 'area'
    constructor%varName(25) = 'condition'
    !constructor%varName(20) = 'particulate'
    constructor%varName(26) = 'initial_volume'
    constructor%varName(27) = 'temp'
    constructor%varName(28) = 'salt'
    end function constructor

    end module tracerdetritus_mod
