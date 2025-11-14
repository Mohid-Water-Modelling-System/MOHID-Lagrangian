    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group,  modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         :  Model
    ! PROJECT       :  Lagrangian Tracer
    ! MODULE        : tracer_waterQuality
    ! URL           : http://www..com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a Lagrangian tracer class for waterQuality modelling and related methods.
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods
    !------------------------------------------------------------------------------

    module tracerwaterQuality_mod

    use tracerBase_mod
    use common_modules
    use sources_mod

    implicit none
    private

    type :: waterQuality_par_class               !<Type - parameters of a Lagrangian tracer object representing a waterQuality material
        integer    :: particulate                   !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec) :: size                          !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
    end type waterQuality_par_class

    type :: waterQuality_state_class             !<Type - State variables of a tracer object representing a waterQuality material
        real(prec) :: density                       !< density of the material
        real(prec) :: radius                        !< Tracer radius (m)
        real(prec) :: volume                        !< Tracer volume (m3)
        real(prec) :: area                          !< Tracer area (m2)
        real(prec) :: condition                     !< Material condition (1-0)
        real(prec) :: concentration                 !< Particle concentration
        real(prec) :: oxygen                        !< oxygen of the tracer
        real(prec) :: ammonia                       !< ammonia of the tracer
        real(prec) :: nitrate                       !< nitrate of the tracer
        real(prec) :: nitrite                       !< nitrite of the tracer
        real(prec) :: inorganic_phosphorus          !< inorganic_phosphorus of the tracer
        real(prec) :: DON_NonRefractory             !< DON_NonRefractory of the tracer
        real(prec) :: DOP_NonRefractory             !< DOP_NonRefractory of the tracer
        real(prec) :: DON_Refractory                !< DON_Refractory of the tracer
        real(prec) :: DOP_Refractory                !< DOP_Refractory of the tracer
        real(prec) :: partOrgNit                !< partOrgNit of the tracer
        real(prec) :: partOrgPho                !< partOrgPho of the tracer
        real(prec) :: phytoplankton                 !< phytoplankton of the tracer
        real(prec) :: zooplankton                   !< zooplankton of the tracer
        real(prec) :: initial_volume                !< Tracer initial volume (m3)
        real(prec) :: temperature                   !< temperature of the tracer
        real(prec) :: salinity                      !< salinity of the tracer
        real(prec) :: radius_cr_min                 !< Tracer min critical radius (m)
        real(prec) :: radius_cr_max                 !< Tracer max critical radius (m)
    end type waterQuality_state_class

    type, extends(tracer_class) :: waterQuality_class    !<Type - The waterQuality material Lagrangian tracer class
        type(waterQuality_par_class)   :: mpar     !<To access material parameters
        type(waterQuality_state_class) :: mnow     !<To access material state variables
    contains
    procedure :: getNumVars
    procedure :: getStateArray
    procedure :: setStateArray
    end type waterQuality_class

    !Public access vars
    public :: waterQuality_class

    !Public access routines
    public :: waterQualityTracer

    interface waterQualityTracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com		
    !> @brief
    !> Method that returns the number of variables used by this tracer
    !---------------------------------------------------------------------------
    integer function getNumVars(self)
    class(waterQuality_class), intent(in) :: self
    getNumVars = 46
    end function getNumVars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Method that returns the state array of this tracer
    !---------------------------------------------------------------------------
    function getStateArray(self)
    class(waterQuality_class), intent(in) :: self
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
    getStateArray(16) = self%now%rugosityVar
    getStateArray(17) = self%now%D50Var
    getStateArray(18) = self%now%dwz
    getStateArray(19) = self%now%dist2bottom
    getStateArray(20) = self%now%beachPeriod
    getStateArray(21) = self%now%beachAreaId
    getStateArray(22) = self%now%beachedWaterLevel
    getStateArray(23) = self%mnow%density
    getStateArray(24) = self%mnow%radius
    getStateArray(25) = self%mnow%volume
    getStateArray(26) = self%mnow%area
    getStateArray(27) = self%mnow%condition
    getStateArray(28) = self%mnow%concentration
    getStateArray(29) = self%mnow%oxygen
    getStateArray(30) = self%mnow%ammonia
    getStateArray(31) = self%mnow%nitrate
    getStateArray(32) = self%mnow%nitrite
    getStateArray(33) = self%mnow%inorganic_phosphorus
    getStateArray(34) = self%mnow%DON_NonRefractory
    getStateArray(35) = self%mnow%DOP_NonRefractory
    getStateArray(36) = self%mnow%DON_Refractory
    getStateArray(37) = self%mnow%DOP_Refractory
    getStateArray(38) = self%mnow%partOrgNit
    getStateArray(39) = self%mnow%partOrgPho
    getStateArray(40) = self%mnow%phytoplankton
    getStateArray(41) = self%mnow%zooplankton
    getStateArray(42) = self%mnow%initial_volume
    getStateArray(43) = self%mnow%temperature
    getStateArray(44) = self%mnow%salinity
    getStateArray(45) = self%mnow%radius_cr_min
    getStateArray(46) = self%mnow%radius_cr_max	
    end function getStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com		
    !> @brief
    !> Method that sets the state array of this tracer
    !---------------------------------------------------------------------------
    subroutine setStateArray(self, stateArray)
    class(waterQuality_class), intent(inout) :: self
    real(prec), dimension(:), intent(in) :: stateArray
    !if(size(stateArray)<self%getNumVars())
    self%now%pos%x 					= StateArray(1)
    self%now%pos%y 					= StateArray(2)
    self%now%pos%z 					= StateArray(3)
    self%now%vel%x 					= StateArray(4)
    self%now%vel%y 					= StateArray(5)
    self%now%vel%z 					= StateArray(6)
    self%now%diffusionVel%x 		= StateArray(7)
    self%now%diffusionVel%y 		= StateArray(8)
    self%now%diffusionVel%z 		= StateArray(9)
    self%now%usedMixingLenght 		= StateArray(10)
    self%now%VelStandardDeviation 	= StateArray(11)
    self%now%TPathHor 				= StateArray(12)
    self%now%age   					= StateArray(13)
    self%mpar%particulate 			= StateArray(14)
    self%now%bathymetry   			= StateArray(15)
    self%now%rugosityVar   			= StateArray(16)
    self%now%D50Var   				= StateArray(17)
    self%now%dwz          			= StateArray(18)
    self%now%dist2bottom 			= StateArray(19)
    self%now%beachPeriod 			= StateArray(20)
    self%now%beachAreaId 			= StateArray(21)
    self%now%beachedWaterLevel 		= StateArray(22)
    self%mnow%density 				= StateArray(23)
    self%mnow%radius				= StateArray(24)
    self%mnow%volume 				= StateArray(25)
    self%mnow%area 					= StateArray(26)
    self%mnow%condition 			= StateArray(27)
    self%mnow%concentration 		= StateArray(28)
    self%mnow%oxygen 				= StateArray(29)
    self%mnow%ammonia 				= StateArray(30)
    self%mnow%nitrate 				= StateArray(31) 
    self%mnow%nitrite 				= StateArray(32)
    self%mnow%inorganic_phosphorus 	= StateArray(33)
    self%mnow%DON_NonRefractory 	= StateArray(34)
    self%mnow%DOP_NonRefractory 	= StateArray(35)
    self%mnow%DON_Refractory 		= StateArray(36)
    self%mnow%DOP_Refractory 		= StateArray(37)
    self%mnow%partOrgNit 			= StateArray(38)
    self%mnow%partOrgPho 			= StateArray(39)
    self%mnow%phytoplankton 		= StateArray(40)
    self%mnow%zooplankton 			= StateArray(41)
    self%mnow%initial_volume 		= StateArray(42)	
    self%mnow%temperature 			= StateArray(43)
    self%mnow%salinity 				= StateArray(44)
	self%mnow%radius_cr_min 		= StateArray(45)
    self%mnow%radius_cr_max 		= StateArray(46)
    end subroutine setStateArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com		
    !> @brief
    !> waterQuality Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    type(waterQuality_class) :: constructor
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
    constructor%par%ttype = Globals%Types%waterQuality
    constructor%mpar%particulate = src%prop%particulate
    constructor%mpar%size = src%prop%radius
    !material state
    constructor%mnow%density = src%prop%density
    constructor%mnow%radius = src%prop%radius
    constructor%mnow%volume = src%prop%volume
    constructor%mnow%area = src%prop%area
    
	!default values
    constructor%mnow%condition = 1.0
    constructor%mnow%concentration = 1000000		! TODO: the value should comaptible with a correct value!
    constructor%mnow%oxygen = 5.0
    constructor%mnow%ammonia = 0.01
    constructor%mnow%nitrate = 0.01
    constructor%mnow%nitrite = 0.001
    constructor%mnow%inorganic_phosphorus = 0.01
    constructor%mnow%DON_NonRefractory = 0.001
    constructor%mnow%DOP_NonRefractory = 0.001
    constructor%mnow%DON_Refractory = 0.001
    constructor%mnow%DOP_Refractory = 0.001
    constructor%mnow%partOrgNit = 0.001
    constructor%mnow%partOrgPho = 0.001
    constructor%mnow%phytoplankton = 0.001    
    constructor%mnow%zooplankton = 0.001

	constructor%mnow%initial_volume = src%prop%volume
    constructor%mnow%temperature = 15.0					! TODO: the value should comaptible with a correct value!
    constructor%mnow%salinity = 36.0					! TODO: the value should comaptible with a correct value!
	
    constructor%mnow%radius_cr_min = 1.0e-4_prec * src%prop%radius	! TODO: the value should comaptible with a correct value!
    constructor%mnow%radius_cr_max = 1.0e+4_prec * src%prop%radius	! TODO: the value should comaptible with a correct value!
    
    !try to find value from material types files
    tag = 'condition'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%condition = src%prop%propValue(idx)
    end if

    tag = 'concentration'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%concentration = src%prop%propValue(idx)
    end if
	
    tag = 'dissolved_oxygen'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%oxygen = src%prop%propValue(idx)
    end if

    tag = 'ammonia'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%ammonia = src%prop%propValue(idx)
    end if

    tag = 'nitrate'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%nitrate = src%prop%propValue(idx)
    end if

    tag = 'nitrite'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%nitrite = src%prop%propValue(idx)
    end if

    tag = 'inorganic_phosphorus'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%inorganic_phosphorus = src%prop%propValue(idx)
    end if

    tag = 'DON_NonRefractory'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%DON_NonRefractory = src%prop%propValue(idx)
    end if

    tag = 'DOP_NonRefractory'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%DOP_NonRefractory = src%prop%propValue(idx)
    end if
	
    tag = 'DON_Refractory'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%DON_Refractory = src%prop%propValue(idx)
    end if

    tag = 'DOP_Refractory'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%DOP_Refractory = src%prop%propValue(idx)
    end if

    tag = 'partOrgNit'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%partOrgNit = src%prop%propValue(idx)
    end if

    tag = 'partOrgPho'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%partOrgPho = src%prop%propValue(idx)
    end if

    tag = 'phytoplankton'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%phytoplankton = src%prop%propValue(idx)
    end if

    tag = 'zooplankton'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%zooplankton = src%prop%propValue(idx)
    end if

    tag = 'initial_volume'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%initial_volume = src%prop%propValue(idx)
    end if
    
    tag = 'temp'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%temperature = src%prop%propValue(idx)
    end if
	
    tag = 'salt'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%salinity = src%prop%propValue(idx)
    end if

    tag = 'radius_cr_min'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%radius_cr_min = src%prop%propValue(idx)
    end if
	
    tag = 'radius_cr_max'
    idx = Utils%find_str(src%prop%propName, tag, .false.)
    if (idx /= MV_INT) then
        constructor%mnow%radius_cr_max = src%prop%propValue(idx)
    end if
	
    if (constructor%mpar%particulate==1) then
        !constructor%mpar%size = src%prop%pt_radius !correcting size to now mean particle size, not tracer size
        !constructor%mnow%concentration = src%prop%ini_concentration
    end if
    
    !filling the rest of the varName list
    constructor%varName(23) = Globals%Var%density
    constructor%varName(24) = 'radius'
    constructor%varName(25) = 'volume'
    constructor%varName(26) = 'area'
    constructor%varName(27) = 'condition'
    constructor%varName(28) = 'concentration'
    !constructor%varName(19) = 'particulate'
    constructor%varName(29) = 'dissolved_oxygen'
    constructor%varName(30) = 'ammonia'
    constructor%varName(31) = 'nitrate'
    constructor%varName(32) = 'nitrite'
    constructor%varName(33) = 'inorganic_phosphorus'
    constructor%varName(34) = 'DON_NonRefractory'
    constructor%varName(35) = 'DOP_NonRefractory'
    constructor%varName(36) = 'DON_Refractory'
    constructor%varName(37) = 'DOP_Refractory'
    constructor%varName(38) = 'partOrgNit'
    constructor%varName(39) = 'partOrgPho'
    constructor%varName(40) = 'phytoplankton'
    constructor%varName(41) = 'zooplankton'
    constructor%varName(42) = 'initial_volume'
    constructor%varName(43) = 'temp'
    constructor%varName(44) = 'salt'
    constructor%varName(45) = 'radius_cr_min'
    constructor%varName(46) = 'radius_cr_max' 
	
    end function constructor

    end module tracerwaterQuality_mod
