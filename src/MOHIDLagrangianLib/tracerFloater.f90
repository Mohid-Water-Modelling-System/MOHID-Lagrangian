    !------------------------------------------------------------------------------
    !        Colab+Atlantic, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_floater
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : Colab+Atlantic, Marine Modelling Group
    ! DATE          : May 2026
    !> @author
    !> Mohsen Shabani Email:shabani.mohsen@outlook.com	
    !
    ! DESCRIPTION:
    !> Module that defines a Lagrangian tracer class for floater modelling and related methods.
    !> The type is defined as a derived type from the pule Lagrangian tracer, and hence inherits all
    !> of it's data and methods
    !------------------------------------------------------------------------------

    module tracerFloater_mod

    use tracerBase_mod
    use common_modules
    use sources_mod

    implicit none
    private

    type :: floater_par_class               !<Type - parameters of a Lagrangian tracer object representing a floater material
        integer    :: particulate                   !< flag to indicate if the material is a particle (false) or a collection of particles (true)
        real(prec) :: size                          !< Size (radius) of the particles (equals to the tracer radius if particulate==false)
        real(prec)  :: Mass
		
        integer  	:: Sphere_N
        real(prec)  :: Sphere_ratio_dry
        real(prec)  :: Sphere_CD
        real(prec)  :: Sphere_Radius
		
        integer  	:: Rectangular_N
        real(prec)  :: Rectangular_ratio_dry
        real(prec)  :: Rectangular_CD
        real(prec)  :: Rectangular_Lx
        real(prec)  :: Rectangular_Ly
        real(prec)  :: Rectangular_Lz
		
        integer  	:: Porous_N
        real(prec)  :: Porous_ratio_dry
        real(prec)  :: Porous_CD
        real(prec)  :: Porous_Lx
        real(prec)  :: Porous_Ly
        real(prec)  :: Porous_Lz
        real(prec)  :: Porous_Depth		
		
        real(prec)  :: Initial_ufloater
        real(prec)  :: Initial_vfloater
        real(prec)  :: Initial_wfloater
    end type floater_par_class

    type :: floater_state_class             !<Type - State variables of a tracer object representing a floater material

        real(prec)  :: temperature                   !< temperature of the tracer
        real(prec)  :: salinity                      !< salinity of the tracer
		real(prec)  :: ufloater						!< x velocity of the tracer
		real(prec)  :: vfloater						!< y velocity of the tracer
		real(prec)  :: wfloater						!< z velocity of the tracer
		
    end type floater_state_class

    type, extends(tracer_class) :: floater_class    !<Type - The floater material Lagrangian tracer class
        type(floater_par_class)   :: mpar     !<To access material parameters
        type(floater_state_class) :: mnow     !<To access material state variables
    contains
    procedure :: getNumVars
    procedure :: getStateArray
    procedure :: setStateArray
    end type floater_class

    !Public access vars
    public :: floater_class

    !Public access routines
    public :: floaterTracer

    interface floaterTracer !< Constructor
    procedure constructor
    end interface

    contains

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com		
    !> @brief
    !> Method that returns the number of variables used by this tracer
    !---------------------------------------------------------------------------
    integer function getNumVars(self)
    class(floater_class), intent(in) :: self
    getNumVars = 27
    end function getNumVars

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Method that returns the state array of this tracer
    !---------------------------------------------------------------------------
    function getStateArray(self)
    class(floater_class), intent(in) :: self
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
    getStateArray(23) = self%mnow%temperature
    getStateArray(24) = self%mnow%salinity
    getStateArray(25) = self%mnow%ufloater
    getStateArray(26) = self%mnow%vfloater
    getStateArray(27) = self%mnow%wfloater
    end function getStateArray

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com		
    !> @brief
    !> Method that sets the state array of this tracer
    !---------------------------------------------------------------------------
    subroutine setStateArray(self, stateArray)
    class(floater_class), intent(inout) :: self
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
    self%mpar%particulate			= StateArray(14)
    self%now%bathymetry   			= StateArray(15)
    self%now%rugosityVar   			= StateArray(16)
    self%now%D50Var   				= StateArray(17)
    self%now%dwz          			= StateArray(18)
    self%now%dist2bottom 			= StateArray(19)
    self%now%beachPeriod 			= StateArray(20)
    self%now%beachAreaId 			= StateArray(21)
    self%now%beachedWaterLevel		= StateArray(22)
    self%mnow%temperature 			= StateArray(23)
    self%mnow%salinity 				= StateArray(24)
    self%mnow%ufloater 				= StateArray(25)
    self%mnow%vfloater 				= StateArray(26)
    self%mnow%wfloater 				= StateArray(27)
    end subroutine setStateArray

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com		
    !> @brief
    !> Floater Tracer constructor
    !> @param[in] id, src, time, p
    !---------------------------------------------------------------------------
    function constructor(id, src, time, p)
    type(floater_class) :: constructor
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
    constructor%par%ttype = Globals%Types%floater
    constructor%mpar%particulate = src%prop%particulate
    ! constructor%mpar%size = src%prop%radius
    ! !material state
    ! constructor%mnow%density = src%prop%density
    ! constructor%mnow%radius = src%prop%radius
    ! constructor%mnow%volume = src%prop%volume
    ! constructor%mnow%area = src%prop%area
    
	!default values
    ! constructor%mnow%condition = 1.0
    ! constructor%mnow%degradation_rate = 1/(100*365*24*3600)
    ! constructor%mnow%concentration = 1000000						! TODO: the value should comaptible with a correct value!  
	! constructor%mnow%Initial_volume = src%prop%volume
    constructor%mnow%temperature = 15.0								! TODO: the value should comaptible with a correct value!
    constructor%mnow%salinity = 36.0								! TODO: the value should comaptible with a correct value!
	
    ! constructor%mnow%radius_cr_min = 1.0e-4_prec * src%prop%radius	! TODO: the value should comaptible with a correct value!
    ! constructor%mnow%radius_cr_max = 1.0e+4_prec * src%prop%radius	! TODO: the value should comaptible with a correct value!
    
	constructor%mpar%Mass = 1.0
	constructor%mpar%Sphere_N = 1
	constructor%mpar%Sphere_ratio_dry = 0.5
	constructor%mpar%Sphere_CD = MV
	constructor%mpar%Sphere_Radius = 0.5

	constructor%mpar%Rectangular_N = 1
	constructor%mpar%Rectangular_ratio_dry = 0.5
	constructor%mpar%Rectangular_CD = MV
	constructor%mpar%Rectangular_Lx = 1.0
	constructor%mpar%Rectangular_Ly = 1.0
	constructor%mpar%Rectangular_Lz = 1.0

	constructor%mpar%Porous_N = 1
	constructor%mpar%Porous_ratio_dry = 1.0
	constructor%mpar%Porous_CD = MV
	constructor%mpar%Porous_Lx = 1.0
	constructor%mpar%Porous_Ly = 1.0
	constructor%mpar%Porous_Lz = 1.0  
	constructor%mpar%Porous_Depth = 5.0

	constructor%mpar%Initial_ufloater = 0.0
	constructor%mpar%Initial_vfloater = 0.0
	constructor%mpar%Initial_wfloater = 0.0 

    !default values
    constructor%mnow%temperature = 15.0
    constructor%mnow%salinity = 36.0
	

    
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

	tag = 'Mass'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Mass = src%prop%propValue(idx)
	end if

	tag = 'Sphere_N'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Sphere_N = src%prop%propValue(idx)
	end if

	tag = 'Sphere_ratio_dry'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Sphere_ratio_dry = src%prop%propValue(idx)
	end if

	tag = 'Sphere_CD'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Sphere_CD = src%prop%propValue(idx)
	end if

	tag = 'Sphere_Radius'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Sphere_Radius = src%prop%propValue(idx)
	end if

	tag = 'Rectangular_N'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Rectangular_N = src%prop%propValue(idx)
	end if

	tag = 'Rectangular_ratio_dry'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Rectangular_ratio_dry = src%prop%propValue(idx)
	end if

	tag = 'Rectangular_CD'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Rectangular_CD = src%prop%propValue(idx)
	end if

	tag = 'Rectangular_Lx'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Rectangular_Lx = src%prop%propValue(idx)
	end if

	tag = 'Rectangular_Ly'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Rectangular_Ly = src%prop%propValue(idx)
	end if

	tag = 'Rectangular_Lz'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Rectangular_Lz = src%prop%propValue(idx)
	end if

	tag = 'Porous_N'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_N = src%prop%propValue(idx)
	end if

	tag = 'Porous_ratio_dry'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_ratio_dry = src%prop%propValue(idx)
	end if

	tag = 'Porous_CD'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_CD = src%prop%propValue(idx)
	end if

	tag = 'Porous_Lx'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_Lx = src%prop%propValue(idx)
	end if

	tag = 'Porous_Ly'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_Ly = src%prop%propValue(idx)
	end if

	tag = 'Porous_Lz'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_Lz = src%prop%propValue(idx)
	end if

	tag = 'Porous_Depth'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Porous_Depth = src%prop%propValue(idx)
	end if
	
	tag = 'Initial_ufloater'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Initial_ufloater = src%prop%propValue(idx)
	end if

	tag = 'Initial_vfloater'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Initial_vfloater = src%prop%propValue(idx)
	end if

	tag = 'Initial_wfloater'
	idx = Utils%find_str(src%prop%propName, tag, .false.)
	if (idx /= MV_INT) then
		constructor%mpar%Initial_wfloater = src%prop%propValue(idx)
	end if

    if (constructor%mpar%particulate==1) then
        !constructor%mpar%size = src%prop%pt_radius !correcting size to now mean particle size, not tracer size
        !constructor%mnow%concentration = src%prop%ini_concentration
    end if

	constructor%mnow%ufloater = constructor%mpar%Initial_ufloater
	constructor%mnow%vfloater = constructor%mpar%Initial_vfloater
	constructor%mnow%wfloater = constructor%mpar%Initial_wfloater

    !filling the rest of the varName list
    constructor%varName(23) = 'temp'
    constructor%varName(24) = 'salt'
    constructor%varName(25) = 'ufloater'
    constructor%varName(26) = 'vfloater'
    constructor%varName(27) = 'wfloater'

    
    end function constructor

    end module tracerFloater_mod
