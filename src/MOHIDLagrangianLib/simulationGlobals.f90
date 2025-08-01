    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_globals
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the simulation global parameter classes and their methods
    !------------------------------------------------------------------------------

    module simulationGlobals_mod

    use penf
    use vecfor_r8p
    !use vecfor_r4p
    use stringifor
    use datetime_module
    use FoX_dom

    use simulationPrecision_mod
    use simulationParallel_omp_mod
    use simulationLogger_mod
    use simulationMemory_mod
    use utilities_mod
    use xmlParser_mod
    use abstract_LinkedList_mod
    use geometry_mod
    use simulationMemory_mod

    implicit none
    private
    
    type, extends(linkedlist) :: stringList_class !< List of strings class
    contains
    procedure :: print => print_stringList
    procedure :: printCurrent => print_stringListCurrent
    procedure :: notRepeated
    procedure :: makeUnion
    procedure :: toArray
    procedure :: copy
    end type stringList_class

    type :: parameters_t   !< Parameters class
        integer         :: Integrator = 1            !< Integration Algorithm 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)
        integer         :: IntegratorIndexes(3)      !< Index list for the integrator selector
        type(string)    :: IntegratorNames(3)        !< Names list for the integrator selector
        integer         :: numOPMthreads             !< number of openMP threads to be used
        real(prec)      :: WarmUpTime                !< Time to freeze the tracers at simulation start (warmup) (s) (default=0.0)
        real(prec)      :: TimeMax                   !< Simulation duration (s)
        real(prec)      :: OutputWriteTime           !< Output write time (1/Hz)
        type(datetime)  :: StartTime                 !< Start date of the simulation
        type(datetime)  :: EndTime                   !< End date of the simulation
        type(datetime)  :: BaseDateTime              !< Base date for time stamping results
        real(prec)      :: BufferSize                !< Controls the frequency of consumption and ammout of input data kept in memory
        integer         :: OutputFormat              !< Format of the output files (default=2) NetCDF=1, VTK=2
        integer         :: OutputFormatIndexes(2)    !< Index list for the output file format selector
        type(string)    :: OutputFormatNames(2)      !< Names list for the output file format selector
        real(prec)      :: FillValueReal = -9.9e15
    contains
    procedure :: setParam
    procedure :: check
    procedure :: print => printsimparameters
    end type parameters_t

    type :: simdefs_t  !< Simulation definitions class
        type(vector)    ::  Dp              !< Initial particle spacing at emission
        real(prec)      ::  dt = MV         !< Timestep for fixed step integrators (s)
        type(vector)    ::  Pointmin        !< Point that defines the lowest corner of the simulation bounding box
        type(vector)    ::  Pointmax        !< Point that defines the upper corner of the simulation bounding box
        type(vector)    ::  Center          !< Point that defines the center of the simulation bounding box
        logical         ::  autoblocksize   !< Flag for automatic Block sizing
        type(vector)    ::  blocksize       !< Size (xyz) of a Block (sub-domain)
        integer         ::  numblocks       !< Number of blocks in the simulation
        integer         ::  numblocksx, numblocksy  !<Number of blocks along x and y
        integer         ::  VerticalVelMethod !< Vertical velocity method
        integer         ::  RemoveLandTracer !< Vertical velocity method
        integer         ::  bathyminNetcdf !< bathymetry is a property inside the netcdf
        real(prec)      ::  tracerMaxAge !< removes tracers with age greater than maxTracerAge
        real(prec)      ::  Temperature_add_offset !< adds offset to temperature from netcdf
        real(prec)      ::  WqDt !< water quality process time step
        logical         ::  inputFromHDF5 = .false.
        integer         ::  FreeLitterAtBeaching = 0
    contains
    procedure :: setdp
    procedure :: setdt
    procedure :: setboundingbox
    procedure :: setCenter
    procedure :: setblocksize
    procedure :: setVerticalVelMethod
    procedure :: setRemoveLandTracer
    procedure :: setbathyminNetcdf
    procedure :: settracerMaxAge
    procedure :: setTemperature_add_offset  !Only here becasue some files from SINTEF were not properly converted due to
                                            ! some error in mohid's convert to hdf5 tool
    procedure :: setWqDt
    procedure :: setFreeLitterAtBeaching
    procedure, public :: setInputFileType
    procedure :: print => printsimdefs
    end type simdefs_t

    type :: constants_t    !< Case Constants class
        type(vector) :: Gravity             !< Gravitational acceleration vector (default=(0 0 -9.81)) (m s-2)
        real(prec)   :: Z0 = 0.0            !< Reference local sea level
        real(prec)   :: smallDt             !< Small dt scale, for numeric precision purposes
        real(prec)   :: BeachingLevel       !<Level above which beaching can occur (m)
        real(prec)   :: ResuspensionCoeff   !< Resuspension velocity amplitud factor (default=0) 
        real(prec)   :: BeachingStopProb    !< Probablity of beaching stopping a tracer (-)
        real(prec)   :: DiffusionCoeff      !< Horizontal diffusion coefficient (-)
        real(prec)   :: WindDragCoeff      !< Wind drag coefficient (% of wind force)
        real(prec)   :: MeanDensity         !< mean ocean water density.
        real(prec)   :: MeanKVisco          !< mean ocean water kinematic viscosity
        integer      :: AddBottomCell       !< open the first bottom cell and put the same value as the cell above
        real(prec)   :: Rugosity            !< Bottom rugosity value (for now the model assumes a constant value over the grid)
        real(prec)   :: Critical_Shear_Erosion !< Bottom critical shear erosion value (for now the model assumes a constant value over the grid)
        real(prec)   :: TOptBacteriaMin          !< minimum temperature of the optimal interval for the Bacteria growth
        real(prec)   :: TOptBacteriaMax          !< maximum temperature of the optimal interval for the Bacteria growt
        real(prec)   :: TBacteriaMin          !< minimum tolerable temperature of the interval for the Bacteria growth
        real(prec)   :: TBacteriaMax          !< maximum tolerable temperature of the interval for the Bacteria growth
        real(prec)   :: BK1          !< constant to control temperature response curve shape
        real(prec)   :: BK2          !< constant to control temperature response curve shape
        real(prec)   :: BK3          !< constant to control temperature response curve shape
        real(prec)   :: BK4          !< constant to control temperature response curve shape
        real(prec)   :: MaxDegradationRate !< maximum degradation rate from bacteria
    contains
    procedure :: setgravity
    procedure :: setz0
    procedure :: setBeachingLevel
    procedure :: setBeachingStopProb
    procedure :: setDiffusionCoeff
    procedure :: setWindDragCoeff
    procedure :: setSmallDt
    procedure :: setResuspensionCoeff
    procedure :: setMeanDensity
    procedure :: setMeanKVisco
    procedure :: setAddBottomCell
    procedure :: setRugosity
    procedure :: setCritical_Shear_Erosion
    procedure :: setTOptBacteriaMin
    procedure :: setTOptBacteriaMax
    procedure :: setTBacteriaMin
    procedure :: setTBacteriaMax
    procedure :: setBK1
    procedure :: setBK2
    procedure :: setBK3
    procedure :: setBK4
    procedure :: setMaxDegradationRate
    procedure :: print => printconstants
    end type constants_t

    type :: filenames_t    !<File names class
        type(string) :: mainxmlfilename     !< Input .xml file name
        type(string) :: propsxmlfilename    !< Properties .xml file name
        type(string) :: inputsXmlFilename   !< .xml file name with metadata on input files
        type(string), allocatable, dimension(:) :: inputFile   !< File names of input data
        type(string), allocatable, dimension(:) :: namingfilename   !< File names of the naming convention .xml files
        type(string) :: tempfilename        !< Generic temporary file name
        type(string) :: outpath             !< General output directory
        type(string) :: casename            !< Name of the running case
    end type filenames_t

    type :: sim_t  !<Simulation related counters and others
        private
        integer :: numdt        !<number of the current iteration
        integer :: numTracer    !<Global Tracer number holder. Incremented at tracer construction or first activation time
    contains
    procedure, public :: getnumdt
    procedure, public :: increment_numdt
    procedure, public :: getnumTracer
    procedure, private :: increment_numTracer
    end type sim_t

    type :: output_t  !<Simulation related counters and others
        private
        integer :: lastOutNumDt     !<number of the last outputed iteration
        integer :: numOutFile       !<number of the current output file
        type(stringList_class) :: varToOutput
    contains
    procedure, public :: increment_numOutFile
    procedure, public :: getnumOutFile
    procedure, public :: getlastOutNumDt
    procedure, public :: setlastOutNumDt
    procedure, private :: addToOutputPool
    procedure, public :: getOutputPoolArray
    procedure, public :: setOutputFields
    end type output_t

    type :: var_names_t
        type(string) :: u !< Name of the 'u' variable in the model
        type(string) :: v
        type(string) :: w
        type(string) :: ssh
        type(string) :: temp
        type(string) :: sal
        type(string) :: density
        type(string) :: vsdx
        type(string) :: vsdy
        type(string) :: hs
        type(string) :: ts
        type(string) :: wd
        type(string) :: wl
        type(string) :: u10
        type(string) :: v10
        type(string) :: rad
        type(string) :: dissolved_oxygen
        type(string) :: nitrate
        type(string) :: nitrite
        type(string) :: ammonia
        type(string) :: partOrgNit
        type(string) :: partOrgPho
        type(string) :: DON_NonRefractory
        type(string) :: DOP_NonRefractory
        type(string) :: DON_Refractory
        type(string) :: DOP_Refractory
        type(string) :: phytoplankton
        type(string) :: zooplankton
        type(string) :: inorganic_phosphorus
        type(string) :: openpoints
        type(string) :: lon
        type(string) :: lat
        type(string) :: level
        type(string) :: time
        type(string) :: landIntMask
        type(string) :: resolution
        type(string) :: bathymetry
        type(string) :: surface
        type(string) :: rate
        type(string) :: dwz
        type(stringList_class) :: bathymetryVariants
        type(stringList_class) :: uVariants !< possible names for 'u' in the input files
        type(stringList_class) :: vVariants
        type(stringList_class) :: wVariants
        type(stringList_class) :: sshVariants
        type(stringList_class) :: tempVariants
        type(stringList_class) :: salVariants
        type(stringList_class) :: densityVariants
        type(stringList_class) :: vsdxVariants
        type(stringList_class) :: vsdyVariants
        type(stringList_class) :: hsVariants
        type(stringList_class) :: tsVariants
        type(stringList_class) :: wdVariants
        type(stringList_class) :: wlVariants
        type(stringList_class) :: u10Variants
        type(stringList_class) :: v10Variants
        type(stringList_class) :: radVariants
        type(stringList_class) :: dissolved_oxygenVariants
        type(stringList_class) :: nitrateVariants
        type(stringList_class) :: nitriteVariants
        type(stringList_class) :: ammoniaVariants
        type(stringList_class) :: partOrgNitVariants
        type(stringList_class) :: partOrgPhoVariants
        type(stringList_class) :: DON_NonRefractoryVariants
        type(stringList_class) :: DOP_NonRefractoryVariants
        type(stringList_class) :: DON_RefractoryVariants
        type(stringList_class) :: DOP_RefractoryVariants
        type(stringList_class) :: phytoplanktonVariants
        type(stringList_class) :: zooplanktonVariants
        type(stringList_class) :: inorganic_phosphorusVariants
        type(stringList_class) :: openpointsVariants
        type(stringList_class) :: lonVariants
        type(stringList_class) :: latVariants
        type(stringList_class) :: levelVariants
        type(stringList_class) :: timeVariants
        type(stringList_class) :: rateVariants
        type(stringList_class) :: varToUse !< list of variables used in a simulation
    contains
    procedure, private :: buildvars
    procedure, public  :: addVar
    procedure, public  :: getVarSimName
    procedure, public  :: checkVarSimName
    procedure, public  :: checkDimensionName
    procedure, public  :: checkSpaceDimensionName
    procedure, public  :: checkLatOrLon
    procedure, public  :: checkVarIsLon
    procedure, public  :: checkVarIsLat
    procedure, public  :: checkVarIsDepth
    procedure, public  :: checkVarIsTime
    procedure, public  :: checkVarIsSSH
    procedure, public  :: checkMappingVar
    end type var_names_t
    
    type :: varBackground_t
        integer :: bkgIndex = 0
        type(stringList_class) :: bkgVars 
    contains
    end type varBackground_t
    
    type :: varBackgroundDict_t
        type(varBackground_t), allocatable, dimension(:) :: bkgDict 
    contains
    end type varBackgroundDict_t

    type :: maskVals_t
        real(prec) :: landVal  = 2.0
        real(prec) :: bedVal   = -1.0
        real(prec) :: waterVal = 0.0
        real(prec) :: beachVal = 1.0
    contains
    end type maskVals_t

    type :: tracerTypes_t
        integer :: base  = 0
        integer :: paper   = 1
        integer :: plastic = 2
        integer :: coliform = 3
        integer :: seed = 4
        integer :: detritus = 5
        integer :: waterQuality = 6
    contains
    end type tracerTypes_t

    type :: dataTypes_t
        type(string) :: currents
        type(string) :: winds
        type(string) :: waves
        type(string) :: waterProps
    contains
    end type dataTypes_t

    type :: sim_time_t
        type(datetime)  :: BaseDateTime              !< Base date for time stamping results
        real(prec)      :: TimeMax = MV              !< Simulation duration (s)
        type(datetime)  :: StartDate                 !< Start date of the simulation
        type(datetime)  :: EndDate                   !< End date of the simulation
        type(datetime)  :: CurrDate                  !< Current date of the simulation
        real(prec)      :: CurrTime = 0              !< Current time of the simulation (s)
    contains
    procedure :: print => printDateTime
    procedure :: setCurrDateTime
    procedure :: getDateTimeStamp
    end type sim_time_t
    
    type :: sources_t
        integer, allocatable, dimension(:) :: sourcesID
        real(prec), allocatable, dimension(:) :: bottom_emission_depth
        real(prec), allocatable, dimension(:) :: biofouling_rate
        real(prec), allocatable, dimension(:) :: biofouling_start_after
    contains    
    end type sources_t
    
    type :: extImpFiles_t
        logical :: waterQualityFileName_hasValue = .false.
        type(string) :: waterQualityFileName
    contains
    procedure :: setWaterQualityFileName
    procedure :: print => printWaterQualityFileName !Will need changing when other files are added, such as oil.
    end type extImpFiles_t
    
    type :: beach_par                          !<Type - parameters of a source object
        integer :: id                           !< unique source identification (integer)
        integer :: coastType                    !< Coast type used to define the probability method
        integer :: unbeach                      !< Should tracers unbeach (1 or 0)
        integer :: runUpEffect                  !< Compute Runup effect
        integer :: runUpEffectUnbeach           !< include runup effect for unbeaching
        real(prec) :: probability               !< adjustable probability of beaching
        real(prec) :: waterColumnThreshold      !< minimum water column for beaching to be computed
        real(prec) :: beachTimeScale            !< Time scale use to calculate the beach probability = 1-exp(dt/beach_time_scale). dt = (CurrentTime - LastAtualization)
        real(prec) :: unbeachTimeScale          !< Time scale use to calculate the Unbeach probability = 1-exp(BeachPeriod/unbeach_time_scale)
        real(prec) :: beachSlope                !< Beach slope
        type(string) :: beaching_geometry       !< polygon describing the beach geometry
        class(shape), allocatable :: geometry   !< geometry shapetype (must be polygon)
    end type beach_par
    
    type :: beach_class           !<Type - The beach class
        type(beach_par)   :: par           !<To access parameters
    contains
    procedure :: initialize => initializeBeachAreas
    !procedure :: isParticulate !Sobrinho - Apagar
    !procedure :: check !Sobrinho - Apagar
    end type beach_class

    type :: beachArray_class
        type(beach_class), allocatable, dimension(:) :: beachArea
    contains
    procedure :: initialize => initAreas
    procedure :: finalize => killAreas
    end type beachArray_class

    !Simulation variables
    
    type :: globals_class   !<Globals class - This is a container for every global variable on the simulation
        type(parameters_t)  :: Parameters
        type(simdefs_t)     :: SimDefs
        type(constants_t)   :: Constants
        type(filenames_t)   :: Names
        type(sim_t)         :: Sim
        type(varBackgroundDict_t) :: BackgroundVarDict
        type(output_t)      :: Output
        type(var_names_t)   :: Var
        type(sim_time_t)    :: SimTime
        type(maskVals_t)    :: Mask
        type(tracerTypes_t) :: Types
        type(dataTypes_t)   :: DataTypes
        type(sources_t)     :: Sources
        type(extImpFiles_t) :: ExtImpFiles !Will have to be changed when other types of files are added (oil for example)
        type(beachArray_class) :: BeachingAreas !< Beaching Areas array, used exclusively for building the case from a description file
    contains
    procedure :: initialize => setdefaults
    procedure :: setTimeDate
    procedure, public :: setNamingConventions
    procedure, public :: setInputFileNames
    procedure :: setVarNames
    procedure :: setDimNames
    procedure :: setCurrVar
    procedure :: fillBackgroundDict
    end type globals_class
    
    type(string) :: notRead
    type(string) :: notSet

    !Simulation variables
    type(globals_class) :: Globals

    !Public access vars
    public :: Globals, stringList_class, notRead, notSet, varBackground_t

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Globals default setting routine.
    !> @param[in] self, outpath
    !---------------------------------------------------------------------------
    subroutine setdefaults(self, outpath)
    class(globals_class), intent(inout) :: self
    integer :: sizem
    type(string), optional, intent(in) :: outpath
    
    notRead = 'notRead'
    notSet = 'notSet'
    
    !parameters
    self%Parameters%Integrator = 1
    self%Parameters%IntegratorIndexes = [1,2,3]
    self%Parameters%IntegratorNames(1) = 'Euler'
    self%Parameters%IntegratorNames(2) = 'Multi-Step Euler'
    self%Parameters%IntegratorNames(3) = 'Runge-Kuta 4'
    self%Parameters%numOPMthreads = OMPManager%getThreads()
    self%Parameters%WarmUpTime = 0.0
    self%Parameters%OutputWriteTime = MV
    self%Parameters%StartTime = datetime()
    self%Parameters%EndTime = datetime()
    self%Parameters%BaseDateTime = datetime(1950,1,1,0,0,0)
    self%Parameters%BufferSize = 3600*24*6 !seconds/hour*hours*days
    self%Parameters%OutputFormat = 2
    self%Parameters%OutputFormatIndexes = [1,2]
    !self%Parameters%OutputFormatNames = ['NetCDF','VTK'] !This is not acceptable because FORTRAN
    self%Parameters%OutputFormatNames(1) = 'NetCDF'
    self%Parameters%OutputFormatNames(2) = 'VTK'
    !SimTime
    self%SimTime%StartDate = datetime()
    self%SimTime%EndDate = datetime()
    self%SimTime%CurrDate = datetime()
    self%SimTime%BaseDateTime = datetime(1950,1,1,0,0,0)
    !Simulation definitions
    self%SimDefs%autoblocksize =.true.
    self%SimDefs%blocksize = 0.0
    self%SimDefs%numblocksx = MV
    self%SimDefs%numblocksy = MV
    self%SimDefs%numblocks = OMPManager%getThreads()
    self%SimDefs%Dp = 0.0
    self%SimDefs%dt = MV
    self%SimDefs%Pointmin = 0.0
    self%SimDefs%Pointmax = 0.0
    self%SimDefs%Center = 0.0
    self%SimDefs%VerticalVelMethod = 1
    self%SimDefs%RemoveLandTracer = 0
    self%SimDefs%bathyminNetcdf = 0
    self%SimDefs%tracerMaxAge = 0
    self%SimDefs%Temperature_add_offset = 0
    !simulation constants
    self%Constants%Gravity= 0.0*ex + 0.0*ey -9.81*ez
    self%Constants%Z0 = 0.0
    self%Constants%BeachingLevel = -3.0
    self%Constants%BeachingStopProb = 0.50
    self%Constants%DiffusionCoeff = 1.0
    self%Constants%WindDragCoeff = 0.03
    self%Constants%smallDt = 0.0
    self%Constants%ResuspensionCoeff = 0.0
    self%Constants%MeanDensity = 1027.0
    self%Constants%MeanKVisco = 1.09E-3
    self%Constants%AddBottomCell = 0
    self%Constants%Rugosity = 0.0025
    self%Constants%Critical_Shear_Erosion = 0.4
    self%Constants%TOptBacteriaMin = 24.8
    self%Constants%TOptBacteriaMax = 25.1
    self%Constants%TBacteriaMin = 5
    self%Constants%TBacteriaMax = 35
    self%Constants%BK1 = 0.05
    self%Constants%BK2 = 0.98
    self%Constants%BK3 = 0.98
    self%Constants%BK4 = 0.02
    self%Constants%MaxDegradationRate = 0.03/86400
    !filenames
    self%Names%mainxmlfilename = notSet
    self%Names%propsxmlfilename = notSet
    self%Names%tempfilename = notSet
    if (present(outpath)) then
        self%Names%outpath = outpath
    else
        self%Names%outpath = notSet
    end if
    self%Names%casename = notSet
    !global time
    self%SimTime%CurrTime = 0.0
    !global counters
    self%Sim%numdt = 0
    self%Sim%numTracer = 0
    !output control variables
    self%Output%lastOutNumDt = 0
    self%Output%numoutfile = 0
    !data types
    self%DataTypes%currents = 'hydrodynamic'
    self%DataTypes%winds = 'meteorology'
    self%DataTypes%waves = 'waves'
    self%DataTypes%waterProps = 'waterProperties'
    !Variable names
    call self%Var%buildvars()

    sizem=sizeof(self)
    call SimMemory%adddef(sizem)

    end subroutine setdefaults
    
    subroutine setInputFileType (self, isInputFileHDF5)
    class(simdefs_t), intent(inout) :: self
    logical isInputFileHDF5
    !Begin------------------------------------------
    
    self%inputFromHDF5 = isInputFileHDF5
        
    end subroutine setInputFileType

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Builds variable list names.
    !---------------------------------------------------------------------------
    subroutine buildvars(self)
    class(var_names_t), intent(inout) :: self
    self%u       = 'u'
    self%v       = 'v'
    self%w       = 'w'
    self%ssh     = 'ssh'
    self%temp    = 'temp'
    self%sal     = 'sal'
    self%density = 'density'
    self%vsdx    = 'vsdx'
    self%vsdy    = 'vsdy'
    self%hs      = 'hs'
    self%ts      = 'ts'
    self%wd      = 'wd'
    self%wd      = 'wl'
    self%u10     = 'u10'
    self%v10     = 'v10'
    self%rad     = 'rad'
    self%lon     = 'lon'
    self%lat     = 'lat'
    self%level   = 'level'
    self%time    = 'time'
    self%landIntMask = 'landIntMask'
    self%resolution = 'resolution'
    self%bathymetry = 'bathymetry'
    self%surface    = 'surface'
    self%rate = 'rate'
    !DWZ is the distance between two vertical faces of a cube (ex: thichness of the vertical layer)
    self%dwz = 'dwz'
    !adding variables to variable pool - PLACEHOLDER, this should come from tracer constructors
    call self%addVar(self%u)
    call self%addVar(self%v)
    call self%addVar(self%w)
    call self%addVar(self%ssh)
    call self%addVar(self%temp)
    call self%addVar(self%sal)
    call self%addVar(self%density)
    call self%addVar(self%bathymetry)
    !call self%addVar(self%lon)
    !call self%addVar(self%lat)
    !call self%addVar(self%level)
    !call self%addVar(self%time)
    end subroutine buildvars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> adding variables to variable pool of a simulation
    !---------------------------------------------------------------------------
    subroutine addVar(self, var)
    class(var_names_t), intent(inout) :: self
    type(string), intent(in) :: var
    if (self%varToUse%notRepeated(var)) call self%varToUse%add(var)
    end subroutine addVar

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns the simulation name of a given variabe name, by brute force searching
    !> trough the naming variable lists...
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    type(string) function getVarSimName(self, var)
    class(var_names_t), intent(inout) :: self
    type(string), intent(in) :: var
    logical value

    !searching for bathymetry
    value = self%bathymetryVariants%notRepeated(var)
    if (var == self%bathymetry .or. .not. value) then
        getVarSimName = self%bathymetry
        return
    end if
    
    !searching for u
    if (var == self%u .or. .not.self%uVariants%notRepeated(var)) then
        getVarSimName = self%u
        return
    end if
    !searching for v
    if (var == self%v .or. .not.self%vVariants%notRepeated(var)) then
        getVarSimName = self%v
        return
    end if
    !searching for w
    if (var == self%w .or. .not.self%wVariants%notRepeated(var)) then
        getVarSimName = self%w
        return
    end if
    !searching for ssh
    if (var == self%ssh .or. .not.self%sshVariants%notRepeated(var)) then
        getVarSimName = self%ssh
        return
    end if
    !searching for temp
    if (var == self%temp .or. .not.self%tempVariants%notRepeated(var)) then
        getVarSimName = self%temp
        return
    end if
    !searching for sal
    if (var == self%sal .or. .not.self%salVariants%notRepeated(var)) then
        getVarSimName = self%sal
        return
    end if
    !searching for density
    if (var == self%density .or. .not.self%densityVariants%notRepeated(var)) then
        getVarSimName = self%density
        return
    end if
    !searching for vsdx
    if (var == self%vsdx .or. .not.self%vsdxVariants%notRepeated(var)) then
        getVarSimName = self%vsdx
        return
    end if
    !searching for vsdy
    if (var == self%vsdy .or. .not.self%vsdyVariants%notRepeated(var)) then
        getVarSimName = self%vsdy
        return
    end if
    !searching for hs
    if (var == self%hs .or. .not.self%hsVariants%notRepeated(var)) then
        getVarSimName = self%hs
        return
    end if
    !searching for ts
    if (var == self%ts .or. .not.self%tsVariants%notRepeated(var)) then
        getVarSimName = self%ts
        return
    end if
    !searching for wd
    if (var == self%wd .or. .not.self%wdVariants%notRepeated(var)) then
        getVarSimName = self%wd
        return
    end if
    !searching for wl
    if (var == self%wl .or. .not.self%wlVariants%notRepeated(var)) then
        getVarSimName = self%wl
        return
    end if
    !searching for u10
    if (var == self%u10 .or. .not.self%u10Variants%notRepeated(var)) then
        getVarSimName = self%u10
        return
    end if
    !searching for v10
    if (var == self%v10 .or. .not.self%v10Variants%notRepeated(var)) then
        getVarSimName = self%v10
        return
    end if
    !searching for rad
    if (var == self%rad .or. .not.self%radVariants%notRepeated(var)) then
        getVarSimName = self%rad
        return
    end if
    !searching for dissolved_oxygen
    if (var == self%dissolved_oxygen .or. .not.self%dissolved_oxygenVariants%notRepeated(var)) then
        getVarSimName = self%dissolved_oxygen
        return
    end if
    !searching for nitrate
    if (var == self%nitrate .or. .not.self%nitrateVariants%notRepeated(var)) then
        getVarSimName = self%nitrate
        return
    end if
    !searching for nitrite
    if (var == self%nitrite .or. .not.self%nitriteVariants%notRepeated(var)) then
        getVarSimName = self%nitrite
        return
    end if
    !searching for ammonia
    if (var == self%ammonia .or. .not.self%ammoniaVariants%notRepeated(var)) then
        getVarSimName = self%ammonia
        return
    end if
    !searching for partOrgNit
    if (var == self%partOrgNit .or. .not.self%partOrgNitVariants%notRepeated(var)) then
        getVarSimName = self%partOrgNit
        return
    end if
    !searching for partOrgPho
    if (var == self%partOrgPho .or. .not.self%partOrgPhoVariants%notRepeated(var)) then
        getVarSimName = self%partOrgPho
        return
    end if
    !searching for DON_NonRefractory
    if (var == self%DON_NonRefractory .or. .not.self%DON_NonRefractoryVariants%notRepeated(var)) then
        getVarSimName = self%DON_NonRefractory
        return
    end if
    !searching for DOP_NonRefractory
    if (var == self%DOP_NonRefractory .or. .not.self%DOP_NonRefractoryVariants%notRepeated(var)) then
        getVarSimName = self%DOP_NonRefractory
        return
    end if
    !searching for DON_Refractory
    if (var == self%DON_Refractory .or. .not.self%DON_RefractoryVariants%notRepeated(var)) then
        getVarSimName = self%DON_Refractory
        return
    end if
    !searching for DOP_Refractory
    if (var == self%DOP_Refractory .or. .not.self%DOP_RefractoryVariants%notRepeated(var)) then
        getVarSimName = self%DOP_Refractory
        return
    end if
    !searching for phytoplankton
    if (var == self%phytoplankton .or. .not.self%phytoplanktonVariants%notRepeated(var)) then
        getVarSimName = self%phytoplankton
        return
    end if
    !searching for zooplankton
    if (var == self%zooplankton .or. .not.self%zooplanktonVariants%notRepeated(var)) then
        getVarSimName = self%zooplankton
        return
    end if
    !searching for inorganic_phosphorus
    if (var == self%inorganic_phosphorus .or. .not.self%inorganic_phosphorusVariants%notRepeated(var)) then
        getVarSimName = self%inorganic_phosphorus
        return
    end if
    !searching for openpoints
    if (var == self%openpoints .or. .not.self%openpointsVariants%notRepeated(var)) then
        getVarSimName = self%openpoints
        return
    end if
    !searching for lon
    if (var == self%lon .or. .not.self%lonVariants%notRepeated(var)) then
        getVarSimName = self%lon
        return
    end if
    !searching for lat
    if (var == self%lat .or. .not.self%latVariants%notRepeated(var)) then
        getVarSimName = self%lat
        return
    end if
    !searching for level
    if (var == self%level .or. .not.self%levelVariants%notRepeated(var)) then
        getVarSimName = self%level
        return
    end if
    !searching for time
    if (var == self%time .or. .not.self%timeVariants%notRepeated(var)) then
        getVarSimName = self%time
        return
    end if
    !searching for rate
    if (var == self%rate .or. .not.self%rateVariants%notRepeated(var)) then
        getVarSimName = self%rate
        return
    end if

    end function getVarSimName
    
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> returns true is the input variable name is defined in the model
    !> @param[in] var
    !---------------------------------------------------------------------------
    logical function checkVarSimName(self, var)
    class(var_names_t), intent(inout) :: self
    type(string), intent(in) :: var
    logical value
    !Begin-------------------------------------------------------------
    checkVarSimName = .false.
    !searching for bathymetry
    value = self%bathymetryVariants%notRepeated(var)
    if (var == self%bathymetry .or. .not. value) then
        checkVarSimName = .true.
        return
    end if
    
    !searching for u
    if (var == self%u .or. .not.self%uVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for v
    if (var == self%v .or. .not.self%vVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for w
    if (var == self%w .or. .not.self%wVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for ssh
    if (var == self%ssh .or. .not.self%sshVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for temp
    if (var == self%temp .or. .not.self%tempVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for sal
    if (var == self%sal .or. .not.self%salVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for density
    if (var == self%density .or. .not.self%densityVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for vsdx
    if (var == self%vsdx .or. .not.self%vsdxVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for vsdy
    if (var == self%vsdy .or. .not.self%vsdyVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for hs
    if (var == self%hs .or. .not.self%hsVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for ts
    if (var == self%ts .or. .not.self%tsVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for wd
    if (var == self%wd .or. .not.self%wdVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for wl
    if (var == self%wl .or. .not.self%wlVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for u10
    if (var == self%u10 .or. .not.self%u10Variants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for v10
    if (var == self%v10 .or. .not.self%v10Variants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for rad
    if (var == self%rad .or. .not.self%radVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for dissolved_oxygen
    if (var == self%dissolved_oxygen .or. .not.self%dissolved_oxygenVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for nitrate
    if (var == self%nitrate .or. .not.self%nitrateVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for nitrite
    if (var == self%nitrite .or. .not.self%nitriteVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for ammonia
    if (var == self%ammonia .or. .not.self%ammoniaVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for partOrgNit
    if (var == self%partOrgNit .or. .not.self%partOrgNitVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for partOrgPho
    if (var == self%partOrgPho .or. .not.self%partOrgPhoVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for DON_NonRefractory
    if (var == self%DON_NonRefractory .or. .not.self%DON_NonRefractoryVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for DOP_NonRefractory
    if (var == self%DOP_NonRefractory .or. .not.self%DOP_NonRefractoryVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for DON_Refractory
    if (var == self%DON_Refractory .or. .not.self%DON_RefractoryVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for DOP_Refractory
    if (var == self%DOP_Refractory .or. .not.self%DOP_RefractoryVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for phytoplankton
    if (var == self%phytoplankton .or. .not.self%phytoplanktonVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for zooplankton
    if (var == self%zooplankton .or. .not.self%zooplanktonVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for inorganic_phosphorus
    if (var == self%inorganic_phosphorus .or. .not.self%inorganic_phosphorusVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for openpoints
    if (var == self%openpoints .or. .not.self%openpointsVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for lon
    if (var == self%lon .or. .not.self%lonVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for lat
    if (var == self%lat .or. .not.self%latVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for level
    if (var == self%level .or. .not.self%levelVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for time
    if (var == self%time .or. .not.self%timeVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if
    !searching for rate
    if (var == self%rate .or. .not.self%rateVariants%notRepeated(var)) then
        checkVarSimName = .true.
        return
    end if

    end function checkVarSimName
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> returns true is the input variable name is defined in the model as a dimension variable
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkDimensionName(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkDimensionName = .false.
    
    !searching for lon
    if (var == self%lon .or. .not.self%lonVariants%notRepeated(var)) then
        checkDimensionName = .true.
        return
    end if
    !searching for lat
    if (var == self%lat .or. .not.self%latVariants%notRepeated(var)) then
        checkDimensionName = .true.
        return
    end if
    !searching for level
    if (var == self%level .or. .not.self%levelVariants%notRepeated(var)) then
        checkDimensionName = .true.
        return
    end if
    !searching for time
    if (var == self%time .or. .not.self%timeVariants%notRepeated(var)) then
        checkDimensionName = .true.
        return
    end if

    end function checkDimensionName
    
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> returns true is the input variable name is defined in the model as a space dimension variable
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkSpaceDimensionName(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkSpaceDimensionName = .false.
    
    !searching for lon
    if (var == self%lon .or. .not.self%lonVariants%notRepeated(var)) then
        checkSpaceDimensionName = .true.
        return
    end if
    !searching for lat
    if (var == self%lat .or. .not.self%latVariants%notRepeated(var)) then
        checkSpaceDimensionName = .true.
        return
    end if
    !searching for level
    if (var == self%level .or. .not.self%levelVariants%notRepeated(var)) then
        checkSpaceDimensionName = .true.
        return
    end if

    end function checkSpaceDimensionName

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> returns true is the input variable name is defined in the model as a dimension variable
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkLatOrLon(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkLatOrLon = .false.
    !searching for lon
    if (var == self%lon .or. .not.self%lonVariants%notRepeated(var)) then
        checkLatOrLon = .true.
        return
    end if
    !searching for lat
    if (var == self%lat .or. .not.self%latVariants%notRepeated(var)) then
        checkLatOrLon = .true.
        return
    end if
    end function checkLatOrLon
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> returns true if the input variable name is defined in the model as longitude
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkVarIsLon(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkVarIsLon = .false.
    !searching for lon
    if (var == self%lon .or. .not.self%lonVariants%notRepeated(var)) then
        checkVarIsLon = .true.
        return
    end if
    end function checkVarIsLon
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> returns true if the input variable name is defined in the model as latitude
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkVarIsLat(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkVarIsLat = .false.
    !searching for lat
    if (var == self%lat .or. .not.self%latVariants%notRepeated(var)) then
        checkVarIsLat = .true.
        return
    end if
    end function checkVarIsLat
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> returns true if the input variable name is defined in the model as depth
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkVarIsDepth(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkVarIsDepth = .false.
    !searching for level
    if (var == self%level .or. .not.self%levelVariants%notRepeated(var)) then
        checkVarIsDepth = .true.
        return
    end if
    end function checkVarIsDepth
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> returns true if the input variable name is defined in the model as time
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkVarIsTime(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkVarIsTime = .false.
    !searching for time
    if (var == self%time .or. .not.self%timeVariants%notRepeated(var)) then
        checkVarIsTime = .true.
        return
    end if
    end function checkVarIsTime
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> returns true if the input variable name is defined in the model as ssh
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkVarIsSSH(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkVarIsSSH = .false.
    !searching for ssh
    if (var == self%ssh .or. .not.self%sshVariants%notRepeated(var)) then
        checkVarIsSSH = .true.
        return
    end if
    end function checkVarIsSSH
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> returns true if the input variable name is defined in the model as openPoints
    !> @param[in] self, var
    !---------------------------------------------------------------------------
    logical function checkMappingVar(self, var)
    type(string), intent(in) :: var
    class(var_names_t), intent(inout) :: self
    !Begin-----------------------------------------------------
    checkMappingVar = .false.
    !searching for time
    if (var == self%openpoints .or. .not.self%openpointsVariants%notRepeated(var)) then
        checkMappingVar = .true.
        return
    end if
    end function checkMappingVar
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> set the name of the input data files
    !---------------------------------------------------------------------------
    subroutine setInputFileNames(self, filenames)
    class(globals_class), intent(inout) :: self
    type(string), dimension(:), intent(in) :: filenames
    allocate(self%Names%inputFile, source = filenames)
    end subroutine setInputFileNames

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> set the naming conventions. Imports names from given .xml
    !> naming files
    !---------------------------------------------------------------------------
    subroutine setNamingConventions(self, filename)
    class(globals_class), intent(inout) :: self
    type(string), dimension(:), intent(in) :: filename
    integer :: i
    type(string) :: outext, tag
    type(Node), pointer :: xmldoc                   !< .xml file handle
    type(Node), pointer :: varNode, dimNode

    allocate(self%Names%namingfilename, source = filename)
    do i=1, size(self%Names%namingfilename)
        call XMLReader%getFile(xmldoc,self%Names%namingfilename(i))
        outext='-->Setting naming conventions from '//self%Names%namingfilename(i)
        call Log%put(outext)
        tag="naming"          !base document node
        call XMLReader%gotoNode(xmldoc,xmldoc,tag)
        tag="variables"
        call XMLReader%gotoNode(xmldoc,varNode,tag)
        tag="dimensions"
        call XMLReader%gotoNode(xmldoc,dimNode,tag)
        call self%setVarNames(varNode)
        call self%setDimNames(dimNode)
    end do

    end subroutine setNamingConventions

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> set the variables naming conventions.
    !---------------------------------------------------------------------------
    subroutine setVarNames(self, varNode)
    class(globals_class), intent(inout) :: self
    type(Node), pointer, intent(in) :: varNode
    type(string) :: tag

    tag="bathymetry"
    call self%setCurrVar(tag, self%Var%bathymetry, self%Var%bathymetryVariants, varNode)
    tag="eastward_sea_water_velocity"
    call self%setCurrVar(tag, self%Var%u, self%Var%uVariants, varNode)
    tag="northward_sea_water_velocity"
    call self%setCurrVar(tag, self%Var%v, self%Var%vVariants, varNode)
    tag="upward_sea_water_velocity"
    call self%setCurrVar(tag, self%Var%w, self%Var%wVariants, varNode)
    tag="sea_surface_height"
    call self%setCurrVar(tag, self%Var%ssh, self%Var%sshVariants, varNode)
    tag="sea_water_temperature"
    call self%setCurrVar(tag, self%Var%temp, self%Var%tempVariants, varNode)
    tag="sea_water_salinity"
    call self%setCurrVar(tag, self%Var%sal, self%Var%salVariants, varNode)
    tag="sea_water_density"
    call self%setCurrVar(tag, self%Var%density, self%Var%densityVariants, varNode)
    tag="sea_surface_wave_stokes_drift_x_velocity"
    call self%setCurrVar(tag, self%Var%vsdx, self%Var%vsdxVariants, varNode)
    tag="sea_surface_wave_stokes_drift_y_velocity"
    call self%setCurrVar(tag, self%Var%vsdy, self%Var%vsdyVariants, varNode)
    tag="sea_surface_significant_wave_height"
    call self%setCurrVar(tag, self%Var%hs, self%Var%hsVariants, varNode)
    tag="sea_surface_wave_mean_period"
    call self%setCurrVar(tag, self%Var%ts, self%Var%tsVariants, varNode)
    tag="sea_surface_wave_mean_direction"
    call self%setCurrVar(tag, self%Var%wd, self%Var%wdVariants, varNode)
    tag="eastward_wind"
    call self%setCurrVar(tag, self%Var%u10, self%Var%u10Variants, varNode)
    tag="northward_wind"
    call self%setCurrVar(tag, self%Var%v10, self%Var%v10Variants, varNode)
    tag="surface_radiation"
    call self%setCurrVar(tag, self%Var%rad, self%Var%radVariants, varNode)
    tag="dissolved_oxygen"
    call self%setCurrVar(tag, self%Var%dissolved_oxygen, self%Var%dissolved_oxygenVariants, varNode)
    tag="nitrate"
    call self%setCurrVar(tag, self%Var%nitrate, self%Var%nitrateVariants, varNode)
    tag="nitrite"
    call self%setCurrVar(tag, self%Var%nitrite, self%Var%nitriteVariants, varNode)
    tag="ammonia"
    call self%setCurrVar(tag, self%Var%ammonia, self%Var%ammoniaVariants, varNode)
    tag="partOrgNit"
    call self%setCurrVar(tag, self%Var%partOrgNit, self%Var%partOrgNitVariants, varNode)
    tag="partOrgPho"
    call self%setCurrVar(tag, self%Var%partOrgPho, self%Var%partOrgPhoVariants, varNode)
    tag="DON_NonRefractory"
    call self%setCurrVar(tag, self%Var%DON_NonRefractory, self%Var%DON_NonRefractoryVariants, varNode)
    tag="DOP_NonRefractory"
    call self%setCurrVar(tag, self%Var%DOP_NonRefractory, self%Var%DOP_NonRefractoryVariants, varNode)
    tag="DON_Refractory"
    call self%setCurrVar(tag, self%Var%DON_Refractory, self%Var%DON_RefractoryVariants, varNode)
    tag="DOP_Refractory"
    call self%setCurrVar(tag, self%Var%DOP_Refractory, self%Var%DOP_RefractoryVariants, varNode)
    tag="phytoplankton"
    call self%setCurrVar(tag, self%Var%phytoplankton, self%Var%phytoplanktonVariants, varNode)
    tag="zooplankton"
    call self%setCurrVar(tag, self%Var%zooplankton, self%Var%zooplanktonVariants, varNode)
    tag="inorganic_phosphorus"
    call self%setCurrVar(tag, self%Var%inorganic_phosphorus, self%Var%inorganic_phosphorusVariants, varNode)
    tag="openpoints"
    call self%setCurrVar(tag, self%Var%openpoints, self%Var%openpointsVariants, varNode)
    
    end subroutine setVarNames

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> set the dimensions naming conventions.
    !---------------------------------------------------------------------------
    subroutine setDimNames(self, dimNode)
    class(globals_class), intent(inout) :: self
    type(Node), pointer, intent(in) :: dimNode
    type(string) :: tag

    tag="longitude"
    call self%setCurrVar(tag, self%Var%lon, self%Var%lonVariants, dimNode)
    tag="latitude"
    call self%setCurrVar(tag, self%Var%lat, self%Var%latVariants, dimNode)
    tag="vertical"
    call self%setCurrVar(tag, self%Var%level, self%Var%levelVariants, dimNode)
    tag="time"
    call self%setCurrVar(tag, self%Var%time, self%Var%timeVariants, dimNode)

    end subroutine setDimNames

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> set the current variable naming conventions. Imports the variable names from
    !> given .xml naming file.
    !> If the same variable is given again two things happen:
    !> a) the system name is overwriten;
    !> b) new variants are added to the name variant list
    !---------------------------------------------------------------------------
    subroutine setCurrVar(self, tag, currVar, currVarNameList, varNode)
    class(globals_class), intent(inout) :: self
    type(string), intent(in) :: tag
    type(string), intent(inout) :: currVar
    type(stringList_class), intent(inout) :: currVarNameList
    type(Node), pointer, intent(in) :: varNode
    integer :: i
    type(string) :: attValue, attName
    type(Node), pointer :: tempNode, variantNode
    type(NodeList), pointer :: varNameList

    call XMLReader%gotoNode(varNode, tempNode, tag, mandatory = .false.)
    if (associated(tempNode)) then !variable description exists in file
        attName="name"
        call XMLReader%getNodeAttribute(varNode, tag, attName, attValue, mandatory = .true.)
        currVar = attValue
        varNameList => getElementsByTagname(tempNode, "variant")
        do i = 0, getLength(varNameList) - 1
            variantNode => item(varNameList, i)
            call XMLReader%getLeafAttribute(variantNode,attName,attValue)
            if (currVarNameList%notRepeated(attValue)) call currVarNameList%add(attValue)
        end do
    end if

    end subroutine setCurrVar
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> 
    !---------------------------------------------------------------------------
    subroutine fillBackgroundDict(self, idx, variables)
    class(globals_class), intent(inout) :: self
    integer, intent(in) :: idx
    type(stringList_class), intent(inout) :: variables
    integer :: i
    !write(*,*)"Entrei no fillBackgroundDict"
    !write(*,*)"tamanho do dicionario = ", size(self%BackgroundVarDict%bkgDict)
    do i=1, size(self%BackgroundVarDict%bkgDict)
        if (self%BackgroundVarDict%bkgDict(i)%bkgIndex == 0) then
            self%BackgroundVarDict%bkgDict(i)%bkgIndex = idx
            !write(*,*)"Index e i = ", idx, i
            call self%BackgroundVarDict%bkgDict(i)%bkgVars%copy(variables)
            exit
        else if(self%BackgroundVarDict%bkgDict(i)%bkgIndex == idx) then
            !write(*,*)"Index e i = ", idx, i
            call self%BackgroundVarDict%bkgDict(i)%bkgVars%copy(variables)
            exit
        end if
    end do
    
    end subroutine fillBackgroundDict

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> initializes time stamp and date stamp variables
    !---------------------------------------------------------------------------
    subroutine setTimeDate(self)
    class(globals_class), intent(inout) :: self
    self%SimTime%TimeMax = self%Parameters%TimeMax
    self%SimTime%StartDate = self%Parameters%StartTime
    self%SimTime%EndDate = self%Parameters%EndTime
    self%SimTime%BaseDateTime = self%Parameters%BaseDateTime
    self%SimTime%CurrDate = self%SimTime%StartDate
    self%SimTime%CurrTime = 0
    end subroutine setTimeDate

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> sets the current time stamp and date stamp given a dt increment
    !---------------------------------------------------------------------------
    subroutine setCurrDateTime(self, dt)
    class(sim_time_t), intent(inout) :: self
    real(prec), intent(in) :: dt
    type(timedelta) :: step
    step = timedelta(seconds=floor(dt), milliseconds=floor((dt-floor(dt))*1000))
    self%CurrDate = self%CurrDate + step
    self%CurrTime = self%CurrTime + dt
    end subroutine setCurrDateTime

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns the time stamp for the current time/date. Time Stamp is days
    !> since epoch, given by BaseDateTime. Defaults to Copernicus standard.
    !---------------------------------------------------------------------------
    real(prec) function getDateTimeStamp(self)
    class(sim_time_t), intent(inout) :: self
    type(timedelta) :: step, day
    day = timedelta(days=1)
    step = self%CurrDate - self%BaseDateTime
    getDateTimeStamp = step%total_seconds()/day%total_seconds()
    end function getDateTimeStamp

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> prints date & time information.
    !---------------------------------------------------------------------------
    subroutine printDateTime(self)
    class(sim_time_t), intent(in) :: self
    type(string) :: outext
    type(string) :: temp_str
    character(len=23) :: temp_char
    temp_char = self%StartDate%isoformat(' ')
    temp_str = temp_char
    outext = outext//'      Start date                = '//temp_str//new_line('a')
    temp_char = self%EndDate%isoformat(' ')
    temp_str = temp_char
    outext = outext//'       End date                  = '//temp_str//new_line('a')
    temp_char = self%BaseDateTime%isoformat(' ')
    temp_str = temp_char
    outext = outext//'       Base date (time stamping) = '//temp_str//new_line('a')
    temp_str=self%TimeMax
    outext = outext//'       Simulation will run for '//temp_str//' s'
    call Log%put(outext,.false.)
    end subroutine printDateTime

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Increments Tracer count. This routine MUST be ATOMIC.
    !---------------------------------------------------------------------------
    subroutine increment_numTracer(self)
    class(sim_t), intent(inout) :: self
    !ATOMIC pragma here please
    self%numTracer = self%numTracer + 1
    end subroutine increment_numTracer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns a new ID for a Tracer.
    !---------------------------------------------------------------------------
    integer function getnumTracer(self)
    class(sim_t), intent(inout) :: self
    call self%increment_numTracer()
    getnumTracer = self%numTracer
    end function getnumTracer

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> incrementing time step count.
    !---------------------------------------------------------------------------
    subroutine increment_numdt(self)
    class(sim_t), intent(inout) :: self
    self%numdt = self%numdt + 1
    end subroutine increment_numdt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the number of time steps.
    !---------------------------------------------------------------------------
    integer function getnumdt(self)
    class(sim_t), intent(inout) :: self
    getnumdt = self%numdt
    end function getnumdt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the number of time steps from last ouptut.
    !---------------------------------------------------------------------------
    integer function getlastOutNumDt(self)
    class(output_t), intent(inout) :: self
    getlastOutNumDt = self%lastOutNumDt
    end function getlastOutNumDt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Sets the number of time steps from last ouptut.
    !---------------------------------------------------------------------------
    subroutine setlastOutNumDt(self, num)
    class(output_t), intent(inout) :: self
    integer, intent(in) :: num
    self%lastOutNumDt = num
    end subroutine setlastOutNumDt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> incrementing output file count.
    !---------------------------------------------------------------------------
    subroutine increment_numOutFile(self)
    class(output_t), intent(inout) :: self
    self%numOutFile = self%numOutFile + 1
    end subroutine increment_numOutFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the number of output files written.
    !---------------------------------------------------------------------------
    integer function getnumOutFile(self)
    class(output_t), intent(inout) :: self
    getnumOutFile = self%numOutFile
    end function getnumOutFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> adds a variable name to the output pool.
    !---------------------------------------------------------------------------
    subroutine addToOutputPool(self, var)
    class(output_t), intent(inout) :: self
    type(string), intent(in) :: var
    if (self%varToOutput%notRepeated(var)) call self%varToOutput%add(var)
    end subroutine addToOutputPool
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> adds variable fields one by one to the output pool, given a string and 
    !> boolean list.
    !---------------------------------------------------------------------------
    subroutine setOutputFields(self, fields, toOutput)
    class(output_t), intent(inout) :: self
    type(string), dimension(:), intent(in) :: fields
    logical, dimension(:), intent(in) :: toOutput
    integer :: i
    if (size(fields) > 0) then
        do i= 1, size(fields)
            if (toOutput(i)) call self%addToOutputPool(fields(i))
        end do
    end if
    end subroutine setOutputFields

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> adds a variable name to the output pool.
    !---------------------------------------------------------------------------
    subroutine getOutputPoolArray(self, array)
    class(output_t), intent(inout) :: self
    type(string), intent(out), allocatable, dimension(:) :: array
    class(*), pointer :: curr
    integer :: i
    type(string) :: outext
    i=1
    allocate(array(self%varToOutput%getSize()))
    call self%varToOutput%reset()               ! reset list iterator
    do while(self%varToOutput%moreValues())     ! loop while there are values to print
        curr => self%varToOutput%currentValue() ! get current value
        select type(curr)
        class is (string)
            array(i) = curr
            class default
            outext = '[Globals::getOutputPoolArray] Unexepected type of content, not a string'
            call Log%put(outext)
            stop
        end select
        call self%varToOutput%next()            ! increment the list iterator
        i = i + 1
    end do
    call self%varToOutput%reset()               ! reset list iterator
    end subroutine getOutputPoolArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns true if a given var is in the output pool.
    !---------------------------------------------------------------------------
    logical function inOutputPool(self, var)
    class(output_t), intent(inout) :: self
    type(string), intent(in) :: var
    inOutputPool = .not.self%varToOutput%notRepeated(var)
    end function inOutputPool

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private parameter setting method. Builds the simulation parametric space
    !> from the input case file.
    !> @param[in] self, parmkey, parmvalue
    !---------------------------------------------------------------------------
    subroutine setParam(self,parmkey,parmvalue)
    class(parameters_t), intent(inout) :: self
    type(string), intent(in) :: parmkey
    type(string), intent(in) :: parmvalue
    type(string), allocatable :: dc(:)
    integer :: i
    integer, dimension(6) :: date
    integer :: sizem
    !add new parameters to this search
    if (parmkey%chars()=="Integrator") then
        self%Integrator=parmvalue%to_number(kind=1_I1P)
        sizem=sizeof(self%Integrator)
    elseif(parmkey%chars()=="Threads") then
        if (parmvalue /= 'auto') then
            self%numOPMthreads=parmvalue%to_number(kind=1_I1P)
            call OMPManager%setThreads(self%numOPMthreads)
        end if
        sizem=sizeof(self%WarmUpTime)
    elseif(parmkey%chars()=="WarmUpTime") then
        self%WarmUpTime=parmvalue%to_number(kind=1._R8P)
        sizem=sizeof(self%WarmUpTime)
    elseif(parmkey%chars()=="OutputWriteTime") then
        self%OutputWriteTime=parmvalue%to_number(kind=1._R8P)
        sizem=sizeof(self%OutputWriteTime)
    elseif(parmkey%chars()=="Start") then
        date = Utils%getDateFromISOString(parmvalue)
        self%StartTime = Utils%getDateTimeFromDate(date)
        if (.not. self%StartTime%isValid()) self%StartTime = datetime()
        sizem=sizeof(self%StartTime)
    elseif(parmkey%chars()=="End") then
        date = Utils%getDateFromISOString(parmvalue)
        self%EndTime = Utils%getDateTimeFromDate(date)
        if (.not. self%EndTime%isValid()) self%EndTime = datetime()
        sizem=sizeof(self%EndTime)
    elseif(parmkey%chars()=="BaseDateTime") then
        date = Utils%getDateFromISOString(parmvalue)
        self%BaseDateTime = Utils%getDateTimeFromDate(date)
        if (.not. self%BaseDateTime%isValid()) self%BaseDateTime = datetime()
        sizem=sizeof(self%BaseDateTime)
    elseif(parmkey%chars()=="BufferSize") then
        self%BufferSize=parmvalue%to_number(kind=1_I_P)
        sizem=sizeof(self%BufferSize)
    elseif(parmkey%chars()=="OutputFormat") then
        self%OutputFormat=parmvalue%to_number(kind=1_I1P)
        sizem=sizeof(self%OutputFormat)
    end if
    call SimMemory%adddef(sizem)

    end subroutine setParam

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Parameter checking method. Checks if mandatory parameters were set
    !---------------------------------------------------------------------------
    subroutine check(self)
    class(parameters_t), intent(inout) :: self
    type(string) :: outext
    type(datetime) :: temp
    type(timedelta) :: simtime

    if ( any(self%IntegratorIndexes == self%Integrator)) then
    else
        outext = '[Globals::parameters::check]: Integrator not recognized, stoping'
        call Log%put(outext)
        stop
    end if
    if ( any(self%OutputFormatIndexes == self%OutputFormat)) then
        if (self%OutputFormat == 1) then
            outext = '[Globals::parameters::check]: NetCDF is not implemented yet, try something nicer like VTK, stoping'
            call Log%put(outext)
            stop
        end if
    else
        outext = '[Globals::parameters::check]: OutputFormat not recognized, stoping'
        call Log%put(outext)
        stop
    end if
    temp = datetime() !default initialization
    !add new parameters to this search
    if (self%OutputWriteTime==MV) then
        outext = '[Globals::parameters::check]: Output sampling rate parameter (OutputWriteTime) is not set, stoping'
        call Log%put(outext)
        stop
    elseif (self%StartTime==temp) then
        outext = '[Globals::parameters::check]: start time parameter (Start) is not set or invalid, stoping'
        call Log%put(outext)
        stop
    elseif (self%EndTime==temp) then
        outext = '[Globals::parameters::check]: end time parameter (End) is not set or invalid, stoping'
        call Log%put(outext)
        stop
    end if
    !Build timemax from the difference between start and end time
    simtime = self%EndTime - self%StartTime
    self%TimeMax = simtime%total_seconds()
    end subroutine check

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Parameter printing method.
    !---------------------------------------------------------------------------
    subroutine printsimparameters(self)
    implicit none
    class(parameters_t), intent(inout) :: self
    type(string) :: outext
    type(string) :: temp_str
    character(len=23) :: temp_char
    outext = '      Integrator scheme is '//self%IntegratorNames(self%Integrator)//new_line('a')
    temp_str=self%numOPMthreads
    temp_str=OMPManager%getThreads()
    outext = outext//'       OMP threads = '//temp_str//new_line('a')
    temp_str=self%WarmUpTime
    outext = outext//'       WarmUpTime = '//temp_str//' s'//new_line('a')
    temp_str=self%OutputWriteTime
    outext = outext//'       OutputWriteTime = '//temp_str//' Hz'//new_line('a')
    outext = outext//'       Output file format is '//self%OutputFormatNames(self%OutputFormat)//new_line('a')
    temp_str=self%BufferSize
    outext = outext//'       Input buffer size is '//temp_str//' s'
    call Log%put(outext,.false.)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Gravity setting routine.
    !> @param[in] self, grav
    !---------------------------------------------------------------------------
    subroutine setgravity (self,grav)
    implicit none
    class(constants_t), intent(inout) :: self
    type(vector), intent(in) :: grav
    integer :: sizem
    type(string) :: outext
    self%Gravity= grav
    if (grav%x==MV) then !Gravity was not read, setting default
        self%Gravity= -9.81*ez
        outext = '       Gravity not specified, setting to default value = (0,0,-9.81)'
        call Log%put(outext,.false.)
    endif
    sizem=sizeof(self%Gravity)
    call SimMemory%adddef(sizem)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Z0 setting routine.
    !> @param[in] self, read_z0
    !---------------------------------------------------------------------------
    subroutine setz0(self,read_z0)
    implicit none
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_z0
    integer :: sizem
    self%Z0=read_z0%to_number(kind=1._R8P)
    sizem = sizeof(self%Z0)
    call SimMemory%adddef(sizem)
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Beaching Level setting routine.
    !> @param[in] self, read_BeachingLevel
    !---------------------------------------------------------------------------
    subroutine setBeachingLevel(self, read_BeachingLevel)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_BeachingLevel
    type(string) :: outext
    integer :: sizem
    if (read_BeachingLevel%to_number(kind=1._R8P) > 0.0) then
        outext='Beaching level must be negative, assuming default value'
        call Log%put(outext)
    else
        self%BeachingLevel=read_BeachingLevel%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%BeachingLevel)
    call SimMemory%adddef(sizem)
    end subroutine setBeachingLevel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Beaching stop probability setting routine.
    !> @param[in] self, read_BeachingStopProb
    !---------------------------------------------------------------------------
    subroutine setBeachingStopProb(self, read_BeachingStopProb)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_BeachingStopProb
    type(string) :: outext
    integer :: sizem
    if (read_BeachingStopProb%to_number(kind=1._R8P) < 0.0) then
        outext='Beaching stopping probability must be zero or positive, assuming default value'
        call Log%put(outext)
    else
        self%BeachingStopProb =read_BeachingStopProb%to_number(kind=1._R8P)*0.01 !user input is in %
    endif
    sizem = sizeof(self%BeachingStopProb)
    call SimMemory%adddef(sizem)
    end subroutine setBeachingStopProb

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Horizontal Diffusion coefficient setting routine.
    !> @param[in] self, read_DiffusionCoeff
    !---------------------------------------------------------------------------
    subroutine setDiffusionCoeff(self, read_DiffusionCoeff)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_DiffusionCoeff
    type(string) :: outext
    integer :: sizem
    if (read_DiffusionCoeff%to_number(kind=1._R8P) < 0.0) then
        outext='Diffusion coefficient must be zero or positive, assuming default value'
        call Log%put(outext)
    else
        self%DiffusionCoeff =read_DiffusionCoeff%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%DiffusionCoeff)
    call SimMemory%adddef(sizem)
    end subroutine setDiffusionCoeff
    !------------------------------------------------------------------------------
    
    !> @author Joao Sobrinho
    !> @brief
    !> Wind drag coefficient setting routine.
    !> @param[in] self, read_WindDragCoeff
    !---------------------------------------------------------------------------
    subroutine setWindDragCoeff(self, read_WindDragCoeff)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_WindDragCoeff
    type(string) :: outext
    integer :: sizem
    if (read_WindDragCoeff%to_number(kind=1._R8P) < 0.0) then
        outext='Diffusion coefficient must be zero or positive, assuming default value'
        call Log%put(outext)
    else
        self%WindDragCoeff =read_WindDragCoeff%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%WindDragCoeff)
    call SimMemory%adddef(sizem)
    end subroutine setWindDragCoeff

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Resuspension setting routine or not
    !> @param[in] self, read_ResuspensionCoeff
    !---------------------------------------------------------------------------
    subroutine setResuspensionCoeff(self, read_ResuspensionCoeff)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_ResuspensionCoeff
    type(string) :: outext
    integer :: sizem
    if (read_ResuspensionCoeff%to_number(kind=1._R8P) < 0.0) then
        outext='Resuspension factor must be zero or positive, assuming default value'
        call Log%put(outext)
    else
        self%ResuspensionCoeff=read_ResuspensionCoeff%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%ResuspensionCoeff)
    call SimMemory%adddef(sizem)
    end subroutine setResuspensionCoeff
    
    !> @author Joao Sobrinho - ColabAtlantic
    !> @brief
    !> Rugosity setting routine or not
    !> @param[in] self, read_Rugosity
    !---------------------------------------------------------------------------
    subroutine setRugosity(self, read_Rugosity)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_Rugosity
    type(string) :: outext
    integer :: sizem
    if (read_Rugosity%to_number(kind=1._R8P) < 0.0) then
        outext='Rugosity must be zero or positive, assuming default value'
        call Log%put(outext)
    else
        self%Rugosity=read_Rugosity%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%Rugosity)
    call SimMemory%adddef(sizem)
    end subroutine setRugosity
    
    !> @author Joao Sobrinho - ColabAtlantic
    !> @brief
    !> Critical shear erosion setting routine or not
    !> @param[in] self, read_Critical_Shear_Erosion
    !---------------------------------------------------------------------------
    subroutine setCritical_Shear_Erosion(self, read_Critical_Shear_Erosion)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_Critical_Shear_Erosion
    type(string) :: outext
    integer :: sizem
    if (read_Critical_Shear_Erosion%to_number(kind=1._R8P) < 0.0) then
        outext='Critical_Shear_Erosion must be zero or positive, assuming default value'
        call Log%put(outext)
    else
        self%Critical_Shear_Erosion=read_Critical_Shear_Erosion%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%Critical_Shear_Erosion)
    call SimMemory%adddef(sizem)
    end subroutine setCritical_Shear_Erosion
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> smallDt setting routine.
    !> @param[in] self, dt
    !---------------------------------------------------------------------------
    subroutine setSmallDt(self, dt)
    class(constants_t), intent(inout) :: self
    real(prec), intent(in) :: dt
    type(string) :: outext
    integer :: sizem
    self%smallDt=dt/3.0
    if (self%smallDt <= 0.0) then
        outext='dt must be positive and non-zero, stopping'
        call Log%put(outext)
        stop
    endif
    sizem = sizeof(self%smallDt)
    call SimMemory%adddef(sizem)
    end subroutine setSmallDt
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Beaching stop probability setting routine.
    !> @param[in] self, read_MeanDensity
    !---------------------------------------------------------------------------
    subroutine setMeanDensity(self, read_MeanDensity)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_MeanDensity
    type(string) :: outext
    integer :: sizem
    if (read_MeanDensity%to_number(kind=1._R8P) <= 0.0) then
        outext='Mean medium density must be positive, assuming default value'
        call Log%put(outext)
    else
        self%MeanDensity =read_MeanDensity%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%MeanDensity)
    call SimMemory%adddef(sizem)
    end subroutine setMeanDensity
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Beaching stop probability setting routine.
    !> @param[in] self, read_MeanDensity
    !---------------------------------------------------------------------------
    subroutine setMeanKVisco(self, read_MeanKVisco)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_MeanKVisco
    type(string) :: outext
    integer :: sizem
    if (read_MeanKVisco%to_number(kind=1._R8P) <= 0.0) then
        outext='Mean medium kinematic viscosity must be positive, assuming default value'
        call Log%put(outext)
    else
        self%MeanKVisco =read_MeanKVisco%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%MeanKVisco)
    call SimMemory%adddef(sizem)
    end subroutine setMeanKVisco
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Add bottom cell to enable interpolation until the bottom.
    !> @param[in] self, read_AddBottomCell
    !---------------------------------------------------------------------------
    subroutine setAddBottomCell(self, read_AddBottomCell)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_AddBottomCell
    integer :: sizem
    self%AddBottomCell = read_AddBottomCell%to_number(kind=1_I1P)
    sizem = sizeof(self%AddBottomCell)
    call SimMemory%adddef(sizem)
    end subroutine setAddBottomCell
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for TOptBacteriaMin.
    !> @param[in] self, read_TOptBacteriaMin
    !---------------------------------------------------------------------------
    subroutine setTOptBacteriaMin(self, read_TOptBacteriaMin)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_TOptBacteriaMin
    type(string) :: outext
    integer :: sizem
    if (read_TOptBacteriaMin%to_number(kind=1._R8P) <= 0.0) then
        outext='TOptBacteriaMin must be positive, assuming default value'
        call Log%put(outext)
    else
        self%TOptBacteriaMin =read_TOptBacteriaMin%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%TOptBacteriaMin)
    call SimMemory%adddef(sizem)
    end subroutine setTOptBacteriaMin
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for TOptBacteriaMax.
    !> @param[in] self, read_TOptBacteriaMax
    !---------------------------------------------------------------------------
    subroutine setTOptBacteriaMax(self, read_TOptBacteriaMax)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_TOptBacteriaMax
    type(string) :: outext
    integer :: sizem
    if (read_TOptBacteriaMax%to_number(kind=1._R8P) <= 0.0) then
        outext='TOptBacteriaMax must be positive, assuming default value'
        call Log%put(outext)
    else
        self%TOptBacteriaMax =read_TOptBacteriaMax%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%TOptBacteriaMax)
    call SimMemory%adddef(sizem)
    end subroutine setTOptBacteriaMax
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for TBacteriaMin.
    !> @param[in] self, read_TBacteriaMin
    !---------------------------------------------------------------------------
    subroutine setTBacteriaMin(self, read_TBacteriaMin)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_TBacteriaMin
    type(string) :: outext
    integer :: sizem
    if (read_TBacteriaMin%to_number(kind=1._R8P) <= 0.0) then
        outext='TBacteriaMin must be positive, assuming default value'
        call Log%put(outext)
    else
        self%TBacteriaMin =read_TBacteriaMin%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%TBacteriaMin)
    call SimMemory%adddef(sizem)
    end subroutine setTBacteriaMin
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for TBacteriaMax.
    !> @param[in] self, read_TBacteriaMax
    !---------------------------------------------------------------------------
    subroutine setTBacteriaMax(self, read_TBacteriaMax)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_TBacteriaMax
    type(string) :: outext
    integer :: sizem
    if (read_TBacteriaMax%to_number(kind=1._R8P) <= 0.0) then
        outext='TBacteriaMax must be positive, assuming default value'
        call Log%put(outext)
    else
        self%TBacteriaMax =read_TBacteriaMax%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%TBacteriaMax)
    call SimMemory%adddef(sizem)
    end subroutine setTBacteriaMax
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for BK1.
    !> @param[in] self, read_BK1
    !---------------------------------------------------------------------------
    subroutine setBK1(self, read_BK1)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_BK1
    type(string) :: outext
    integer :: sizem
    if (read_BK1%to_number(kind=1._R8P) <= 0.0) then
        outext='BK1 must be positive, assuming default value'
        call Log%put(outext)
    else
        self%BK1 =read_BK1%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%BK1)
    call SimMemory%adddef(sizem)
    end subroutine setBK1
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for BK2.
    !> @param[in] self, read_BK2
    !---------------------------------------------------------------------------
    subroutine setBK2(self, read_BK2)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_BK2
    type(string) :: outext
    integer :: sizem
    if (read_BK2%to_number(kind=1._R8P) <= 0.0) then
        outext='BK2 must be positive, assuming default value'
        call Log%put(outext)
    else
        self%BK2 =read_BK2%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%BK2)
    call SimMemory%adddef(sizem)
    end subroutine setBK2
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for BK3.
    !> @param[in] self, read_BK3
    !---------------------------------------------------------------------------
    subroutine setBK3(self, read_BK3)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_BK3
    type(string) :: outext
    integer :: sizem
    if (read_BK3%to_number(kind=1._R8P) <= 0.0) then
        outext='BK3 must be positive, assuming default value'
        call Log%put(outext)
    else
        self%BK3 =read_BK3%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%BK3)
    call SimMemory%adddef(sizem)
    end subroutine setBK3
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for BK4.
    !> @param[in] self, read_BK4
    !---------------------------------------------------------------------------
    subroutine setBK4(self, read_BK4)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_BK4
    type(string) :: outext
    integer :: sizem
    if (read_BK4%to_number(kind=1._R8P) <= 0.0) then
        outext='BK4 must be positive, assuming default value'
        call Log%put(outext)
    else
        self%BK4 =read_BK4%to_number(kind=1._R8P)
    endif
    sizem = sizeof(self%BK4)
    call SimMemory%adddef(sizem)
    end subroutine setBK4
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> defines value for bacterial maximum growth rate (1/d).
    !> @param[in] self, read_MaxDegradationRate
    !---------------------------------------------------------------------------
    subroutine setMaxDegradationRate(self, read_MaxDegradationRate)
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_MaxDegradationRate
    type(string) :: outext
    integer :: sizem
    if (read_MaxDegradationRate%to_number(kind=1._R8P) <= 0.0) then
        outext='MaxDegradationRate must be positive, assuming default value'
        call Log%put(outext)
    else
        self%MaxDegradationRate =read_MaxDegradationRate%to_number(kind=1._R8P) / 86400.0
    endif
    sizem = sizeof(self%MaxDegradationRate)
    call SimMemory%adddef(sizem)
    end subroutine setMaxDegradationRate

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public constants printing routine.
    !---------------------------------------------------------------------------
    subroutine printconstants(self)
    implicit none
    class(constants_t), intent(in) :: self
    type(string) :: outext
    type(string) :: temp_str(3)

    temp_str(1)=self%Gravity%x
    temp_str(2)=self%Gravity%y
    temp_str(3)=self%Gravity%z
    outext = '      Gravity is '//new_line('a')//&
        '       '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')
    temp_str(1)=self%Z0
    outext = outext//'       Z0 = '//temp_str(1)//' m'//new_line('a')
    temp_str(1)=self%BeachingLevel
    outext = outext//'       BeachingLevel = '//temp_str(1)//' m'//new_line('a')
    temp_str(1)=self%BeachingStopProb
    outext = outext//'       BeachingStopProb = '//temp_str(1)//' -'//new_line('a')
    temp_str(1)=self%DiffusionCoeff
    outext = outext//'       DiffusionCoeff = '//temp_str(1)//' -'//new_line('a')  
    temp_str(1)=self%WindDragCoeff
    outext = outext//'       WindDragCoeff = '//temp_str(1)//' -'//new_line('a')  
    temp_str(1)=self%ResuspensionCoeff
    outext = outext//'       ResuspensionCoeff = '//temp_str(1)//new_line('a')
    temp_str(1)=self%Critical_Shear_Erosion
    outext = outext//'       CriticalShearErosion = '//temp_str(1)//new_line('a')
    temp_str(1)=self%AddBottomCell
    outext = outext//'       AddBottomCell = '//temp_str(1)//new_line('a')
    temp_str(1)=self%MeanDensity
    outext = outext//'       MeanDensity = '//temp_str(1)//new_line('a')
    temp_str(1)=self%MeanKVisco
    outext = outext//'       MeanKVisco = '//temp_str(1)//new_line('a')
    temp_str(1)=self%TOptBacteriaMin
    outext = outext//'       TOptBacteriaMin = '//temp_str(1)//new_line('a')
    temp_str(1)=self%TOptBacteriaMax
    outext = outext//'       TOptBacteriaMax = '//temp_str(1)//new_line('a')
    temp_str(1)=self%TBacteriaMin
    outext = outext//'       TBacteriaMin = '//temp_str(1)//new_line('a')
    temp_str(1)=self%TBacteriaMax
    outext = outext//'       TBacteriaMax = '//temp_str(1)//new_line('a')
    temp_str(1)=self%BK1
    outext = outext//'       BK1 = '//temp_str(1)//new_line('a')
    temp_str(1)=self%BK2
    outext = outext//'       BK2 = '//temp_str(1)//new_line('a')
    temp_str(1)=self%BK3
    outext = outext//'       BK3 = '//temp_str(1)//new_line('a')
    temp_str(1)=self%BK4
    outext = outext//'       BK4 = '//temp_str(1)//new_line('a')
    temp_str(1)=self%MaxDegradationRate
    outext = outext//'       MaxDegradationRate = '//temp_str(1)//''
    call Log%put(outext,.false.)
    end subroutine printconstants

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Dp setting routine.
    !> @param[in] self, read_dp
    !---------------------------------------------------------------------------
    subroutine setdp(self,read_dp)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(vector), intent(in) :: read_dp
    integer :: sizem
    self%Dp=read_dp
    sizem = sizeof(self%Dp)
    call SimMemory%adddef(sizem)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Dt setting routine
    !> @param[in] self, read_dt
    !---------------------------------------------------------------------------
    subroutine setdt(self,read_dt)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_dt
    type(string) :: outext
    integer :: sizem
    self%dt=read_dt%to_number(kind=1._R8P)
    if (self%dt <= 0.0) then
        outext='dt must be positive and non-zero, stopping'
        call Log%put(outext)
        stop
    endif
    sizem = sizeof(self%dt)
    call SimMemory%adddef(sizem)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Water quality Dt setting routine
    !> @param[in] self, read_WqDt
    !---------------------------------------------------------------------------
    subroutine setWqDt(self, read_WqDt)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(string) :: outext
    integer :: sizem
    real(prec) :: read_WqDt
    self%WqDt=read_WqDt
    if (self%WqDt <= 0.0) then
        outext='water quality Dt (defined in waterquality.dat - DT_SECONDS) must be positive and non-zero, stopping'
        call Log%put(outext)
        stop
    endif
    sizem = sizeof(self%WqDt)
    call SimMemory%adddef(sizem)
    end subroutine setWqDt
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Bounding box setting routine
    !> @param[in] self, point_, coords
    !---------------------------------------------------------------------------
    subroutine setboundingbox(self,point_, coords)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: point_
    type(vector) :: coords
    integer :: sizem
    if (point_%chars() == "BoundingBoxMin") then
        self%Pointmin= coords
    elseif (point_%chars() == "BoundingBoxMax") then
        self%Pointmax= coords
    endif
    sizem=sizeof(coords)
    call SimMemory%adddef(sizem)
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Set the center of the simulation
    !> @param[in] self, center
    !---------------------------------------------------------------------------
    subroutine setCenter(self, center)
    class(simdefs_t), intent(inout) :: self
    type(vector) :: center
    integer :: sizem
    self%Center= center
    sizem=sizeof(center)
    call SimMemory%adddef(sizem)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> blocksize box setting routine
    !> @param[in] self, bsize
    !---------------------------------------------------------------------------
    subroutine setblocksize(self, bsize)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(vector) :: bsize
    integer :: sizem
    self%blocksize = bsize
    sizem = sizeof(bsize)
    call SimMemory%adddef(sizem)
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Choose if tracers must be removed when they reach land or not
    !> @param[in] self, read_RemoveLandTracer
    !---------------------------------------------------------------------------
    subroutine setRemoveLandTracer(self, read_RemoveLandTracer)
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_RemoveLandTracer
    type(string) :: outext
    integer :: sizem
    if ((read_RemoveLandTracer%to_number(kind=1._I8P) < 0) .OR. (read_RemoveLandTracer%to_number(kind=1._I8P) > 1)) then
        outext='Remove land tracers option must be 0:No or 1:yes, assuming default value'
        call Log%put(outext)
    else
        self%RemoveLandTracer=read_RemoveLandTracer%to_number(kind=1._I8P)
    end if
    sizem = sizeof(self%RemoveLandTracer)
    call SimMemory%adddef(sizem)
    end subroutine setRemoveLandTracer
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Set bathymetry construct from the bathymetric netcdf property
    !> @param[in] self, read_bathyminNetcdf
    !---------------------------------------------------------------------------
    subroutine setbathyminNetcdf(self, read_bathyminNetcdf)
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_bathyminNetcdf
    type(string) :: outext
    integer :: sizem
    if ((read_bathyminNetcdf%to_number(kind=1._I8P) < 0) .OR. (read_bathyminNetcdf%to_number(kind=1._I8P) > 1)) then
        outext='read bathymetry from netcdf file must be 0:no or 1:yes, assuming default value'
        call Log%put(outext)
    else
        self%bathyminNetcdf=read_bathyminNetcdf%to_number(kind=1._I8P)
    end if
    sizem = sizeof(self%bathyminNetcdf)
    call SimMemory%adddef(sizem)
    end subroutine setbathyminNetcdf
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Set tracer max age
    !> @param[in] self, read_tracerMaxAge
    !---------------------------------------------------------------------------
    subroutine settracerMaxAge(self, read_tracerMaxAge)
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_tracerMaxAge
    type(string) :: outext
    integer :: sizem
    if ((read_tracerMaxAge%to_number(kind=1._R8P) < 0)) then
        outext='read tracerMaxAge must be greater or equal to 0: assuming default value 0'
        call Log%put(outext)
    else
        self%tracerMaxAge=read_tracerMaxAge%to_number(kind=1._R8P)
    end if
    sizem = sizeof(self%tracerMaxAge)
    call SimMemory%adddef(sizem)
    end subroutine settracerMaxAge

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Set Temperature_add_offset
    !> @param[in] self, read_Temperature_add_offset
    !---------------------------------------------------------------------------
    subroutine setTemperature_add_offset(self, read_Temperature_add_offset)
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_Temperature_add_offset
    integer :: sizem
    self%Temperature_add_offset=read_Temperature_add_offset%to_number(kind=1._R8P)
    
    sizem = sizeof(self%Temperature_add_offset)
    call SimMemory%adddef(sizem)
    end subroutine setTemperature_add_offset
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Resuspension setting routine or not
    !> @param[in] self, read_VerticalVelMethod
    !---------------------------------------------------------------------------
    subroutine setVerticalVelMethod(self, read_VerticalVelMethod)
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_VerticalVelMethod
    type(string) :: outext
    integer :: sizem
    if ((read_VerticalVelMethod%to_number(kind=1._I8P) < 0) .OR. (read_VerticalVelMethod%to_number(kind=1._I8P) > 3)) then
        outext='Vertical velocity method must be 1:From velocity fields, 2:Divergence or 3:Disabled, assuming default value'
        call Log%put(outext)
    else
        self%VerticalVelMethod=read_VerticalVelMethod%to_number(kind=1._I8P)
    endif
    sizem = sizeof(self%VerticalVelMethod)
    call SimMemory%adddef(sizem)
    end subroutine setVerticalVelMethod
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Is beaching method based on project FreeLitterAt enabled
    !> @param[in] self, read_FreeLitterAtBeaching
    !---------------------------------------------------------------------------
    subroutine setFreeLitterAtBeaching(self, read_FreeLitterAtBeaching)
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_FreeLitterAtBeaching
    type(string) :: outext
    integer :: sizem
    if ((read_FreeLitterAtBeaching%to_number(kind=1._I8P) < 0) .OR. (read_FreeLitterAtBeaching%to_number(kind=1._I8P) > 1)) then
        outext='FreeLitterAtBeaching must be either 0 or 1'
        call Log%put(outext)
    else
        self%FreeLitterAtBeaching=read_FreeLitterAtBeaching%to_number(kind=1._I8P)
    endif
    sizem = sizeof(self%FreeLitterAtBeaching)
    call SimMemory%adddef(sizem)
    end subroutine setFreeLitterAtBeaching
    
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> source inititialization proceadure - initializes Beaching variables
    !> @param[in] beach, id, coastType, probability, waterColumnThreshold, beachTimeScale, unbeach, unbeachTimeScale, runUpEffect,&
    !> beachSlope, runUpEffectUnbeach, beaching_geometry, shapetype
    !---------------------------------------------------------------------------
    subroutine initializeBeachAreas(beach, id, coastType, probability, waterColumnThreshold, &
                                beachTimeScale, unbeach, unbeachTimeScale, runUpEffect, &
                                beachSlope, runUpEffectUnbeach, beaching_geometry, shapetype)
    class(beach_class) :: beach
    integer, intent(in) :: id
    integer, intent(in) :: coastType
    integer, intent(in) :: unbeach
    integer, intent(in) :: runUpEffect
    integer, intent(in) :: runUpEffectUnbeach
    real(prec), intent(in) :: probability
    real(prec), intent(in) :: waterColumnThreshold
    real(prec), intent(in) :: beachTimeScale
    real(prec), intent(in) :: unbeachTimeScale
    real(prec), intent(in) :: beachSlope
    type(string), intent(in) :: beaching_geometry
    class(shape), intent(in) :: shapetype
    
    !Locals-----------------------------------------
    
    integer :: sizem, i
    type(string) :: outext
    integer :: err
    !begin------------------------------------------
    
    !Setting parameters
    beach%par%id=id
    beach%par%coastType=coastType
    beach%par%probability = probability
    beach%par%waterColumnThreshold = waterColumnThreshold
    beach%par%beachTimeScale = beachTimeScale
    beach%par%unbeach = unbeach
    beach%par%unbeachTimeScale = unbeachTimeScale
    beach%par%runUpEffect = runUpEffect
    !Setting possible variable position
    beach%par%beachSlope = beachSlope
    beach%par%runUpEffectUnbeach = runUpEffectUnbeach
    beach%par%beaching_geometry = beaching_geometry
    allocate(beach%par%geometry, source=shapetype)

    sizem = sizeof(beach)
    call SimMemory%addbeachArea(sizem)
    !TODO - Add print like it is done in sources.

    end subroutine initializeBeachAreas
    
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> source allocation routine - allocates beaching areas objects
    !> @param[in] self,nareas
    !---------------------------------------------------------------------------
    subroutine initAreas(self,nareas)
    implicit none
    class(beachArray_class), intent(inout) :: self
    integer, intent(in) :: nareas
    integer err
    type(string) :: outext, temp
    allocate(self%beachArea(nareas), stat=err)
    if(err/=0)then
        outext='[initAreas]: Cannot allocate Beach Areas, stoping'
        call Log%put(outext)
        stop
    else
        temp = nareas
        outext = 'Allocated '// temp // ' Areas.'
        call Log%put(outext)
    endif
    end subroutine initAreas

    !---------------------------------------------------------------------------
    !> @author Sobrinho
    !> @brief
    !> source group destructor - deallocates beach areas objects
    !---------------------------------------------------------------------------
    subroutine killAreas(self)
    implicit none
    class(beachArray_class), intent(inout) :: self
    integer err
    type(string) :: outext
    if (allocated(self%beachArea)) deallocate(self%beachArea, stat=err)
    if(err/=0)then
        outext='[killSources]: Cannot deallocate Sources, stoping'
        call Log%put(outext)
        stop
    endif
    end subroutine killAreas

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation definitions printing routine.
    !---------------------------------------------------------------------------
    subroutine printsimdefs(self)
    class(simdefs_t), intent(in) :: self
    type(string) :: outext
    type(string) :: temp_str(3)

    temp_str(1)=self%Dp%normL2()
    outext = '      Initial resolution is '//temp_str(1)//' m'//new_line('a')
    temp_str(1)=self%dt
    outext = '      Timestep is '//temp_str(1)//' s'//new_line('a')
    temp_str(1)=self%Pointmin%x
    temp_str(2)=self%Pointmin%y
    temp_str(3)=self%Pointmin%z
    outext = outext//'       Pointmin (BB) is '//new_line('a')//&
        '       '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')
    temp_str(1)=self%Pointmax%x
    temp_str(2)=self%Pointmax%y
    temp_str(3)=self%Pointmax%z
    outext = outext//'       Pointmax (BB) is '//new_line('a')//&
        '       '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')
    if (self%autoblocksize) then
        outext = outext//'       Blocks are automatically sized'//new_line('a')
    else
        temp_str(1)=self%blocksize%x
        temp_str(2)=self%blocksize%y
        outext = outext//'       Blocks are sized '//new_line('a')//&
            '       '//temp_str(1)//' X '//temp_str(2)//new_line('a')
    end if
    temp_str(1)=self%VerticalVelMethod
    outext = outext//'       VerticalVelMethod = '//temp_str(1)//new_line('a')
    temp_str(1)=self%RemoveLandTracer
    outext = outext//'       RemoveLandTracer = '//temp_str(1)//new_line('a')
    temp_str(1)=self%bathyminNetcdf
    outext = outext//'       bathyminNetcdf = '//temp_str(1)//new_line('a')
    temp_str(1)=self%tracerMaxAge
    outext = outext//'       tracerMaxAge = '//temp_str(1)//new_line('a')
    temp_str(1)=self%Temperature_add_offset
    outext = outext//'       Temperature_add_offset = '//temp_str(1)//""
    call Log%put(outext,.false.)
    end subroutine printsimdefs

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the links of the list
    !---------------------------------------------------------------------------
    subroutine print_stringList(this)
    class(stringList_class), intent(in) :: this
    class(*), pointer :: curr
    call this%reset()               ! reset list iterator
    do while(this%moreValues())     ! loop while there are values to print
        call this%printCurrent()
        call this%next()            ! increment the list iterator
    end do
    call this%reset()               ! reset list iterator
    end subroutine print_stringList

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the current link of the list
    !---------------------------------------------------------------------------
    subroutine print_stringListCurrent(this)
    class(stringList_class), intent(in) :: this
    class(*), pointer :: curr
    type(string) :: outext
    curr => this%currentValue() ! get current value
    select type(curr)
    class is (string)
        outext = curr
        call Log%put(outext, .false.)
        class default
        outext = '[stringList_class::print] Unexepected type of content, not a string'
        call Log%put(outext)
        stop
    end select
    end subroutine print_stringListCurrent

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that checks if an element is already on the list. Returns false if
    !> it is.
    !---------------------------------------------------------------------------
    logical function notRepeated(this, str)
    class(stringList_class), intent(in) :: this
    class(string), intent(in) :: str
    class(*), pointer :: curr
    type(string) :: outext
    notRepeated = .true.
    call this%reset()               ! reset list iterator
    !write(*,*)"entrada do while"
    do while(this%moreValues())     ! loop while there are values to print
        curr => this%currentValue() ! get current value
        select type(curr)
        class is (string)
            if (curr == str) then
                notRepeated = .false.
                return
            end if
            class default
            outext = '[stringList_class::notRepeated] Unexepected type of content, not a string'
            call Log%put(outext)
            stop
        end select
        call this%next()            ! increment the list iterator
    end do
    !write(*,*)"saida do while"
    call this%reset()               ! reset list iterator
    end function notRepeated
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to make the union of two string lists.
    !---------------------------------------------------------------------------
    subroutine makeUnion(this, other)
    class(stringList_class), intent(inout) :: this
    class(stringList_class), intent(inout) :: other
    type(string), allocatable, dimension(:) :: otherList
    integer :: i
    
    call other%toArray(otherList)
    do i = 1, size(otherList)
        if (this%notRepeated(otherList(i))) then
            call this%add(otherList(i))
        end if
    end do
        
    end subroutine makeUnion
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to copy string lists.
    !---------------------------------------------------------------------------
    subroutine copy(this, other)
    class(stringList_class), intent(inout) :: this
    class(stringList_class), intent(inout ) :: other
    type(string), allocatable, dimension(:) :: otherList
    integer :: i
    !write(*,*)"Entrei na copia"
    call other%toArray(otherList)
    call this%finalize()
    !write(*,*)"tamanho da lista = ", size(otherList)
    do i = 1, size(otherList)
        !write(*,*)"variavel adicionada = ", trim(otherList(i))
        call this%add(otherList(i))
    end do
    !write(*,*)"Sai da copia"
    end subroutine copy
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that gets an array with the contents of the list.
    !---------------------------------------------------------------------------
    subroutine toArray(this, currList)
    class(stringList_class), intent(inout) :: this   
    class(*), pointer :: curr
    type(string) :: outext
    type(string), allocatable, dimension(:), intent(inout) :: currList
    integer :: i
    
    i =1
    allocate(currList(this%getSize()))
    
    call this%reset()               ! reset list iterator
    do while(this%moreValues())     ! loop while there are values to print
        curr => this%currentValue() ! get current value
        select type(curr)
        class is (string)
            currList(i) = curr
            i = i + 1
            class default
            outext = '[stringList_class::toArray] Unexepected type of content, not a string'
            call Log%put(outext)
            stop
        end select
        call this%next()            ! increment the list iterator
    end do
    call this%reset()    ! reset list iterator
    
    end subroutine toArray
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that gets a MOHID type waterquality module implementation file.
    !---------------------------------------------------------------------------
    subroutine setWaterQualityFileName(self, fileName)
    implicit none
    class(extImpFiles_t), intent(inout) :: self
    type(string), intent(in) :: fileName
    type(string) :: outext
    integer :: sizem
    logical :: exists
    type(string) :: fileNamePointer
    !Verifies if file exits
    fileNamePointer = fileName
    
    inquire(FILE = trim(adjustl(fileNamePointer%chars())), EXIST = exists)
    if (.not. exists) then
        outext='WaterQuality input file not found. stopping'
        call Log%put(outext)
        stop
    else
        self%waterQualityFileName = fileName
        self%waterQualityFileName_hasValue = .true.
    endif
    
    sizem = sizeof(self%WaterQualityFileName)
    call SimMemory%adddef(sizem)
    sizem = sizeof(self%WaterQualityFileName)
    call SimMemory%adddef(sizem)
    end subroutine setWaterQualityFileName

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that prints external implementation files.
    !---------------------------------------------------------------------------
    subroutine printWaterQualityFileName(self)
    class(extImpFiles_t), intent(in) :: self
    type(string) :: outext
    type(string) :: temp_str(1)

    temp_str(1)=self%WaterQualityFileName%chars()
    outext = '      Water quality file name is '//temp_str(1)//new_line('a')
    
    end subroutine
    
    end module simulationGlobals_mod
