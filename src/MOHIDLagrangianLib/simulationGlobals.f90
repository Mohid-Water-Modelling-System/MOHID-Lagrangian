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

    implicit none
    private

    type, extends(linkedlist) :: stringList_class !< List of strings class
    contains
    procedure :: print => print_stringList
    procedure :: printCurrent => print_stringListCurrent
    procedure :: notRepeated
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
        integer         :: VerticalVelMethod !< Vertical velocity method
        integer         :: RemoveLandTracer !< Vertical velocity method
    contains
    procedure :: setdp
    procedure :: setdt
    procedure :: setboundingbox
    procedure :: setCenter
    procedure :: setblocksize
    procedure :: setVerticalVelMethod
    procedure :: setRemoveLandTracer
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
        real(prec)   :: MeanDensity         !< mean ocean water density.
        real(prec)   :: MeanKVisco          !< mean ocean water kinematic viscosity
    contains
    procedure :: setgravity
    procedure :: setz0
    procedure :: setBeachingLevel
    procedure :: setBeachingStopProb
    procedure :: setDiffusionCoeff
    procedure :: setSmallDt
    procedure :: setResuspensionCoeff
    procedure :: setMeanDensity
    procedure :: setMeanKVisco    
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
        type(string) :: temp
        type(string) :: sal
        type(string) :: density
        type(string) :: vsdx
        type(string) :: vsdy
        type(string) :: u10
        type(string) :: v10
        type(string) :: lon
        type(string) :: lat
        type(string) :: level
        type(string) :: time
        type(string) :: landIntMask
        type(string) :: resolution
        type(string) :: bathymetry
        type(string) :: surface
        type(string) :: rate
        type(stringList_class) :: uVariants !< possible names for 'u' in the input files
        type(stringList_class) :: vVariants
        type(stringList_class) :: wVariants
        type(stringList_class) :: tempVariants
        type(stringList_class) :: salVariants
        type(stringList_class) :: densityVariants
        type(stringList_class) :: vsdxVariants
        type(stringList_class) :: vsdyVariants
        type(stringList_class) :: u10Variants
        type(stringList_class) :: v10Variants
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
    end type var_names_t

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
    contains    
    end type sources_t

    type :: globals_class   !<Globals class - This is a container for every global variable on the simulation
        type(parameters_t)  :: Parameters
        type(simdefs_t)     :: SimDefs
        type(constants_t)   :: Constants
        type(filenames_t)   :: Names
        type(sim_t)         :: Sim
        type(output_t)      :: Output
        type(var_names_t)   :: Var
        type(sim_time_t)    :: SimTime
        type(maskVals_t)    :: Mask
        type(tracerTypes_t) :: Types
        type(dataTypes_t)   :: DataTypes
        type(sources_t)     :: Sources
    contains
    procedure :: initialize => setdefaults
    procedure :: setTimeDate
    procedure, public :: setNamingConventions
    procedure, public :: setInputFileNames
    procedure :: setVarNames
    procedure :: setDimNames
    procedure :: setCurrVar
    end type globals_class
    
    type(string) :: notRead
    type(string) :: notSet

    !Simulation variables
    type(globals_class) :: Globals

    !Public access vars
    public :: Globals, stringList_class, notRead, notSet

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
    self%SimDefs%RemoveLandTracer = 1
    !simulation constants
    self%Constants%Gravity= 0.0*ex + 0.0*ey -9.81*ez
    self%Constants%Z0 = 0.0
    self%Constants%BeachingLevel = -3.0
    self%Constants%BeachingStopProb = 0.50
    self%Constants%DiffusionCoeff = 1.0
    self%Constants%smallDt = 0.0
    self%Constants%ResuspensionCoeff = 0.0
     self%Constants%MeanDensity = 1027.0
     self%Constants%MeanKVisco = 1.09E-3
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
    self%temp    = 'temp'
    self%sal     = 'sal'
    self%density = 'density'
    self%vsdx    = 'vsdx'
    self%vsdy    = 'vsdy'
    self%u10     = 'u10'
    self%v10     = 'v10'
    self%lon     = 'lon'
    self%lat     = 'lat'
    self%level   = 'level'
    self%time    = 'time'
    self%landIntMask = 'landIntMask'
    self%resolution = 'resolution'
    self%bathymetry = 'bathymetry'
    self%surface    = 'surface'
    self%rate = 'rate'
    !adding variables to variable pool - PLACEHOLDER, this should come from tracer constructors
    call self%addVar(self%u)
    call self%addVar(self%v)
    call self%addVar(self%w)
    call self%addVar(self%temp)
    call self%addVar(self%sal)
    call self%addVar(self%density)
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

    tag="eastward_sea_water_velocity"
    call self%setCurrVar(tag, self%Var%u, self%Var%uVariants, varNode)
    tag="northward_sea_water_velocity"
    call self%setCurrVar(tag, self%Var%v, self%Var%vVariants, varNode)
    tag="upward_sea_water_velocity"
    call self%setCurrVar(tag, self%Var%w, self%Var%wVariants, varNode)
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
    tag="eastward_wind"
    call self%setCurrVar(tag, self%Var%u10, self%Var%u10Variants, varNode)
    tag="northward_wind"
    call self%setCurrVar(tag, self%Var%v10, self%Var%v10Variants, varNode)

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

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Resuspension setting routine or not
    !> @param[in] self, read_BeachingLevel
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
    temp_str(1)=self%ResuspensionCoeff
    outext = outext//'       ResuspensionCoeff = '//temp_str(1)//new_line('a')
    temp_str(1)=self%MeanDensity
    outext = outext//'       MeanDensity = '//temp_str(1)//new_line('a')
    temp_str(1)=self%MeanKVisco
    outext = outext//'       MeanKVisco = '//temp_str(1)//''
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
    !> @param[in] self, read_BeachingLevel
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
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Resuspension setting routine or not
    !> @param[in] self, read_BeachingLevel
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
    outext = outext//'       RemoveLandTracer = '//temp_str(1)//''
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
    call this%reset()               ! reset list iterator
    end function notRepeated

    end module simulationGlobals_mod
