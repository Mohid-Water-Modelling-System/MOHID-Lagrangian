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

    module simulation_globals_mod

    use penf
    use vecfor_r8p
    !use vecfor_r4p
    use stringifor
    use datetime_module
    
    use simulation_precision_mod
    use simulation_parallel_omp_mod
    use simulation_logger_mod
    use simulation_memory_mod
    use utilities_mod

    implicit none
    private

    type :: parameters_t   !< Parameters class
        integer         :: Integrator = 1            !< Integration Algorithm 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)
        integer         :: IntegratorIndexes(3)      !< Index list for the integrator selector
        type(string)    :: IntegratorNames(3)        !< Names list for the integrator selector
        integer         :: numOPMthreads             !< number of openMP threads to be used
        real(prec)      :: WarmUpTime = 0.0          !< Time to freeze the tracers at simulation start (warmup) (s) (default=0.0)
        real(prec)      :: TimeMax = MV              !< Simulation duration (s)
        real(prec)      :: TimeOut = MV              !< Time out data (1/Hz)
        type(datetime)  :: StartTime                 !< Start date of the simulation
        type(datetime)  :: EndTime                   !< End date of the simulation
        type(datetime)  :: BaseDateTime              !< Base date for time stamping results
        integer         :: OutputFormat = 2          !< Format of the output files (default=2) NetCDF=1, VTK=2
        integer         :: InputFormat = 1           !< Format of the Input files (default=1) NetCDF=1, VTK=2
        integer         :: OutputFormatIndexes(2)    !< Index list for the output file format selector
        type(string)    :: OutputFormatNames(2)      !< Names list for the output file format selector
    contains
    procedure :: setparameter
    procedure :: check
    procedure :: print => printsimparameters
    end type parameters_t

    type :: simdefs_t  !< Simulation definitions class
        real(prec)      ::  Dp              !< Initial particle spacing at emission
        real(prec)      ::  dt = MV         !< Timestep for fixed step integrators (s)
        type(vector)    ::  Pointmin        !< Point that defines the lowest corner of the simulation bounding box
        type(vector)    ::  Pointmax        !< Point that defines the upper corner of the simulation bounding box
        logical         ::  autoblocksize = .true.   !< Flag for automatic Block sizing
        type(vector)    ::  blocksize       !< Size (xyz) of a Block (sub-domain)
        integer         ::  numblocks       !< Number of blocks in the simulation
        integer         ::  numblocksx, numblocksy  !<Number of blocks along x and y
    contains
    procedure :: setdp
    procedure :: setdt
    procedure :: setboundingbox
    procedure :: setblocksize
    procedure :: print => printsimdefs
    end type simdefs_t

    type :: constants_t    !< Case Constants class
        type(vector) :: Gravity             !< Gravitational acceleration vector (default=(0 0 -9.81)) (m s-2)
        real(prec)   :: Z0 = 0.0            !< Reference local sea level
        real(prec)   :: Rho_ref = 1000.0    !< Reference density of the medium (default=1000.0) (kg m-3)
    contains
    procedure :: setgravity
    procedure :: setz0
    procedure :: setrho
    procedure :: print => printconstants
    end type constants_t

    type :: filenames_t    !<File names class
        type(string) :: mainxmlfilename     !< Input .xml file name
        type(string) :: propsxmlfilename    !< Properties .xml file name
        type(string) :: tempfilename        !< Generic temporary file name
        type(string) :: outpath             !< General output directory
        type(string) :: casename            !< Name of the running case
    end type filenames_t

    type src_parm_t   !<Lists for Source parameters
        type(string), allocatable, dimension(:) :: baselist !<Lists for base tracer parameters
        type(string), allocatable, dimension(:) :: particulatelist  !<List for parameters of particulate type tracers
    contains
    procedure :: buildlists
    end type

    type :: sim_t  !<Simulation related counters and others
        private
        integer :: numdt        !<number of the current iteration
        integer :: numoutfile   !<number of the current output file
        integer :: numTracer    !<Global Tracer number holder. Incremented at tracer construction or first activation time
    contains
    procedure, public :: increment_numdt
    procedure, public :: increment_numoutfile
    procedure, public :: getnumdt
    procedure, public :: getnumoutfile    
    procedure, public :: getnumTracer
    procedure, private :: increment_numTracer
    end type sim_t
    
    type :: var_names_t
        type(string) :: u
        type(string) :: v
        type(string) :: w
        type(string) :: temp
        type(string) :: sal
        type(string) :: density
        type(string) :: lon
        type(string) :: lat
        type(string) :: depth
        type(string) :: time
    contains
    procedure, private :: buildvars
    end type var_names_t
    
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
    !procedure :: getTimeStamp
    end type sim_time_t

    type :: globals_class   !<Globals class - This is a container for every global variable on the simulation
        type(parameters_t)  :: Parameters
        type(simdefs_t)     :: SimDefs
        type(constants_t)   :: Constants
        type(filenames_t)   :: Names
        type(src_parm_t)    :: SrcProp
        type(sim_t)         :: Sim
        type(var_names_t)   :: Var
        type(sim_time_t)    :: SimTime
    contains
    procedure :: initialize => setdefaults
    procedure :: setTimeDate
    end type globals_class

    !Simulation variables
    type(globals_class) :: Globals

    !Public access vars
    public :: Globals

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Globals default setting routine.
    !> @param[in] self, outpath
    !---------------------------------------------------------------------------
    subroutine setdefaults(self, outpath)
    implicit none
    class(globals_class), intent(inout) :: self
    integer :: sizem
    type(string), optional, intent(in) :: outpath
    !parameters
    self%Parameters%Integrator = 1
    self%Parameters%IntegratorIndexes = [1,2,3]
    self%Parameters%IntegratorNames(1) = 'Euler'
    self%Parameters%IntegratorNames(2) = 'Multi-Step Euler'
    self%Parameters%IntegratorNames(3) = 'Runge-Kuta 4'
    self%Parameters%numOPMthreads = OMPManager%getThreads()
    self%Parameters%WarmUpTime = 0.0
    self%Parameters%TimeOut = MV
    self%Parameters%StartTime = datetime()
    self%Parameters%EndTime = datetime()
    self%Parameters%BaseDateTime = datetime(1950,1,1,0,0,0)
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
    self%SimDefs%numblocks = 16  !placeholder number, should be numThreads or numProcesses or computed by user dimensions
    self%SimDefs%Dp = MV
    self%SimDefs%dt = MV
    self%SimDefs%Pointmin = 0.0
    self%SimDefs%Pointmax = 0.0
    !simulation constants
    self%Constants%Gravity= 0.0*ex + 0.0*ey -9.81*ez
    self%Constants%Z0 = 0.0
    self%Constants%Rho_ref = 1000.0
    !filenames
    self%Names%mainxmlfilename = 'not_set'
    self%Names%propsxmlfilename = 'not_set'
    self%Names%tempfilename = 'not_set'
    if (present(outpath)) then
        self%Names%outpath = outpath
    else
        self%Names%outpath = 'not_set'
    end if
    self%Names%casename = 'not_set'
    !global time
    self%SimTime%CurrTime = 0.0
    !global counters
    self%Sim%numdt = 0
    self%Sim%numoutfile = 0
    self%Sim%numTracer = 0
    !Source parameters list
    call self%SrcProp%buildlists()
    !Variable names
    call self%Var%buildvars()

    sizem=sizeof(self)
    call SimMemory%adddef(sizem)

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Builds variable list names.
    !---------------------------------------------------------------------------
    subroutine buildvars(self)
    implicit none
    class(var_names_t), intent(inout) :: self
    self%u       = 'vel_x'
    self%v       = 'vel_y'
    self%w       = 'vel_z'
    self%temp    = 'temp'
    self%sal     = 'sal'
    self%density = 'density'
    self%lon     = 'lon'
    self%lat     = 'lat'
    self%depth   = 'depth'
    self%time    = 'time'
    end subroutine buildvars
        
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> initializes time stamp and date stamp variables
    !---------------------------------------------------------------------------
    subroutine setTimeDate(self)
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
    class(sim_t), intent(inout) :: self
    self%numdt = self%numdt + 1
    end subroutine increment_numdt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the number of time steps.
    !---------------------------------------------------------------------------
    integer function getnumdt(self)
    implicit none
    class(sim_t), intent(inout) :: self
    getnumdt = self%numdt
    end function getnumdt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> incrementing output file count.
    !---------------------------------------------------------------------------
    subroutine increment_numoutfile(self)
    implicit none
    class(sim_t), intent(inout) :: self
    self%numoutfile = self%numoutfile + 1
    end subroutine increment_numoutfile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the number of output files written.
    !---------------------------------------------------------------------------
    integer function getnumoutfile(self)
    implicit none
    class(sim_t), intent(inout) :: self
    getnumoutfile = self%numoutfile
    end function getnumoutfile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to build the parameters list of the Sources
    !---------------------------------------------------------------------------
    subroutine buildlists(self)
    implicit none
    class(src_parm_t), intent(inout) :: self
    allocate(self%baselist(5))
    self%baselist(1) = 'particulate'
    self%baselist(2) = 'density'
    self%baselist(3) = 'radius'
    self%baselist(4) = 'condition'
    self%baselist(5) = 'degradation_rate'
    allocate(self%particulatelist(2))
    self%particulatelist(1) = 'intitial_concentration'
    self%particulatelist(2) = 'particle_radius'
    end subroutine buildlists

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Private parameter setting method. Builds the simulation parametric space 
    !> from the input case file.
    !> @param[in] self, parmkey, parmvalue
    !---------------------------------------------------------------------------
    subroutine setparameter(self,parmkey,parmvalue)
    implicit none
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
            self%numOPMthreads=parmvalue%to_number(kind=1._R4P)
            call OMPManager%setThreads(self%numOPMthreads)
        end if
        sizem=sizeof(self%WarmUpTime)
    elseif(parmkey%chars()=="WarmUpTime") then
        self%WarmUpTime=parmvalue%to_number(kind=1._R4P)
        sizem=sizeof(self%WarmUpTime)
    elseif(parmkey%chars()=="TimeOut") then
        self%TimeOut=parmvalue%to_number(kind=1._R4P)
        sizem=sizeof(self%TimeOut)
    elseif(parmkey%chars()=="StartTime") then
        date = Utils%getDateFromISOString(parmvalue)
        self%StartTime = datetime(date(1),date(2),date(3),date(4),date(5),date(6))
        if (.not. self%StartTime%isValid()) self%StartTime = datetime()        
        sizem=sizeof(self%StartTime)        
    elseif(parmkey%chars()=="EndTime") then
        date = Utils%getDateFromISOString(parmvalue)
        self%EndTime = datetime(date(1),date(2),date(3),date(4),date(5),date(6))
        if (.not. self%EndTime%isValid()) self%EndTime = datetime()        
        sizem=sizeof(self%EndTime)
    elseif(parmkey%chars()=="BaseDateTime") then
        date = Utils%getDateFromISOString(parmvalue)
        self%BaseDateTime = datetime(date(1),date(2),date(3),date(4),date(5),date(6))
        if (.not. self%BaseDateTime%isValid()) self%BaseDateTime = datetime()        
        sizem=sizeof(self%BaseDateTime)
    elseif(parmkey%chars()=="OutputFormat") then
        self%OutputFormat=parmvalue%to_number(kind=1_I1P)
        sizem=sizeof(self%OutputFormat)
    end if
    call SimMemory%adddef(sizem)

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Parameter checking method. Checks if mandatory parameters were set
    !---------------------------------------------------------------------------
    subroutine check(self)
    implicit none
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
    if (self%TimeOut==MV) then
        outext = '[Globals::parameters::check]: sampling rate parameter (TimeOut) is not set, stoping'
        call Log%put(outext)
        stop
    elseif (self%StartTime==temp) then
        outext = '[Globals::parameters::check]: start time parameter (StartTime) is not set or invalid, stoping'
        call Log%put(outext)
        stop
    elseif (self%EndTime==temp) then
        outext = '[Globals::parameters::check]: end time parameter (EndTime) is not set or invalid, stoping'
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
    temp_str=self%TimeOut
    outext = outext//'       TimeOut = '//temp_str//' Hz'//new_line('a')
    outext = outext//'       Output file format is '//self%OutputFormatNames(self%OutputFormat)
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
    self%Z0=read_z0%to_number(kind=1._R4P)
    sizem = sizeof(self%Z0)
    call SimMemory%adddef(sizem)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Rho_Ref setting routine.
    !> @param[in] self, read_rho
    !---------------------------------------------------------------------------
    subroutine setrho(self,read_rho)
    implicit none
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_rho
    type(string) :: outext
    integer :: sizem
    self%Rho_ref=read_rho%to_number(kind=1._R4P)
    if (self%Rho_ref.le.0.0) then
        outext='Rho_ref must be positive and non-zero, stopping'
        call Log%put(outext)
        stop
    endif
    sizem = sizeof(self%Rho_ref)
    call SimMemory%adddef(sizem)
    end subroutine

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
    temp_str(1)=self%Rho_ref
    outext = outext//'       Rho_ref = '//temp_str(1)//' kg/m^3'

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
    type(string), intent(in) :: read_dp
    type(string) :: outext
    integer :: sizem
    self%Dp=read_dp%to_number(kind=1._R4P)
    if (self%Dp.le.0.0) then
        outext='Dp must be positive and non-zero, stopping'
        call Log%put(outext)
        stop
    endif
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
    self%dt=read_dt%to_number(kind=1._R4P)
    if (self%dt.le.0.0) then
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
    if (point_%chars() == "pointmin") then
        self%Pointmin= coords
    elseif (point_%chars() == "pointmax") then
        self%Pointmax= coords
    endif
    sizem=sizeof(coords)
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
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Public simulation definitions printing routine.
    !---------------------------------------------------------------------------
    subroutine printsimdefs(self)
    implicit none
    class(simdefs_t), intent(in) :: self
    type(string) :: outext
    type(string) :: temp_str(3)

    temp_str(1)=self%Dp
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
        outext = outext//'       Blocks are automatically sized'
    else
        temp_str(1)=self%blocksize%x
        temp_str(2)=self%blocksize%y
        outext = outext//'       Blocks are sized '//new_line('a')//&
            '       '//temp_str(1)//' X '//temp_str(2)
    end if
    call Log%put(outext,.false.)
    end subroutine printsimdefs

    end module simulation_globals_mod
