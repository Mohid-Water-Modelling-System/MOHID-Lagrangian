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
    use vecfor_r4p !Should include a preprocessor switch
    !use vecfor
    use stringifor
    use datetime_module

    use simulation_precision_mod
    use simulation_logger_mod
    use simulation_memory_mod

    implicit none
    private

    type parameters_t   !< Parameters class
        integer    :: Integrator = 1            !< Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
        real(prec) :: CFL = 0.5                 !< Courant Friedrichs Lewy condition number
        real(prec_time) :: WarmUpTime = 0.0     !< Time to freeze the tracers at simulation start (warmup) (s) (default=0.0)
        real(prec_time) :: TimeMax = MV         !< Simulation duration (s)
        real(prec) :: TimeOut = MV              !< Time out data (1/Hz)
        type(datetime) :: StartTime
        type(datetime) :: EndTime
    contains
    procedure :: setparameter
    procedure :: check
    procedure :: print => printsimparameters
    end type

    type simdefs_t  !< Simulation definitions class
        real(prec)      ::  Dp = MV         !< Initial particle spacing at source generation
        real(prec_time) ::  dt = MV         !< Timestep for fixed step integrators (s)
        type(vector)    ::  Pointmin        !< Point that defines the lowest corner of the simulation bounding box
        type(vector)    ::  Pointmax        !< Point that defines the upper corner of the simulation bounding box
        logical         ::  autoblocksize = .true.   !< Flag for automatic Block sizing
        type(vector)    ::  blocksize       !< Size (width & heigth) of a Block (sub-domain)
        integer         ::  numblocks       !< Number of blocks in the simulation
    contains
    procedure :: setdp
    procedure :: setdt
    procedure :: setboundingbox
    procedure :: setblocksize
    procedure :: print => printsimdefs
    end type

    type constants_t    !< Case Constants class
        type(vector) :: Gravity             !< Gravitational acceleration vector (default=(0 0 -9.81)) (m s-2)
        real(prec)   :: Z0 = 0.0            !< Reference local sea level
        real(prec)   :: Rho_ref = 1000.0    !< Reference density of the medium (default=1000.0) (kg m-3)
    contains
    procedure :: setgravity
    procedure :: setz0
    procedure :: setrho
    procedure :: print => printconstants
    end type

    type filenames_t    !<File names class
        type(string) :: mainxmlfilename     !< Input .xml file name
        type(string) :: propsxmlfilename    !< Properties .xml file name
        type(string) :: tempfilename        !< Generic temporary file name
    end type

    type globals_class   !<Globals class - This is a container for every global variable on the simulation
        type(parameters_t)  :: Parameters
        type(simdefs_t)     :: SimDefs
        type(constants_t)   :: Constants
        type(filenames_t)   :: FileNames
        real(prec_time)     :: SimTime
      contains
      procedure :: initialize => setdefaults
    end type

    !Simulation variables
    type(globals_class) :: Globals

    !Public access vars
    public :: Globals

    contains

      !---------------------------------------------------------------------------
      !> @author Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Globals default setting routine.
      !---------------------------------------------------------------------------
      subroutine setdefaults(self)
      implicit none
      class(globals_class), intent(inout) :: self
      integer :: sizem
      !parameters
      self%Parameters%Integrator = 1
      self%Parameters%CFL = 0.5
      self%Parameters%WarmUpTime = 0.0
      self%Parameters%TimeOut = MV
      self%Parameters%TimeOut = MV
      self%Parameters%StartTime = datetime()
      self%Parameters%EndTime = datetime()
      !Simulation definitions
      self%SimDefs%autoblocksize =.true.
      self%SimDefs%blocksize = 0.0
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
      self%FileNames%mainxmlfilename = 'not_set'
      self%FileNames%propsxmlfilename = 'not_set'
      self%FileNames%tempfilename = 'not_set'
      !global time
      self%SimTime = 0.0

      sizem=sizeof(self)
      call SimMemory%adddef(sizem)

      end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private parameter setting method. Builds the simulation parametric space from the input case file.
    !
    !> @param[in] parmkey, parmvalue
    !---------------------------------------------------------------------------
    subroutine setparameter(self,parmkey,parmvalue)
    implicit none
    class(parameters_t), intent(inout) :: self
    type(string), intent(in) :: parmkey
    type(string), intent(in) :: parmvalue
    type(string), allocatable :: dc(:)
    integer :: i, date(6)
    integer :: sizem
    !add new parameters to this search
    if (parmkey%chars()=="Integrator") then
        self%Integrator=parmvalue%to_number(kind=1_I1P)
        sizem=sizeof(self%Integrator)
    elseif(parmkey%chars()=="CFL") then
        self%CFL=parmvalue%to_number(kind=1._R4P)
        sizem=sizeof(self%CFL)
    elseif(parmkey%chars()=="WarmUpTime") then
        self%WarmUpTime=parmvalue%to_number(kind=1._R4P)
        sizem=sizeof(self%WarmUpTime)
    elseif(parmkey%chars()=="TimeOut") then
        self%TimeOut=parmvalue%to_number(kind=1._R4P)
        sizem=sizeof(self%TimeOut)
    elseif(parmkey%chars()=="StartTime") then
        call parmvalue%split(tokens=dc, sep=' ')
        if (size(dc) == 6) then
            do i=1, size(dc)
                date(i) = dc(i)%to_number(kind=1._R4P)
            end do
            self%StartTime = datetime(date(1),date(2),date(3),date(4),date(5),date(6))
            if (self%StartTime%isValid()) then
            else
                self%StartTime = datetime() !reseting to default so it is caught later on
            end if
            sizem=sizeof(self%StartTime)
        else
            stop '[Globals::setparameter] StartTime parameter not in correct format. Eg. "2009 3 1 0 0 0"'
        end if
    elseif(parmkey%chars()=="EndTime") then
        call parmvalue%split(tokens=dc, sep=' ')
        if (size(dc) == 6) then
            do i=1, size(dc)
                date(i) = dc(i)%to_number(kind=1._R4P)
            end do
            self%EndTime = datetime(date(1),date(2),date(3),date(4),date(5),date(6))
            if (self%EndTime%isValid()) then
            else
                self%EndTime = datetime() !reseting to default so it is caught later on
            end if
            sizem=sizeof(self%EndTime)
        else
            stop '[Globals::setparameter] EndTime parameter not in correct format. Eg. "2009 3 1 0 0 0"'
        end if
    endif
    call SimMemory%adddef(sizem)

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Parameter checking method. Checks if mandatory parameters were set
    !---------------------------------------------------------------------------
    subroutine check(self)
    implicit none
    class(parameters_t), intent(inout) :: self
    type(string) :: outext
    type(datetime) :: temp
    type(timedelta) :: simtime

    temp = datetime() !default initialization
    !add new parameters to this search
    if (self%TimeOut==MV) then
        outext = 'Simulation sampling rate parameter (TimeOut) is not set, stoping'
        call Log%put(outext)
        stop
    elseif (self%StartTime==temp) then
        outext = 'Simulation start time parameter (StartTime) is not set or invalid, stoping'
        call Log%put(outext)
        stop
    elseif (self%EndTime==temp) then
        outext = 'Simulation end time parameter (EndTime) is not set or invalid, stoping'
        call Log%put(outext)
        stop
    endif
    !Build timemax from the difference between start and end time
    simtime = self%EndTime - self%StartTime
    self%TimeMax = simtime%total_seconds()
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Parameter printing method.
    !---------------------------------------------------------------------------
    subroutine printsimparameters(self)
    implicit none
    class(parameters_t), intent(inout) :: self
    type(string) :: outext
    type(string) :: temp_str
    character(len=23) :: temp_char
    call getintegratorname(temp_str,self%Integrator)
    outext = '      Integrator scheme is '//temp_str//new_line('a')
    temp_str=self%CFL
    outext = outext//'       CFL = '//temp_str//new_line('a')
    temp_str=self%WarmUpTime
    outext = outext//'       WarmUpTime = '//temp_str//' s'//new_line('a')
    temp_str=self%TimeOut
    outext = outext//'       TimeOut = '//temp_str//' Hz'//new_line('a')
    temp_char = self%StartTime%isoformat(' ')
    temp_str = temp_char
    outext = outext//'       StartTime = '//temp_str//new_line('a')
    temp_char = self%EndTime%isoformat(' ')
    temp_str = temp_char
    outext = outext//'       EndTime   = '//temp_str//new_line('a')
    temp_str=self%TimeMax
    outext = outext//'       Simulation will run for '//temp_str//' s'
    call Log%put(outext,.false.)
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Routine to get integrator scheme name
    !---------------------------------------------------------------------------
    subroutine getintegratorname(name,code)
    implicit none
    type(string), intent(inout) :: name
    integer, intent(in) :: code
    if (code==1) then
        name='Verlet'
    elseif(code==2)then
        name='Symplectic'
    elseif(code==3)then
        name='Runge-Kuta 4'
    endif
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Gravity setting routine.
    !
    !> @param[in] grav
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Z0 setting routine.
    !
    !> @param[in] read_z0
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Rho_Ref setting routine.
    !
    !> @param[in] read_rho
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
    ! Routine Author Name and Affiliation.
    !
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Dp setting routine.
    !
    !> @param[in] read_dp
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Dt setting routine.
    !
    !> @param[in] read_dt
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Bounding box setting routine.
    !
    !> @param[in] point_, coords
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> blocksize box setting routine.
    !
    !> @param[in] bsize
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
    ! Routine Author Name and Affiliation.
    !
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
    temp_str(1)=self%blocksize%x
    temp_str(2)=self%blocksize%y
    outext = outext//'       Blocks are sized '//new_line('a')//&
        '       '//temp_str(1)//' X '//temp_str(2)

    call Log%put(outext,.false.)
    end subroutine

  end module simulation_globals_mod
