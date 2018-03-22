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

    module simulation_globals

    use commom_modules

    implicit none
    private

    type parameters_t   !< Parameters class
        integer    :: Integrator = 1        !< Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
        real(prec) :: CFL = 0.5             !< Courant Friedrichs Lewy condition number
        real(prec) :: WarmUpTime = 0.0      !< Time to freeze the tracers at simulation start (warmup) (s) (default=0.0)
        real(prec) :: TimeMax = MV          !< Simulation duration (s)
        real(prec) :: TimeOut = MV          !< Time out data (1/Hz)
    contains
    procedure :: setparameter
    procedure :: printout => printsimparameters
    end type
    
    type simdefs_t  !< Simulation definitions class
        real(prec) :: Dp = MV               !< Initial particle spacing at source generation
        real(prec_time) :: dt               !< Timestep for fixed step integrators (s)
        type(vector)    ::  Pointmin        !< Point that defines the lowest corner of the simulation bounding box
        type(vector)    ::  Pointmax        !< Point that defines the upper corner of the simulation bounding box
    contains
    procedure :: setdp
    procedure :: setdt
    procedure :: setboundingbox
    procedure :: printout => printsimdefs
    end type
    
    type constants_t    !< Case Constants class
        type(vector) :: Gravity             !< Gravitational acceleration vector (default=(0 0 -9.81)) (m s-2)
        real(prec) :: Rho_ref = 1000.0      !< Reference density of the medium (default=1000.0) (kg m-3)
    contains
    procedure :: setgravity
    procedure :: setrho
    end type
    
    real(prec_time) :: SimTime
       
    !Simulation variables
    type(parameters_t)  :: Parameters
    type(simdefs_t)     :: SimDefs
    type(constants_t)   :: Constants
    
    !Public access vars
    public :: SimTime, Parameters, SimDefs, Constants

    contains
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private parameter setting routine. Builds the simulation parametric space from the input case file.
    !
    !> @param[in] parmkey, parmvalue
    !---------------------------------------------------------------------------
    subroutine setparameter(self,parmkey,parmvalue)
    implicit none
    class(parameters_t), intent(inout) :: self
    type(string), intent(in) :: parmkey
    type(string), intent(in) :: parmvalue
    character(80) :: value
    integer :: size
    !add new parameters to this search
    if (parmkey%chars()=="Integrator") then
        self%Integrator=parmvalue%to_number(kind=1_I1P)
        size=sizeof(self%Integrator)
    elseif(parmkey%chars()=="CFL") then
        self%CFL=parmvalue%to_number(kind=1._R4P)
        size=sizeof(self%CFL)
    elseif(parmkey%chars()=="WarmUpTime") then
        self%WarmUpTime=parmvalue%to_number(kind=1._R4P)
        size=sizeof(self%WarmUpTime)
    elseif(parmkey%chars()=="TimeMax") then
        self%TimeMax=parmvalue%to_number(kind=1._R4P)
        size=sizeof(self%TimeMax)
    elseif(parmkey%chars()=="TimeOut") then
        self%TimeOut=parmvalue%to_number(kind=1._R4P)
        size=sizeof(self%TimeOut)        
    endif
    call SimMemory%adddef(size)
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private parameter printing routine.
    !---------------------------------------------------------------------------
    subroutine printsimparameters(self)
    implicit none
    class(parameters_t), intent(inout) :: self
    type(string) :: outext
    type(string) :: temp_str
    call getintegratorname(temp_str,self%Integrator)
    outext = '      Integrator scheme is '//temp_str//new_line('a')
    temp_str=self%CFL
    outext = outext//'       CFL='//temp_str//new_line('a')
    temp_str=self%WarmUpTime
    outext = outext//'       WarmUpTime='//temp_str//' s'//new_line('a')
    temp_str=self%TimeMax
    outext = outext//'       TimeMax='//temp_str//' s'//new_line('a')
    temp_str=self%TimeOut
    outext = outext//'       TimeOut='//temp_str//' Hz'
    call ToLog(outext,.false.)
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine to get integrator scheme name
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
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public Gravity setting routine.
    !
    !> @param[in] grav
    !---------------------------------------------------------------------------
    subroutine setgravity (self,grav)
    implicit none
    class(constants_t), intent(inout) :: self
    type(vector) :: grav
    integer :: size
    self%Gravity= grav
    size=sizeof(self%Gravity)
    call SimMemory%adddef(size)
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Provate Rho_Ref setting routine.
    !
    !> @param[in] read_rho
    !---------------------------------------------------------------------------
    subroutine setrho(self,read_rho)
    implicit none
    class(constants_t), intent(inout) :: self
    type(string), intent(in) :: read_rho
    type(string) :: outext
    integer :: size
    self%Rho_ref=read_rho%to_number(kind=1._R4P)
    if (self%Rho_ref.le.0.0) then
        outext='Rho_ref must be positive and non-zero, stopping'
        call ToLog(outext)
        stop
    endif
    size = sizeof(self%Rho_ref)
    call SimMemory%adddef(size)
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private dp setting routine.
    !
    !> @param[in] read_dp
    !---------------------------------------------------------------------------
    subroutine setdp(self,read_dp)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_dp
    type(string) :: outext
    integer :: size
    self%Dp=read_dp%to_number(kind=1._R4P)
    if (self%Dp.le.0.0) then
        outext='Dp must be positive and non-zero, stopping'
        call ToLog(outext)
        stop
    endif
    size = sizeof(self%Dp)
    call SimMemory%adddef(size)
    return
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private dt setting routine.
    !
    !> @param[in] read_dt
    !---------------------------------------------------------------------------
    subroutine setdt(self,read_dt)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: read_dt
    type(string) :: outext
    integer :: size
    self%dt=read_dt%to_number(kind=1._R4P)
    if (self%dt.le.0.0) then
        outext='dt must be positive and non-zero, stopping'
        call ToLog(outext)
        stop
    endif
    size = sizeof(self%dt)
    call SimMemory%adddef(size)
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private bounding box setting routine.
    !
    !> @param[in] point_, coords
    !---------------------------------------------------------------------------
    subroutine setboundingbox(self,point_, coords)
    implicit none
    class(simdefs_t), intent(inout) :: self
    type(string), intent(in) :: point_
    type(vector) :: coords
    integer :: size
    if (point_%chars() == "pointmin") then
        self%Pointmin= coords
    elseif (point_%chars() == "pointmax") then
        self%Pointmax= coords
    endif
    size=sizeof(coords)
    call SimMemory%adddef(size)
    return
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
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
        '       '//temp_str(1)//new_line('a')//&
        '       '//temp_str(2)//new_line('a')//&
        '       '//temp_str(3)//new_line('a')
    temp_str(1)=self%Pointmax%x
    temp_str(2)=self%Pointmax%y
    temp_str(3)=self%Pointmax%z
    outext = outext//'       Pointmax (BB) is '//new_line('a')//&
        '       '//temp_str(1)//new_line('a')//&
        '       '//temp_str(2)//new_line('a')//&
        '       '//temp_str(3)

    call ToLog(outext,.false.)
    end subroutine

    end module simulation_globals
