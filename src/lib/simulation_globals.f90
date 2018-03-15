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
    !> Module to finalize the simulation. This presents a public routine that is in charge of deallocating all global variables, closing all files and print some simulation-related statistics.
    !------------------------------------------------------------------------------

    module simulation_globals

    use commom_modules

    implicit none
    private

    !Public access vars
    public :: Integrator, CFL
    public :: WarmUpTime
    public :: TimeMax, TimeOut
    public :: Dp, Pointmin, Pointmax
    public :: Gravity
    public :: Rho_ref
    
    !Public access procedures
    public :: setSimParameter, setSimGravity, setSimRho, setSimDp, setSimBounds
    public :: printSimParameters, printSimDefs

    !Parameters
    integer    :: Integrator = 1        !< Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
    real(prec) :: CFL = 0.5             !< Courant Friedrichs Lewy condition number
    real(prec) :: WarmUpTime = 0.0      !< Time to freeze the tracers at simulation start (warmup) (s) (default=0.0)
    real(prec) :: TimeMax = MV          !< Simulation duration (s)
    real(prec) :: TimeOut = MV          !< Time out data (1/Hz)
    !Simulation definitions
    real(prec) :: Dp = MV               !< Initial particle spacing at source generation
    type(vector)    ::  Pointmin        !< Point that defines the lowest corner of the simulation bounding box
    type(vector)    ::  Pointmax        !< Point that defines the upper corner of the simulation bounding box
    !Case Constants
    type(vector) :: Gravity             !< Gravitational acceleration vector (default=(0 0 -9.81)) (m s-2)
    real(prec) :: Rho_ref = 1000.0      !< Reference density of the medium (default=1000.0) (kg m-3)
    
    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public parameter setting routine. Builds the simulation parametric space from the input xml case file.
    !
    !> @param[in] parmkey, parmvalue
    !---------------------------------------------------------------------------
    subroutine setSimParameter(parmkey,parmvalue)
    implicit none
    type(string), intent(in) :: parmkey
    type(string), intent(in) :: parmvalue
    character(80) :: value
    !add new parameters to this search
    if (parmkey%chars()=="Integrator") then
        Integrator=parmvalue%to_number(kind=1_I1P)        
    elseif(parmkey%chars()=="CFL") then
        CFL=parmvalue%to_number(kind=1._R4P)
    elseif(parmkey%chars()=="WarmUpTime") then
        WarmUpTime=parmvalue%to_number(kind=1._R4P)
    elseif(parmkey%chars()=="TimeMax") then
        TimeMax=parmvalue%to_number(kind=1._R4P)
    elseif(parmkey%chars()=="TimeOut") then
        TimeOut=parmvalue%to_number(kind=1._R4P)
    endif
    return
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
    subroutine setSimGravity (grav)
    implicit none    
    type(vector) :: grav
    Gravity= grav
    return
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public Rho_Ref setting routine.
    !
    !> @param[in] read_rho
    !---------------------------------------------------------------------------
    subroutine setSimRho(read_rho)
    implicit none
    type(string), intent(in) :: read_rho  
    type(string) :: outext
    Rho_ref=read_rho%to_number(kind=1._R4P)
    if (Rho_ref.le.0.0) then
        outext='Rho_ref must be positive and non-zero, stopping'
        call ToLog(outext)
        stop
    endif
    return
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public dp setting routine.
    !
    !> @param[in] read_dp
    !---------------------------------------------------------------------------
    subroutine setSimDp(read_dp)
    implicit none
    type(string), intent(in) :: read_dp
    type(string) :: outext
    Dp=read_dp%to_number(kind=1._R4P)
    if (Dp.le.0.0) then
        outext='Dp must be positive and non-zero, stopping'
        call ToLog(outext)
        stop
    endif
    return
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public bounding box setting routine.
    !
    !> @param[in] point_, coords 
    !---------------------------------------------------------------------------
    subroutine setSimBounds(point_, coords)
    implicit none
    type(string), intent(in) :: point_
    type(vector) :: coords
    if (point_%chars() == "pointmin") then
        Pointmin= coords
    elseif (point_%chars() == "pointmax") then
        Pointmax= coords
    endif    
    return
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public simulation definitions printing routine. 
    !---------------------------------------------------------------------------
    subroutine printSimDefs
    implicit none    
    type(string) :: outext    
    type(string) :: temp_str(3)
      
    temp_str(1)=Dp
    outext = '      Initial resolution is '//temp_str(1)//' m'//new_line('a')
    temp_str(1)=Pointmin%x
    temp_str(2)=Pointmin%y
    temp_str(3)=Pointmin%z    
    outext = outext//'       Pointmin (BB) is '//new_line('a')//&
                            '       '//temp_str(1)//new_line('a')//&
                            '       '//temp_str(2)//new_line('a')//&
                            '       '//temp_str(3)//new_line('a')
    temp_str(1)=Pointmax%x
    temp_str(2)=Pointmax%y
    temp_str(3)=Pointmax%z    
    outext = outext//'       Pointmax (BB) is '//new_line('a')//&
                            '       '//temp_str(1)//new_line('a')//&
                            '       '//temp_str(2)//new_line('a')//&
                            '       '//temp_str(3)
    
    call ToLog(outext,.false.)
    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public parameter printing routine. 
    !---------------------------------------------------------------------------
    subroutine printSimParameters
    implicit none    
    type(string) :: outext    
    type(string) :: temp_str
    call getintegratorname(temp_str,Integrator)  
    outext = '      Integrator scheme is '//temp_str//new_line('a')
    temp_str=CFL
    outext = outext//'       CFL='//temp_str//new_line('a')
    temp_str=WarmUpTime
    outext = outext//'       WarmUpTime='//temp_str//' s'//new_line('a')
    temp_str=TimeMax
    outext = outext//'       TimeMax='//temp_str//' s'//new_line('a')
    temp_str=TimeOut
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
    
    
    end module simulation_globals
