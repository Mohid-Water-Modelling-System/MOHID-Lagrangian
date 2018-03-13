    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_parameters
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold all the simulation related parameters, definitions and methods
    !------------------------------------------------------------------------------

    module simulation_parameters

    use commom_modules

    implicit none
    private

    !Public access vars
    public :: Integrator, CFL
    public :: InitFreeze
    public :: TimeMax, TimeOut
    public :: Dp, Pointmin, Pointmax
    public :: Gravity
    public :: Rho_ref

    !Public access procedures
    public :: setSimParameter, setSimGravity, setSimRho, setSimDp, setSimBounds

    !Parameters
    integer    :: Integrator = 1        !> Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
    real(prec) :: CFL = 0.5             !> Courant–Friedrichs–Lewy condition number
    real(prec) :: InitFreeze = 0.0      !> Time to freeze the tracers at simulation start (warmup) (default=0.0)
    real(prec) :: TimeMax = MV          !> Simulation duration
    real(prec) :: TimeOut = MV          !> Time out data (1/Hz)
    !Simulation definitions
    real(prec) :: Dp = MV               !> Initial particle spacing at source generation
    type(vector)    ::  Pointmin        !> Point that defines the lowest corner of the simulation bounding box
    type(vector)    ::  Pointmax        !> Point that defines the upper corner of the simulation bounding box
    !Case Constants
    type(vector) :: Gravity             !> Gravitational acceleration vector (default=(0 0 -9.81)) (m s-2)
    real(prec) :: Rho_ref = 1000.0      !> Reference density of the medium (default=1000.0) (kg m-3)
    

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Public parameter setting routine. Builds the simulation parametric space from the input xml case file.
    !
    !> @param[in] parmkey
    !> @param[in] parmvalue
    !---------------------------------------------------------------------------
    subroutine setSimParameter(parmkey,parmvalue)
    implicit none
    type(string), intent(in) :: parmkey
    type(string), intent(in) :: parmvalue
    !add new parameters to this search three
    if (parmkey%chars()=="Integrator") then
        Integrator=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="CFL") then
        CFL=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="InitFreeze") then
        InitFreeze=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="TimeMax") then
        TimeMax=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="TimeOut") then
        TimeOut=parmvalue%to_number(prec)
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
    !> @param[in] rho
    !---------------------------------------------------------------------------
    subroutine setSimRho(read_rho)
    implicit none
    type(string), intent(in) :: read_rho    
    character(80) :: chars    
    chars = read_rho%chars()
    read(chars,*)Rho_ref    
    if (Rho_ref.le.0.0) then
        stop 'Rho_ref must be positive and non-zero, stopping.'
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
    character(80) :: chars    
    chars = read_dp%chars()
    read(chars,*)Dp    
    if (Dp.le.0.0) then
        stop 'Dp must be positive and non-zero, stopping.'
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
    
    

    end module simulation_parameters
