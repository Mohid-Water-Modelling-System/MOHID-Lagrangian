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
    !> Module to hold all the simulation related parameters and methods
    !------------------------------------------------------------------------------

    module simulation_parameters

    use commom_modules

    implicit none
    private

    !Public access vars
    public :: Integrator
    public :: InitFreeze
    public :: TimeMax, TimeOut

    !Public access procedures
    public :: setSimParameter

    integer    :: Integrator = 1        !> Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
    real(prec) :: InitFreeze = 0.0      !> Time to freeze the tracers at simulation start (warmup) (default=0.0)
    real(prec) :: TimeMax = MV          !> Simulation duration
    real(prec) :: TimeOut = MV          !> Time out data (1/Hz)

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

    !need to add new parameters to this search three
    if (parmkey%chars()=="Integrator") then
        Integrator=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="InitFreeze") then
        InitFreeze=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="TimeMax") then
        TimeMax=parmvalue%to_number(prec)
    elseif(parmkey%chars()=="TimeOut") then
        TimeOut=parmvalue%to_number(prec)
    endif

    end subroutine

    end module simulation_parameters
