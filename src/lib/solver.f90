    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : solver
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an Solver class. This class invokes the different available integration
    !> algorithms as methods, and these invoke the necessary interpolation objects.
    !------------------------------------------------------------------------------

    module solver_mod

    use common_modules
    use AoT_mod
    use background_mod
    use interpolator_mod

    implicit none
    private

    type :: solver_class        !< Solver class
        integer :: solverType = 1   !< Integration Algorithm 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)
        type(string) :: name        !< Name of the Integrator algorithm
        type(interpolator_class) :: Interpolator !< The interpolator object for this Solver
    contains
    procedure :: initialize => initSolver
    procedure :: runStep
    procedure, private :: runStepMSEuler
    procedure, private :: runStepEuler
    procedure, private :: runStepRK4
    procedure :: print => printSolver
    end type solver_class

    !Public access vars
    public :: solver_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, according to
    !> the selected integration algorithm
    !> @param[in] self, aot, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStep(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec_time), intent(in) :: time, dt
    if (self%solverType == 1) call self%runStepEuler(aot, bdata, time, dt)
    if (self%solverType == 2) call self%runStepMSEuler(aot, bdata, time, dt)
    if (self%solverType == 3) call self%runStepRK4(aot, bdata, time, dt)
    end subroutine runStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Velocity Euler integration algorithm. This is a one-shot type
    !> explicit scheme with low computational cost, mostly for quick tests and
    !> debug. Implements
    !> \f$ {\vec {x}}_{t+\Delta t}={\vec {x}}_{t}+{\vec {v}}_{t}\Delta t+{\frac {1}{2}}{\vec {a}}({\vec {x}}_{t})\Delta t^{2}\f$
    !> \f$ {\vec {v}}_{t+\Delta t}={\vec {v}}_{t}+\frac{{\vec {a}}_{t+\Delta t}+{\vec {a}}_{t}}{2}\Delta t\f$
    !> @param[in] self, aot, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStepEuler(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec_time), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v 

    ! print*, 'got here!'
    ! np = size(aot%id) !number of particles
    ! allocate(rand_vel_u(np))
    ! allocate(rand_vel_v(np))
    ! call random_number(rand_vel_u)
    ! call random_number(rand_vel_v)
    ! rand_vel_u = 0.001*(2*rand_vel_u - 1)
    ! rand_vel_v = 0.001*(2*rand_vel_v - 1)
    ! aot%u = rand_vel_u
    ! aot%v = rand_vel_v
    ! !aot%w = rand_vel
    ! !print*, aot%u
    ! !update positions
    ! aot%x = aot%x + aot%u*dt
    ! aot%y = aot%y + aot%v*dt
    ! !!aot%z = aot%z + aot%w*dt

    !interpolate each background
    do bkg = 1, size(bdata)
        np = size(aot%id) !number of particles
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        !print*, 'np=', np, 'nf=', nf, 'nbkg=', size(bdata)
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        !run the interpolator
        call self%Interpolator%run(aot, bdata(bkg), time, var_dt, var_name)
        !update velocities
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        aot%u = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        aot%v = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%w, .true.)
        aot%w = var_dt(:,nf)        
        !update positions
        aot%x = aot%x + Utils%m2geo(aot%u*dt, aot%y, .false.)
        aot%y = aot%y + Utils%m2geo(aot%v*dt, aot%y, .true.)
        aot%z = aot%z + aot%w*dt
        !update other vars...
    end do

    end subroutine runStepEuler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Multi-Step Euler integration algorithm. This is a predictor-corrector type
    !> explicit scheme with excelent conservation properties and average cost
    !> @param[in] self, aot, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStepMSEuler(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec_time), intent(in) :: time, dt
    real(prec_time) :: mstime
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name

    ! interpolate each background
    do bkg = 1, size(bdata)
        np = size(aot%id) !number of particles
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        !Predictor step
        !run the interpolator
        call self%Interpolator%run(aot, bdata(bkg), time, var_dt, var_name)
        !update velocities for the predictor step
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        aot%u = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        aot%v = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%w, .true.)
        aot%w = var_dt(:,nf)
        !update positions for the predictor step
        aot%x = aot%x + Utils%m2geo(aot%u*dt, aot%y, .false.)*0.5
        aot%y = aot%y + Utils%m2geo(aot%v*dt, aot%y, .true.)*0.5
        aot%z = aot%z + aot%w*dt*0.5
        !Corrector step
        !run the interpolator
        mstime = time+0.5*dt
        call self%Interpolator%run(aot, bdata(bkg), mstime, var_dt, var_name)
        !update velocities for the corrector step
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        aot%u = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        aot%v = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%w, .true.)
        aot%w = var_dt(:,nf)
        !update positions for the corrector step
        aot%x = aot%x + Utils%m2geo(aot%u*dt, aot%y, .false.)*0.5
        aot%y = aot%y + Utils%m2geo(aot%v*dt, aot%y, .true.)*0.5
        aot%z = aot%z + aot%w*dt*0.5
        !update other vars...
    end do

    end subroutine runStepMSEuler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Runge-Kuta 4th order integration algorithm. This is an explicit scheme
    !> with medium to high computational cost
    !> @param[in] self, aot, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStepRK4(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec_time), intent(in) :: time, dt

    end subroutine runStepRK4

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializer method for the Solver class. Sets the type of integrator
    !> and name of the algorithm this Solver will call
    !> @param[in] self, flag, name
    !---------------------------------------------------------------------------
    subroutine initSolver(self, flag, name)
    class(solver_class), intent(inout) :: self
    integer, intent(in) :: flag
    type(string), intent(in) :: name
    type(string) :: interpName
    self%solverType = flag
    self%name = name
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initSolver

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the Solver information
    !---------------------------------------------------------------------------
    subroutine printSolver(self)
    class(solver_class), intent(inout) :: self
    type(string) :: outext, t
    outext = 'Solver algorithm is '//self%name
    call Log%put(outext,.false.)
    end subroutine printSolver

    end module solver_mod