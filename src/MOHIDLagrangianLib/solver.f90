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
    !> algorithms as methods, and these invoke the necessary kernel functions.
    !------------------------------------------------------------------------------

    module solver_mod


    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod
    use kernel_mod

    implicit none
    private

    type :: solver_class        !< Solver class
        integer :: solverType = 1   !< Integration Algorithm 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)
        type(string) :: name        !< Name of the Integrator algorithm
        type(kernel_class) :: Kernel
    contains
    procedure :: initialize => initSolver
    procedure :: runStep
    procedure, private :: runStepEuler
    procedure, private :: runStepMSEuler
    procedure, private :: runStepRK4
    procedure :: print => printSolver
    end type solver_class

    !Public access vars
    public :: solver_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer State Vector in one time-step, according to
    !> the selected integration algorithm
    !> @param[in] self, state, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStep(self, state, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(stateVector_class), dimension(:), intent(inout) :: state
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    !so the forward integrators don't overextend beyond calendar time
    if (time+dt < Globals%Parameters%TimeMax) then
        if (self%solverType == 1) call self%runStepEuler(state, bdata, time, dt)
        if (self%solverType == 2) call self%runStepMSEuler(state, bdata, time, dt)
        if (self%solverType == 3) call self%runStepRK4(state, bdata, time, dt)
    end if
    end subroutine runStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Velocity Euler integration algorithm. This is a one-shot type
    !> explicit scheme with low computational cost. Implements
    !> \f$ {\vec {x}}_{t+\Delta t}={\vec {x}}_{t}+{\vec {v}}_{t}\Delta t+{\frac {1}{2}}{\vec {a}}({\vec {x}}_{t})\Delta t^{2}\f$
    !> \f$ {\vec {v}}_{t+\Delta t}={\vec {v}}_{t}+\frac{{\vec {a}}_{t+\Delta t}+{\vec {a}}_{t}}{2}\Delta t\f$
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStepEuler(self, sv, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(stateVector_class), dimension(:), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: i

    do i=1, size(sv)
        sv(i)%state = sv(i)%state + self%Kernel%run(sv(i), bdata, time, dt)*dt
    end do

    end subroutine runStepEuler

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state vector in one time-step, using a
    !> Multi-Step Euler integration algorithm. This is a predictor-corrector type
    !> explicit scheme with excelent conservation properties and average cost
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStepMSEuler(self, sv, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(stateVector_class), dimension(:), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    type(stateVector_class) :: predSv
    real(prec), allocatable, dimension(:,:) :: predKernel
    real(prec):: mstime
    integer :: i

    mstime = time + 0.5*dt

    do i=1, size(sv)
        !creating predictor step state vector
        call sv(i)%copyState(predSv)
        allocate(predKernel(size(sv(i)%state, 1), size(sv(i)%state, 2)))
        !computing predictor step
        predKernel = self%Kernel%run(sv(i), bdata, time, dt)
        predSv%state = sv(i)%state + predKernel*dt
        !computing corrector step
        sv(i)%state = sv(i)%state + (predKernel +  self%Kernel%run(predSv, bdata, mstime, dt))*(dt*0.5)
        !deallocating
        call predSv%finalize()
        deallocate(predKernel)
    end do

    end subroutine runStepMSEuler


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Runge-Kuta 4th order integration algorithm. This is an explicit scheme
    !> with medium to high computational cost
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    subroutine runStepRK4(self, sv, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(stateVector_class), dimension(:), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    type(stateVector_class) :: intSv
    real(prec), allocatable, dimension(:,:) :: k1, k2, k3, k4
    integer :: i

    do i=1, size(sv)
        !creating intermediate step state vectors
        call sv(i)%copyState(intSv)
        allocate(k1(size(sv(i)%state, 1), size(sv(i)%state, 2)))
        allocate(k2(size(sv(i)%state, 1), size(sv(i)%state, 2)))
        allocate(k3(size(sv(i)%state, 1), size(sv(i)%state, 2)))
        allocate(k4(size(sv(i)%state, 1), size(sv(i)%state, 2)))
        !1st step
        k1 = self%Kernel%run(sv(i), bdata, time, dt)*dt
        intSv%state = sv(i)%state + k1*0.5
        !2nd step
        k2 = self%Kernel%run(intSv, bdata, time + 0.5*dt, dt)*dt
        intSv%state = sv(i)%state + k2*0.5
        !3rd step
        k3 = self%Kernel%run(intSv, bdata, time + 0.5*dt, dt)*dt
        intSv%state = sv(i)%state + k3
        !4th step
        k4 = self%Kernel%run(intSv, bdata, time + dt, dt)*dt
        !computing the new state
        sv(i)%state = sv(i)%state + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        !deallocating
        call intSv%finalize()
        deallocate(k4)
        deallocate(k3)
        deallocate(k2)
        deallocate(k1)
        
    end do

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
    call self%Kernel%initialize()
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