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
        integer :: solverType = 1   !< Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
        type(string) :: name        !< Name of the Integrator algorithm
        type(interpolator_class) :: Interpolator !< The interpolator object for this Solver
    contains
    procedure :: initialize => initSolver
    procedure :: runStep
    procedure, private :: runStepSymplectic
    procedure, private :: runStepVerlet
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
    !> @param[in] self, aot, bdata, dt
    !---------------------------------------------------------------------------
    subroutine runStep(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    if (self%solverType == 2) call self%runStepSymplectic(aot, bdata, time, dt)
    if (self%solverType == 1) call self%runStepVerlet(aot, bdata, time, dt)
    if (self%solverType == 3) call self%runStepRK4(aot, bdata, time, dt)
    end subroutine runStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Symplectic integration algorithm. This is a predictor-corrector type
    !> explicit scheme with excelent conservation properties and average cost
    !> @param[in] self, aot, bdata, dt
    !---------------------------------------------------------------------------
    subroutine runStepSymplectic(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
        
    end subroutine runStepSymplectic

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Verlet integration algorithm. This is a one-shot type
    !> explicit scheme with low computational cost, mostly for quick tests and 
    !> debug. Implements 
    !> \f$ {\vec {x}}_{1}={\vec {x}}_{0}+{\vec {v}}_{0}\Delta t+{\frac {1}{2}}{\vec {A}}({\vec {x}}_{0})\Delta t^{2}\f$
    !> @param[in] self, aot, bdata, dt
    !---------------------------------------------------------------------------
    subroutine runStepVerlet(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt

    ! interpolate each background 
    do bkg = 1, size(bdata)
        np = size(aot%id)
        nf = bdata(bkg)%fields%getSize()
        allocate(var_dt(np,nf))
        call self%Interpolator%run(aot, bdata(bkg), time, dt, var_dt)
        !put the interpolated vars from var_dt back into the AoT
    end do 

    !now that we interpolated the variables, we need to know what to do with them. 
    !need to find a way to save velocity to velocity and temperature to temperature.
    !list with vars to interpolate from input xml, and config xml with standard
    !netcdf name overrides is needed for specific institution files?
           
    end subroutine runStepVerlet

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Runge-Kuta 4th order integration algorithm. This is an explicit scheme 
    !> with medium to high computational cost
    !> @param[in] self, aot, bdata, dt
    !---------------------------------------------------------------------------
    subroutine runStepRK4(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
                
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
    self%solverType = flag
    self%name = name
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