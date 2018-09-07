    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : integrator
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an integrator class
    !------------------------------------------------------------------------------

    module integrator_mod

    use common_modules
    use AoT_mod
    use background_mod

    implicit none
    private

    type :: integrator_class        !< Integrator class
        integer :: integratorType = 1   !< Integration Algorithm 1:Verlet, 2:Symplectic, 3:RK4 (default=1)
        type(string) :: name            !< Name of the Integrator algorithm
    contains
    procedure :: initialize => initIntegrator
    procedure :: runStep
    procedure, private :: runStepSymplectic
    procedure, private :: runStepVerlet
    procedure, private :: runStepRK4
    procedure :: print => printIntegrator
    end type integrator_class
    
    !Public access vars
    public :: integrator_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, according to 
    !> the selected integration algorithm
    !> @param[in] self, aot, data, dt)
    !---------------------------------------------------------------------------
    subroutine runStep(self, aot, data, time, dt)
    class(integrator_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: data
    real(prec), intent(in) :: time, dt
    if (self%integratorType == 2) call self%runStepSymplectic(aot, data, time, dt)
    if (self%integratorType == 1) call self%runStepVerlet(aot, data, time, dt)
    if (self%integratorType == 3) call self%runStepRK4(aot, data, time, dt)
    end subroutine runStep

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Symplectic integration algorithm. This is a predictor-corrector type
    !> explicit scheme with excelent conservation properties and average cost
    !> @param[in] self, aot, data, dt)
    !---------------------------------------------------------------------------
    subroutine runStepSymplectic(self, aot, data, time, dt)
    class(integrator_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: data
    real(prec), intent(in) :: time, dt
        
    end subroutine runStepSymplectic

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Verlet integration algorithm. This is a one-shot type
    !> explicit scheme with low computational cost, mostly for quick tests and 
    !> debug
    !> {\vec {x}}_{1}={\vec {x}}_{0}+{\vec {v}}_{0}\Delta t+{\frac {1}{2}}{\vec {A}}({\vec {x}}_{0})\Delta t^{2}
    !> @param[in] self, aot, data, dt)
    !---------------------------------------------------------------------------
    subroutine runStepVerlet(self, aot, data, time, dt)
    class(integrator_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: data
    real(prec), intent(in) :: time, dt
           
    end subroutine runStepVerlet

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that integrates the Tracer state array in one time-step, using a
    !> Runge-Kuta 4th order integration algorithm. This is an explicit scheme 
    !> with medium to high computational cost
    !> @param[in] self, aot, data, dt)
    !---------------------------------------------------------------------------
    subroutine runStepRK4(self, aot, data, time, dt)
    class(integrator_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: data
    real(prec), intent(in) :: time, dt
                
    end subroutine runStepRK4

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializer method for the Integrator class. Sets the type of integrator
    !> and name of the algorithm this Integrator will call
    !> @param[in] self, flag, name
    !---------------------------------------------------------------------------
    subroutine initIntegrator(self, flag, name)
    class(integrator_class), intent(inout) :: self
    integer, intent(in) :: flag
    type(string), intent(in) :: name
    self%integratorType = flag
    self%name = name
    end subroutine initIntegrator

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the Integrator information
    !---------------------------------------------------------------------------
    subroutine printIntegrator(self)
    class(integrator_class), intent(inout) :: self
    type(string) :: outext, t
    outext = 'Integrator algorithm is '//self%name
    call Log%put(outext,.false.)
    end subroutine printIntegrator

    end module integrator_mod