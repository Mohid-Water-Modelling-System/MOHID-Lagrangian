module kernelAbstractions

    type :: kernel_class        !< Solver class
    integer :: kernelType = 1   !< kernel Integrator 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)
    type(string) :: name        !< Name of the Integrator algorithm
    type(interpolator_class) :: Interpolator !< The interpolator object for this Solver
contains
    procedure :: initialize => initKernel
    procedure, private :: runKernel
    procedure, private :: Lagrangian
    procedure, private :: AnisotropicDiffusion
    procedure, private :: Windage
    procedure :: print => printSolver
end type kernel_class

contains 


function Lagrangian(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class), intent(in) :: Lagrangian
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v

do bkg = 1, size(bdata)
        np = size(aot%id) !number of particles
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        call self%Interpolator%run(aot, bdata(bkg), time, var_dt, var_name)
        !update velocities
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        Lagrangian%u = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        Lagrangian%v = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%w, .true.)
        Lagrangian%w = var_dt(:,nf)
end do

end function Lagrangian


function AnisotropicDiffusion(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class), intent(in) :: Lagrangian
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v

do bkg = 1, size(bdata)
        np = size(aot%id) !number of particles
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        call random_number(var_dt)
        !update velocities
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        Diffusion%u = (2*var_dt(:,nf)-1)*sqrt(2*aot%D/dt)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        Diffusion%v = (2*var_dt(:,nf)-1)*sqrt(2*aot%D/dt)
        nf = Utils%find_str(var_name, Globals%Var%w, .true.)
        Diffusion%w = var_dt(:,nf)
end do



end function AnisotropicDiffusion

!---------------------------------------------------------------------------
!> @author Ricardo Birjukovs Canelas - MARETEC
!> @brief
!> Method that integrates the Tracer state array in one time-step, according to
!> the selected integration algorithm
!> @param[in] self, aot, bdata, time, dt
!---------------------------------------------------------------------------
subroutine runKernel(self, aot, bdata, time, dt)
    class(solver_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    if (self%kernelType == 1) call self%Lagrangian(aot, bdata, time, dt)
    if (self%kernelType == 2) call self%Diffusion(aot, bdata, time, dt)
    if (self%solverType == 3) call self%Windage(aot, bdata, time, dt)
end subroutine runKernel


!---------------------------------------------------------------------------
!> @author Ricardo Birjukovs Canelas - MARETEC
!> @brief
!> Initializer method for the Solver class. Sets the type of integrator
!> and name of the algorithm this Solver will call
!> @param[in] self, flag, name
!---------------------------------------------------------------------------
subroutine initKernel(self, flag, name)
    class(solver_class), intent(inout) :: self
    integer, intent(in) :: flag
    type(string), intent(in) :: name
    type(string) :: interpName
    self%kernelType = flag
    self%name = name
    interpName = 'linear'
    call self%Kernel%initialize(1,interpName)
end subroutine initKernel

!---------------------------------------------------------------------------
!> @author Ricardo Birjukovs Canelas - MARETEC
!> @brief
!> Method that prints the Solver information
!---------------------------------------------------------------------------
subroutine printKernel(self)
    class(solver_class), intent(inout) :: self
    type(string) :: outext, t
    outext = 'Kernel type is '//self%name
    call Log%put(outext,.false.)
    end subroutine printKernel

end module kernelAbstractions
