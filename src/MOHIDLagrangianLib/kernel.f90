module kernel
!------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : solver
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa Paz
    !
    ! DESCRIPTION:
    !> Defines an abstraction kernel class. This class split the solver and 
    !> evaluation proces to allow other functions to be evaluated.
    !------------------------------------------------------------------------------
    use common_modules
    use AoT_mod
    use background_mod
    use interpolator_mod

    type :: kernel_class        !< Solver class
        integer :: kernelType = 1   !< kernel Integrator 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)
        type(string) :: name        !< Name of the Integrator algorithm
        type(interpolator_class) :: Interpolator !< The interpolator object for this Solver
    contains
    procedure :: initialize => initKernel
    procedure, private :: runKernel
    procedure, private :: Lagrangian
    procedure, private :: Diffusion
    procedure :: print => printKernel
end type kernel_class

contains 


!---------------------------------------------------------------------------
!> @author Daniel Garaboa Paz - USC
!> @brief
!> Lagrangian Kernel, evaluate the velocities at given points
!> using the interpolants and split the evaluation part from the solver module.
!> @param[in] self, aot, bdata, time, dt
!---------------------------------------------------------------------------
function Lagrangian(self, aot, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class) :: Lagrangian
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



function Diffusion(self, aot, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class) :: Diffusion
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v
    real(prec) :: D = 1.

do bkg = 1, size(bdata)
        np = size(aot%id) !number of particles
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        call random_number(var_dt)
        !update velocities
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        ! For the moment we set D = 1 m/s, then the D parameter
        ! should be part of the array of tracers parameter
        Diffusion%u = (2*var_dt(:,nf)-1)*sqrt(2*D/dt)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        Diffusion%v = (2*var_dt(:,nf)-1)*sqrt(2*D/dt)
        nf = Utils%find_str(var_name, Globals%Var%w, .true.)
        Diffusion%w = var_dt(:,nf)
end do



end function Diffusion

!---------------------------------------------------------------------------
!> @author Ricardo Birjukovs Canelas - MARETEC
!> @brief
!> Method that integrates the Tracer state array in one time-step, according to
!> the selected integration algorithm
!> @param[in] self, aot, bdata, time, dt
!---------------------------------------------------------------------------
function runKernel(self, aot, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(aot_class), intent(inout) :: aot
    type(aot_class) :: runKernel
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    if (self%kernelType == 1) runKernel = self%Lagrangian(aot, bdata, time, dt)
    if (self%kernelType == 2) runKernel = self%Diffusion(aot, bdata, time, dt)
end function runKernel


!---------------------------------------------------------------------------
!> @author Ricardo Birjukovs Canelas - MARETEC
!> @brief
!> Initializer method for the Solver class. Sets the type of integrator
!> and name of the algorithm this Solver will call
!> @param[in] self, flag, name
!---------------------------------------------------------------------------
subroutine initKernel(self, flag, name)
    class(kernel_class), intent(inout) :: self
    integer, intent(in) :: flag
    type(string), intent(in) :: name
    type(string) :: interpName
    self%kernelType = flag
    self%name = name
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
end subroutine initKernel

!---------------------------------------------------------------------------
!> @author Ricardo Birjukovs Canelas - MARETEC
!> @brief
!> Method that prints the Solver information
!---------------------------------------------------------------------------
subroutine printKernel(self)
    class(kernel_class), intent(inout) :: self
    type(string) :: outext, t
    outext = 'Kernel type is '//self%name
    call Log%put(outext,.false.)
    end subroutine printKernel

end module kernel