module kernel_mod
!------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernel
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa Paz
    !
    ! DESCRIPTION:
    !> Defines an abstraction kernel class adpated from solver class. 
    !> This class split the solver and evaluation proces to allow other functions 
    !> to be evaluated.
    !------------------------------------------------------------------------------
    use common_modules
    use AoT_mod
    use background_mod
    use interpolator_mod

    type :: kernel_class        !< Solver class
        integer :: kernelType = 1   !< kernel Integrator 1:Lagrangian, 2:Diffusion
        type(string) :: name        !< Name of the kernel
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernel
    procedure :: runKernel
    procedure, private :: Lagrangian
    procedure, private :: Diffusion
    procedure :: print => printKernel
end type kernel_class

public :: kernel_class
contains 


!---------------------------------------------------------------------------
!> @author Daniel Garaboa Paz - USC
!> @brief
!> Lagrangian Kernel, evaluate the velocities at given points
!> using the interpolants and split the evaluation part from the solver module.
!> @param[in] self, aot, bdata, time, dt
!> @param[out] daot_dt
!---------------------------------------------------------------------------
subroutine Lagrangian(self, aot, bdata, time, dt, daot_dt)
    class(kernel_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class), intent(out) :: daot_dt
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name

    daot_dt = aot
    
    do bkg = 1, size(bdata)
            np = size(aot%id) !number of particles
            nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
            allocate(var_dt(np,nf))
            allocate(var_name(nf))
            call self%Interpolator%run(aot, bdata(bkg), time, var_dt, var_name)
            !update velocities
            nf = Utils%find_str(var_name, Globals%Var%u, .true.)
            daot_dt%u = var_dt(:,nf)
            nf = Utils%find_str(var_name, Globals%Var%v, .true.)
            daot_dt%v = var_dt(:,nf)
            nf = Utils%find_str(var_name, Globals%Var%w, .false.)
            if (nf /= MV_INT) daot_dt%w = var_dt(:,nf)
            if (nf == MV_INT) daot_dt%w = 0.0
            daot_dt%u = Utils%m2geo(daot_dt%u, aot%y, .False.)
            daot_dt%v = Utils%m2geo(daot_dt%v, aot%y, .True.)
    end do

end subroutine Lagrangian


!---------------------------------------------------------------------------
!> @author Daniel Garaboa Paz - USC
!> @brief
!> Diffusion Kernel, computes the anisotropic diffusion assuming a constant
!> diffusion coefficient. D = 1 m/s
!> @param[in] self, aot, bdata, time, dt
!> @param[out] daot_dt
!---------------------------------------------------------------------------
subroutine Diffusion(self, aot, bdata, time, dt, daot_dt)
    class(kernel_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class),intent(out) :: daot_dt
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v,rand_vel_w
    real(prec) :: D = 0.1

    daot_dt = aot

    ! Advection term
    do bkg = 1, size(bdata)
            np = size(aot%id) !number of particles
            nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
            allocate(var_dt(np,nf))
            allocate(var_name(nf))
            call self%Interpolator%run(aot, bdata(bkg), time, var_dt, var_name)
            !update velocities
            nf = Utils%find_str(var_name, Globals%Var%u, .true.)
            daot_dt%u = var_dt(:,nf)
            nf = Utils%find_str(var_name, Globals%Var%v, .true.)
            daot_dt%v = var_dt(:,nf)
            nf = Utils%find_str(var_name, Globals%Var%w, .false.)
            if (nf /= MV_INT) daot_dt%w = var_dt(:,nf)
            if (nf == MV_INT) daot_dt%w = 0.0

    end do

    ! Diffusion term
    allocate(rand_vel_u(np), rand_vel_v(np), rand_vel_w(np))
    call random_number(rand_vel_u)
    call random_number(rand_vel_v)
    call random_number(rand_vel_w)
    !update velocities
    ! For the moment we set D = 0.1 m/s, then the D parameter
    ! should be part of the array of tracers parameter

    ! Advection + Diffusion! (we neglect the w diffusion component)
    ! bounds error from bottom-top layer!
    daot_dt%u = daot_dt%u + (2.*rand_vel_u-1.)*sqrt(2.*D/dt)
    daot_dt%v = daot_dt%v + (2.*rand_vel_v-1.)*sqrt(2.*D/dt)
    !daot_dt%w = daot_dt%w + (2.*rand_vel_w-1.)*sqrt(2.*0.001/dt)


    daot_dt%u = Utils%m2geo(daot_dt%u, aot%y, .False.)
    daot_dt%v = Utils%m2geo(daot_dt%v, aot%y, .True.)


end subroutine Diffusion


!---------------------------------------------------------------------------
!> @author Daniel Garaboa Paz - GFNL
!> @brief
!> Adaptation from runSolver (Ricardo) method that evaluates the specific 
!> kernel, according to the selected kernel
!> @param[in] self, aot, bdata, time, dt
!> @param[out] daot_dt
!---------------------------------------------------------------------------
subroutine runKernel(self, aot, bdata, time, dt,daot_dt)
    class(kernel_class), intent(inout) :: self
    type(aot_class), intent(in) :: aot
    type(aot_class) :: daot_dt
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    if (self%kernelType == 1) call self%Lagrangian(aot, bdata, time, dt, daot_dt)
    if (self%kernelType == 2) call self%Diffusion(aot, bdata, time, dt, daot_dt)
end subroutine runKernel


!---------------------------------------------------------------------------
!> @author Daniel Garaboa Paz - GFNL
!> @brief
!> Initializer method adpated from for solver kernel class. Sets the type of 
!> kernel and the interpolator to evaluate it.
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
!> Method that prints the Kernel information
!---------------------------------------------------------------------------
subroutine printKernel(self)
    class(kernel_class), intent(inout) :: self
    type(string) :: outext, t
    outext = 'Kernel type is '//self%name
    call Log%put(outext,.false.)
    end subroutine printKernel

end module kernel_mod