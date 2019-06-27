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
    use stateVector_mod
    use background_mod
    use interpolator_mod

    type :: kernel_class        !< Solver class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernel
    procedure :: run => runKernel
    procedure, private :: Lagrangian
    !procedure, private :: Diffusion
    end type kernel_class

    public :: kernel_class
    contains

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Adaptation from runSolver (Ricardo) method that evaluates the specific
    !> kernel, according to the selected kernel
    !> @param[in] self, aot, bdata, time, dt
    !> @param[out] daot_dt
    !---------------------------------------------------------------------------
    function runKernel(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self

    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: runKernel

    if (sv%ttype == Globals%Types%base) then
        runKernel = self%Lagrangian(sv, bdata, time, dt)! + self%Diffusion(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%paper) then
        runKernel = self%Lagrangian(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%Lagrangian(sv, bdata, time, dt)
    end if


    end function runKernel

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Lagrangian Kernel, evaluate the velocities at given points
    !> using the interpolants and split the evaluation part from the solver module.
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function Lagrangian(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg, i
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Lagrangian

    Lagrangian = 0.0    
    !interpolate each background
    do bkg = 1, size(bdata)
        np = size(sv%active) !number of particles
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)
        !write dx/dt
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        Lagrangian(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        Lagrangian(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)
        nf = Utils%find_str(var_name, Globals%Var%w, .false.)
        if (nf /= MV_INT) Lagrangian(:,3) = var_dt(:,nf)
        if (nf == MV_INT) Lagrangian(:,3) = 0.0
        !write vel
        do i=1,3
            sv%state(:,i+3) = Lagrangian(:,i)
        end do
        !update land mask status
        nf = Utils%find_str(var_name, Globals%Var%landMask, .false.)
        if (nf /= MV_INT) sv%landMask = nint(var_dt(:,nf))
        if (nf == MV_INT) sv%landMask = Globals%Mask%waterVal
        !marking tracers for deletion because they are in land
        where(sv%landMask == 2) sv%active = .false.
        !update land interaction status
        nf = Utils%find_str(var_name, Globals%Var%landIntMask, .false.)
        if (nf /= MV_INT) sv%landIntMask = nint(var_dt(:,nf))
        if (nf == MV_INT) sv%landIntMask = Globals%Mask%waterVal
        !update other vars...
    end do

    end function Lagrangian


    !!---------------------------------------------------------------------------
    !!> @author Daniel Garaboa Paz - USC
    !!> @brief
    !!> Diffusion Kernel, computes the anisotropic diffusion assuming a constant
    !!> diffusion coefficient. D = 1 m/s
    !!> @param[in] self, aot, bdata, time, dt
    !!> @param[out] daot_dt
    !!---------------------------------------------------------------------------
    !subroutine Diffusion(self, aot, bdata, time, dt, daot_dt)
    !class(kernel_class), intent(inout) :: self
    !type(aot_class), intent(in) :: aot
    !type(aot_class),intent(out) :: daot_dt
    !type(background_class), dimension(:), intent(in) :: bdata
    !real(prec), intent(in) :: time, dt
    !integer :: np, nf, bkg
    !real(prec), dimension(:,:), allocatable :: var_dt
    !type(string), dimension(:), allocatable :: var_name
    !real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v,rand_vel_w
    !real(prec) :: D = 1.
    !
    !np = size(aot%id)
    !allocate(rand_vel_u(np), rand_vel_v(np), rand_vel_w(np))
    !call random_number(rand_vel_u)
    !call random_number(rand_vel_v)
    !call random_number(rand_vel_w)
    !!update velocities
    !! For the moment we set D = 1 m/s, then the D parameter
    !! should be part of the array of tracers parameter
    !daot_dt%u = (2.*rand_vel_u-1.)*sqrt(2.*D/dt)
    !daot_dt%v = (2.*rand_vel_v-1.)*sqrt(2.*D/dt)
    !daot_dt%w = (2.*rand_vel_w-1.)*sqrt(2.*D/dt)
    !
    !
    !end subroutine Diffusion




    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for solver kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernel(self)
    class(kernel_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initKernel

    end module kernel_mod