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
    use stateVector_mod
    use background_mod
    use interpolator_mod

    type :: kernel_class        !< Solver class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernel
    procedure :: run => runKernel
    procedure, private :: LagrangianKinematic
    procedure, private :: DiffusionIsotropic
    end type kernel_class

    public :: kernel_class
    contains

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Adaptation from runSolver (Ricardo) method that evaluates the specific
    !> kernel, according to the selected kernel
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function runKernel(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self

    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: runKernel

    if (sv%ttype == Globals%Types%base) then
        runKernel = self%LagrangianKinematic(sv, bdata, time, dt) + self%DiffusionIsotropic(sv, dt)
    else if (sv%ttype == Globals%Types%paper) then
        runKernel = self%DiffusionIsotropic(sv, dt)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%LagrangianKinematic(sv, bdata, time, dt)
    end if

    end function runKernel

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Lagrangian Kernel, evaluate the velocities at given points
    !> using the interpolants and split the evaluation part from the solver module.
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function LagrangianKinematic(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    integer :: np, nf, bkg, i
    real(prec) :: maxLevel(2)
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: LagrangianKinematic

    LagrangianKinematic = 0.0
    !interpolate each background
    do bkg = 1, size(bdata)
        np = size(sv%active) !number of Tracers
        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
        allocate(var_dt(np,nf))
        allocate(var_name(nf))
        !correcting for maximum admissible level in the background
        maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
        if (maxLevel(2) /= MV) where (sv%state(:,3) > maxLevel(2)) sv%state(:,3) = maxLevel(2)-0.00001
        !interpolating all of the data
        call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)
        !write dx/dt
        nf = Utils%find_str(var_name, Globals%Var%u, .true.)
        LagrangianKinematic(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)
        sv%state(:,4) = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%v, .true.)
        LagrangianKinematic(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)
        sv%state(:,5) = var_dt(:,nf)
        nf = Utils%find_str(var_name, Globals%Var%w, .false.)
        if (nf /= MV_INT) then
            LagrangianKinematic(:,3) = var_dt(:,nf)
            sv%state(:,6) = var_dt(:,nf)
        else if (nf == MV_INT) then
            LagrangianKinematic(:,3) = 0.0
            sv%state(:,6) = 0.0
        end if        
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

    end function LagrangianKinematic

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Diffusion Kernel, computes the anisotropic diffusion assuming a constant
    !> diffusion coefficient. D = 1 m/s
    !> @param[in] self, sv, dt
    !---------------------------------------------------------------------------
    function DiffusionIsotropic(self, sv, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    real(prec), intent(in) :: dt
    integer :: np, i
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v,rand_vel_w
    real(prec) :: D = 1.
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: DiffusionIsotropic

    DiffusionIsotropic = 0.0
    np = size(sv%active) !number of Tracers
    allocate(rand_vel_u(np), rand_vel_v(np), rand_vel_w(np))
    call random_number(rand_vel_u)
    call random_number(rand_vel_v)
    call random_number(rand_vel_w)
    !update velocities
    ! For the moment we set D = 1 m/s, then the D parameter
    ! should be part of the array of tracers parameter
    DiffusionIsotropic(:,1) = Utils%m2geo((2.*rand_vel_u-1.)*sqrt(2.*D/dt), sv%state(:,2), .false.)
    DiffusionIsotropic(:,2) = Utils%m2geo((2.*rand_vel_v-1.)*sqrt(2.*D/dt), sv%state(:,2), .true.)
    !DiffusionIsotropic(:,3) = (2.*rand_vel_w-1.)*sqrt(2.*D*0.0005/dt)
    where (sv%state(:,6) /= 0.0) DiffusionIsotropic(:,3) = (2.*rand_vel_w-1.)*sqrt(2.*D*0.0005/dt)

    end function DiffusionIsotropic

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