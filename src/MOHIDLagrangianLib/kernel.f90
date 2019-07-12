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
    !> Defines an abstract physics kernel class.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach. Different types of state vectors (corresponding to different
    !> types of tracers, with different quantities attached), will be affected by
    !> different processes (some suffer beaching, others don't have diffusion, etc)
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
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
    procedure, private :: StokesDrift
    procedure, private :: Windage
    procedure, private :: hasRequiredVars
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
        runKernel = self%LagrangianKinematic(sv, bdata, time, beaching = .true.) + self%StokesDrift(sv, bdata, time) + self%Windage(sv, bdata, time) !+ self%DiffusionIsotropic(sv, dt)
    else if (sv%ttype == Globals%Types%paper) then
        runKernel = self%LagrangianKinematic(sv, bdata, time, beaching = .true.) + self%StokesDrift(sv, bdata, time) + self%Windage(sv, bdata, time) !+ self%DiffusionIsotropic(sv, dt)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%LagrangianKinematic(sv, bdata, time, beaching = .true.) + self%StokesDrift(sv, bdata, time) + self%Windage(sv, bdata, time) !+ self%DiffusionIsotropic(sv, dt)
    end if

    end function runKernel

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Lagrangian Kernel, evaluate the velocities at given points
    !> using the interpolants and split the evaluation part from the solver module.
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function LagrangianKinematic(self, sv, bdata, time, beaching)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    logical, intent(in), optional :: beaching
    logical :: beach
    integer :: np, nf, bkg, i
    real(prec) :: maxLevel(2)
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: LagrangianKinematic
    real(prec), dimension(size(sv%state,1)) :: beachCoeff
    real(prec), dimension(size(sv%state,1)) :: beachCoeffRand

    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%u
    requiredVars(2) = Globals%Var%v

    beach = .false.
    if(present(beaching)) beach = beaching

    beachCoeff = 1.0
    if (beach) then
        call random_number(beachCoeffRand)
        beachCoeffRand = min(0.0, beachCoeffRand-0.15)
    end if

    LagrangianKinematic = 0.0
    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(self%hasRequiredVars(bdata(bkg)%variables, requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !correcting for maximum admissible level in the background
                maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
                if (maxLevel(2) /= MV) where (sv%state(:,3) > maxLevel(2)) sv%state(:,3) = maxLevel(2)-0.00001
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)

                !update land mask status
                nf = Utils%find_str(var_name, Globals%Var%landMask, .false.)
                if (nf /= MV_INT) sv%landMask = nint(var_dt(:,nf))
                if (nf == MV_INT) sv%landMask = Globals%Mask%waterVal
                !marking tracers for deletion because they are in land
                where(sv%landMask == 2) sv%active = .false.
                !update land interaction status
                nf = Utils%find_str(var_name, Globals%Var%landIntMask, .false.)
                if (nf /= MV_INT) then
                    sv%landIntMask = nint(var_dt(:,nf))
                    if (beach) then
                        where(sv%landIntMask == Globals%Mask%beachVal) beachCoeff = beachCoeffRand
                    end if
                end if
                if (nf == MV_INT) sv%landIntMask = Globals%Mask%waterVal

                !write dx/dt
                nf = Utils%find_str(var_name, Globals%Var%u, .true.)
                LagrangianKinematic(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)*beachCoeff
                sv%state(:,4) = var_dt(:,nf)*beachCoeff
                nf = Utils%find_str(var_name, Globals%Var%v, .true.)
                LagrangianKinematic(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)*beachCoeff
                sv%state(:,5) = var_dt(:,nf)*beachCoeff
                nf = Utils%find_str(var_name, Globals%Var%w, .false.)
                if (nf /= MV_INT) then
                    LagrangianKinematic(:,3) = var_dt(:,nf)*beachCoeff
                    sv%state(:,6) = var_dt(:,nf)*beachCoeff
                else if (nf == MV_INT) then
                    LagrangianKinematic(:,3) = 0.0
                    sv%state(:,6) = 0.0
                end if

                deallocate(var_dt)
                deallocate(var_name)
            end if
        end if
    end do

    end function LagrangianKinematic

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Computes the influence of wave velocity in tracer kinematics
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function StokesDrift(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: np, nf, bkg
    real(prec) :: waveCoeff
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: StokesDrift
    real(prec), dimension(size(sv%state,1)) :: depth

    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%vsdx
    requiredVars(2) = Globals%Var%vsdy

    waveCoeff = 0.01

    StokesDrift = 0.0
    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(self%hasRequiredVars(bdata(bkg)%variables, requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)
                !computing the depth weight
                depth = sv%state(:,3)
                where (depth>=0.0) depth = 0.0
                depth = exp(depth)
                !write dx/dt
                nf = Utils%find_str(var_name, Globals%Var%vsdx, .true.)
                StokesDrift(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)*waveCoeff*depth
                nf = Utils%find_str(var_name, Globals%Var%vsdy, .true.)
                StokesDrift(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)*waveCoeff*depth
                deallocate(var_dt)
                deallocate(var_name)
            end if
        end if
    end do

    end function StokesDrift

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Computes the influence of wind velocity in tracer kinematics
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Windage(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: np, nf, bkg
    real(prec) :: windCoeff
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Windage
    real(prec), dimension(size(sv%state,1)) :: depth

    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%u10
    requiredVars(2) = Globals%Var%v10

    windCoeff = 0.03

    Windage = 0.0
    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(self%hasRequiredVars(bdata(bkg)%variables, requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)
                !computing the depth weight
                depth = sv%state(:,3)
                where (depth>=0.0) depth = 0.0
                depth = exp(10.0*depth)
                !write dx/dt
                nf = Utils%find_str(var_name, Globals%Var%u10, .true.)
                Windage(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)*windCoeff*depth
                nf = Utils%find_str(var_name, Globals%Var%v10, .true.)
                Windage(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)*windCoeff*depth
                deallocate(var_dt)
                deallocate(var_name)
            end if
        end if
    end do

    end function Windage

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
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns true if a given set of strings (array) is contained in a list of strings
    !> @param[in] self, varList, reqVars
    !---------------------------------------------------------------------------
    logical function hasRequiredVars(self, varList, reqVars)
    class(kernel_class), intent(in) :: self
    type(stringList_class) :: varList
    type(string), dimension(:), intent(in) :: reqVars
    integer :: i
    hasRequiredVars = .true.
    do i=1, size(reqVars)
        if (varList%notRepeated(reqVars(i))) then
            hasRequiredVars = .false.
            return
        end if
    end do
    end function hasRequiredVars

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