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

    use kernelLitter_mod
    use kernelVerticalMotion_mod

    type :: kernel_class        !< Kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernel
    procedure :: run => runKernel
    procedure, private :: setCommonProcesses
    procedure, private :: LagrangianKinematic
    procedure, private :: DiffusionMixingLength
    procedure, private :: DiffusionIsotropic
    procedure, private :: StokesDrift
    procedure, private :: Windage
    procedure, private :: Beaching
    procedure, private :: Aging
    end type kernel_class

    type(kernelLitter_class) :: Litter       !< litter kernels
    type(kernelVerticalMotion_class) :: VerticalMotion   !< VerticalMotion kernels
 
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

    !running preparations for kernel lanch
    call self%setCommonProcesses(sv, bdata, time)

    !running kernels for each type of tracer
    if (sv%ttype == Globals%Types%base) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv)
        runKernel = self%Beaching(sv, runKernel)
    else if (sv%ttype == Globals%Types%paper) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + Litter%DegradationLinear(sv) + VerticalMotion%Buoyancy(sv, bdata, time) + &
                    VerticalMotion%Resuspension(sv, bdata, time)
        runKernel = self%Beaching(sv, runKernel)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + Litter%DegradationLinear(sv) + VerticalMotion%Buoyancy(sv, bdata, time) + &
                    VerticalMotion%Resuspension(sv, bdata, time)
        runKernel = self%Beaching(sv, runKernel)
    end if
    
    runKernel = VerticalMotion%CorrectVerticalBounds(sv, runKernel, bdata, time, dt)

    end function runKernel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Sets the state vector land interaction mask values and corrects for
    !> maximum level of tracers.
    !> Accounts for global periodicity.
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    subroutine setCommonProcesses(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: np, nf, bkg
    real(prec) :: maxLevel(2)
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars

    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%landIntMask
    requiredVars(2) = Globals%Var%resolution

    ! global periodicity conditions
    where (sv%state(:,1) > 180.0) sv%state(:,1) = sv%state(:,1) - 360.0
    where (sv%state(:,1) < -180.0) sv%state(:,1) = sv%state(:,1) + 360.0

    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(bdata(bkg)%hasVars(requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !correcting for maximum admissible level in the background
                maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
                if (maxLevel(2) /= MV) where (sv%state(:,3) > maxLevel(2)) sv%state(:,3) = maxLevel(2)-0.00001
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                !update land interaction status
                nf = Utils%find_str(var_name, Globals%Var%landIntMask)
                sv%landIntMask = var_dt(:,nf)
                !marking tracers for deletion because they are in land
                if (Globals%SimDefs%RemoveLandTracer == 1) then
                     where(int(sv%landIntMask + Globals%Mask%landVal*0.05) == Globals%Mask%landVal) sv%active = .false.
                end if
                !update resolution proxy
                nf = Utils%find_str(var_name, Globals%Var%resolution,.true.)
                sv%resolution = var_dt(:,nf)
                deallocate(var_name)
                deallocate(var_dt)
            end if
        end if
    end do

    end subroutine setCommonProcesses

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Lagrangian Kernel, evaluate the velocities at given points
    !> using the interpolants and split the evaluation part from the solver module.
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function LagrangianKinematic(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: np, nf, bkg, i
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1), size(sv%state,2)) :: LagrangianKinematic

    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%u
    requiredVars(2) = Globals%Var%v

    LagrangianKinematic = 0.0
    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(bdata(bkg)%hasVars(requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)
                !write dx/dt
                nf = Utils%find_str(var_name, Globals%Var%u, .true.)
                LagrangianKinematic(:,1) = Utils%m2geo(var_dt(:, nf), sv%state(:,2), .false.)
                sv%state(:,4) = var_dt(:,nf)
                nf = Utils%find_str(var_name, Globals%Var%v, .true.)
                LagrangianKinematic(:,2) = Utils%m2geo(var_dt(:, nf), sv%state(:,2), .true.)
                sv%state(:,5) = var_dt(:,nf)
                nf = Utils%find_str(var_name, Globals%Var%w, .false.)
                if ((nf /= MV_INT) .and. (Globals%SimDefs%VerticalVelMethod == 1)) then
                    LagrangianKinematic(:,3) = var_dt(:, nf)
                    sv%state(:,6) = var_dt(:, nf)
                else if ((nf /= MV_INT) .and. (Globals%SimDefs%VerticalVelMethod == 2)) then
                    LagrangianKinematic(:,3) = VerticalMotion%Divergence(sv, bdata, time)
                    sv%state(:,6) = LagrangianKinematic(:,3)
                else if ((nf == MV_INT) .or. (Globals%SimDefs%VerticalVelMethod == 3)) then
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
            if(bdata(bkg)%hasVars(requiredVars)) then
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
                where(sv%landIntMask < Globals%Mask%landVal) StokesDrift(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)*waveCoeff*depth
                nf = Utils%find_str(var_name, Globals%Var%vsdy, .true.)
                where(sv%landIntMask < Globals%Mask%landVal) StokesDrift(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)*waveCoeff*depth
                deallocate(var_name)
                deallocate(var_dt)
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
            if(bdata(bkg)%hasVars(requiredVars)) then
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
                where(sv%landIntMask < Globals%Mask%landVal) Windage(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)*windCoeff*depth
                nf = Utils%find_str(var_name, Globals%Var%v10, .true.)
                where(sv%landIntMask < Globals%Mask%landVal) Windage(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)*windCoeff*depth
                deallocate(var_name)
                deallocate(var_dt)
            end if
        end if
    end do

    end function Windage

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Beaching Kernel, uses the already updated state vector and determines if
    !> and how beaching occurs. Affects the state vector and state vector derivative.
    !> @param[in] self, sv, svDt
    !---------------------------------------------------------------------------
    function Beaching(self, sv, svDt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)), intent(in) :: svDt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Beaching
    real(prec), dimension(size(sv%state,1)) :: beachCoeff
    real(prec), dimension(size(sv%state,1)) :: beachCoeffRand
    real(prec), dimension(size(sv%state,1)) :: beachWeight
    real(prec) :: lbound, ubound
    integer :: i

    beachCoeff = 1.0
    call random_number(beachCoeffRand) !this is a uniform distribution generator
    beachCoeffRand = max(0.0, beachCoeffRand - Globals%Constants%BeachingStopProb)  !clipping the last % to zero
    where (beachCoeffRand /= 0.0 ) beachCoeffRand = beachCoeffRand*(1.0/maxval(beachCoeffRand))

    Beaching = svDt

    if (Globals%Constants%BeachingStopProb /= 0.0) then !beaching is completely turned off if the stopping propability is zero

        !getting the bounds for the interpolation of the land interaction field that correspond to beaching
        lbound = (Globals%Mask%beachVal + Globals%Mask%waterVal)*0.5
        ubound = (Globals%Mask%beachVal + Globals%Mask%landVal)*0.5

        !beachWeight = 1 - 0.5*(sv%landIntMask - lbound)/(ubound-lbound) !linear distance weight for beaching
        beachWeight = 1 - 0.9*(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound) !quadratic weight

        !replacing 1.0 with a coefficient from beaching where needed
        where(sv%landIntMask <= ubound .and. sv%landIntMask >= lbound) beachCoeff = beachCoeffRand*beachWeight
        do i=1,3
            Beaching(:,i) = svDt(:,i)*beachCoeff !position derivative is affected
            sv%state(:,i+3) = sv%state(:,i+3)*beachCoeff !so are the velocities
        end do
        !zero vertical velocity in beached particles?
        !where (beachCoeff /= 1.0)
        !    Beaching(:,3) = 0.0
        !    sv%state(:,6) = 0.0
        !end where

    end if

    end function Beaching

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Aging kernel. Sets the age variable to be updated by dt by the solver
    !> @param[in] self, sv
    !---------------------------------------------------------------------------
    function Aging(self, sv)
    class(kernel_class), intent(in) :: self
    type(stateVector_class), intent(in) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Aging
    integer :: nf
    type(string) :: tag

    Aging = 0.0
    tag = 'age'
    nf = Utils%find_str(sv%varName, tag, .true.)
    !setting the age variable to be updated by dt by the solver for all tracers
    !effectively, we're just setting the derivative of 'age' to be 1.0 :)
    Aging(:,nf) = 1.0

    end function Aging

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> mixing length diffusion kernel, computes random velocities at given
    !> instants to model diffusion processes. These are valid while the tracer
    !> travels a given mixing length, propotional to the resolution of the
    !> background (and its ability to resove motion scales)
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function DiffusionMixingLength(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), intent(in) :: dt
    integer :: np, nf, bkg
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: DiffusionMixingLength
    real(prec), dimension(:), allocatable :: resolution
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v, rand_vel_w

    allocate(requiredVars(1))
    requiredVars(1) = Globals%Var%resolution

    DiffusionMixingLength = 0.0

    if (Globals%Constants%DiffusionCoeff == 0.0) return

    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(bdata(bkg)%hasVars(requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_name(nf))
                allocate(var_dt(np,nf))
                allocate(resolution(np))
                allocate(rand_vel_u(np), rand_vel_v(np), rand_vel_w(np))
                call random_number(rand_vel_u)
                call random_number(rand_vel_v)
                call random_number(rand_vel_w)
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                !update resolution estimate
                nf = Utils%find_str(var_name, Globals%Var%resolution, .true.)
                resolution = var_dt(:,nf)
                !if we are still in the same path, use the same random velocity, do nothing
                !if we ran the path, new random velocities are generated and placed
                where ((sv%state(:,10) > 2.0*resolution) .and. (sv%landIntMask < Globals%Mask%landVal))
                    DiffusionMixingLength(:,7) = (2.*rand_vel_u-1.)*sqrt(Globals%Constants%DiffusionCoeff*abs(sv%state(:,4))/dt)/dt
                    DiffusionMixingLength(:,8) = (2.*rand_vel_v-1.)*sqrt(Globals%Constants%DiffusionCoeff*abs(sv%state(:,5))/dt)/dt
                    DiffusionMixingLength(:,9) = (2.*rand_vel_w-1.)*sqrt(0.000001*Globals%Constants%DiffusionCoeff*abs(sv%state(:,6))/dt)/dt
                    sv%state(:,10) = 0.0
                    !update system positions
                    DiffusionMixingLength(:,1) = Utils%m2geo(DiffusionMixingLength(:,7), sv%state(:,2), .false.)*dt
                    DiffusionMixingLength(:,2) = Utils%m2geo(DiffusionMixingLength(:,8), sv%state(:,2), .true.)*dt
                    DiffusionMixingLength(:,3) = DiffusionMixingLength(:,9)*dt
                elsewhere
                    !update system positions
                    DiffusionMixingLength(:,1) = Utils%m2geo(sv%state(:,7), sv%state(:,2), .false.)
                    DiffusionMixingLength(:,2) = Utils%m2geo(sv%state(:,8), sv%state(:,2), .true.)
                    DiffusionMixingLength(:,3) = sv%state(:,9)
                end where
                !update system velocities
                !sv%state(:,4) = sv%state(:,4) + DiffusionMixingLength(:,7)*dt
                !sv%state(:,5) = sv%state(:,5) + DiffusionMixingLength(:,8)*dt
                !sv%state(:,6) = sv%state(:,6) + DiffusionMixingLength(:,9)*dt
                !update used mixing length
                DiffusionMixingLength(:,10) = sqrt(sv%state(:,4)*sv%state(:,4) + sv%state(:,5)*sv%state(:,5) + sv%state(:,6)*sv%state(:,6))
                deallocate(rand_vel_u, rand_vel_v, rand_vel_w)
                deallocate(resolution)
                deallocate(var_dt)
                deallocate(var_name)
            end if
        end if
    end do
    end function DiffusionMixingLength

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
    DiffusionIsotropic(:,1) = Utils%m2geo((2.*rand_vel_u-1.)*sqrt(2.*Globals%Constants%DiffusionCoeff/dt), sv%state(:,2), .false.)
    DiffusionIsotropic(:,2) = Utils%m2geo((2.*rand_vel_v-1.)*sqrt(2.*Globals%Constants%DiffusionCoeff/dt), sv%state(:,2), .true.)
    !DiffusionIsotropic(:,3) = (2.*rand_vel_w-1.)*sqrt(2.*D*0.0005/dt)
    where (sv%state(:,6) /= 0.0) DiffusionIsotropic(:,3) = (2.*rand_vel_w-1.)*sqrt(2.*Globals%Constants%DiffusionCoeff*0.0005/dt)

    end function DiffusionIsotropic

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernel(self)
    class(kernel_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    call Litter%initialize()
    call VerticalMotion%initialize()
    end subroutine initKernel

    end module kernel_mod