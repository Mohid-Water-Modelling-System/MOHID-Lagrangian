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
    use kernelUtils_mod
    
    use kernelLitter_mod
    use kernelVerticalMotion_mod
    use kernelColiform_mod
    use kernelDetritus_mod

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
    type(kernelColiform_class) :: Coliform !< coliform kernels
    type(kernelDetritus_class) :: Detritus !< coliform kernels
    type(kernelUtils_class) :: KernelUtils   !< kernel utils
 
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
                    VerticalMotion%Resuspension(sv, bdata, time, dt)
        runKernel = self%Beaching(sv, runKernel)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + Litter%DegradationLinear(sv) + VerticalMotion%Buoyancy(sv, bdata, time) + &
                    VerticalMotion%Resuspension(sv, bdata, time, dt)
        runKernel = self%Beaching(sv, runKernel)
    else if (sv%ttype == Globals%Types%coliform) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + coliform%MortalityT90(sv, bdata, time) + coliform%Dilution(sv, bdata, time, dt)
        runKernel = self%Beaching(sv, runKernel)
    else if (sv%ttype == Globals%Types%seed) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%DiffusionMixingLength(sv, bdata, time, dt) + self%Aging(sv) + &
                    VerticalMotion%Buoyancy(sv, bdata, time) + VerticalMotion%Resuspension(sv, bdata, time, dt)
        runKernel = self%Beaching(sv, runKernel)
    else if (sv%ttype == Globals%Types%detritus) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%DiffusionMixingLength(sv, bdata, time, dt) + self%Aging(sv) + &
                    VerticalMotion%Buoyancy(sv, bdata, time) + VerticalMotion%Resuspension(sv, bdata, time, dt) + &
                    detritus%Degradation(sv, bdata, time, dt)
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
    integer :: i, j, col_age, col_bat, col_bat_sv, col_dwz, col_dwz_sv, col_landintmask, col_res
    real(prec) :: maxLevel(2)
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    type(string) :: tag
    logical bottom_emmission
    !-----------------------------------------------------------
    allocate(requiredVars(4))
    requiredVars(1) = Globals%Var%landIntMask
    requiredVars(2) = Globals%Var%resolution
    requiredVars(3) = Globals%Var%bathymetry
    requiredVars(4) = Globals%Var%dwz
    
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
    bottom_emmission = .false.
    col_bat = Utils%find_str(var_name, Globals%Var%bathymetry, .true.)
    !Set tracers bathymetry
    col_bat_sv = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
    sv%state(:,col_bat_sv) = var_dt(:,col_bat)
    !Set tracers dwz
    col_dwz = Utils%find_str(var_name, Globals%Var%dwz, .true.)
    col_dwz_sv = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    sv%state(:,col_dwz_sv) = var_dt(:,col_dwz)
    
    tag = 'age'
    col_age = Utils%find_str(sv%varName, tag, .true.)
    
    where (sv%state(:,3) < var_dt(:,col_bat)) sv%state(:,3) = var_dt(:,col_bat)
        
    if (size(sv%source) > 0) then
        !if any of the sources defined by the user has the option bottom_emission then the model must check
        !wheter any new tracer needs to be positioned at the bottom
        if (maxval(Globals%Sources%bottom_emission_depth) > 0) then
            bottom_emmission = .true.
        end if
    end if
    ! global periodicity conditions
    where (sv%state(:,1) > 180.0) sv%state(:,1) = sv%state(:,1) - 360.0
    where (sv%state(:,1) < -180.0) sv%state(:,1) = sv%state(:,1) + 360.0
    !interpolate each background
    
    !correcting for maximum admissible level in the background
    maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .false.)   
    if (maxLevel(2) /= MV) where (sv%state(:,3) > maxLevel(2)) sv%state(:,3) = maxLevel(2)-0.00001
    !update land interaction status
    col_landintmask = Utils%find_str(var_name, Globals%Var%landIntMask)
    sv%landIntMask = var_dt(:,col_landintmask)
    !if bottom emission is active, check if tracer age is 0 (has just been added to the simulation)
    ! and if true, place those particles at the bottom (bathymetric value)

    if (bottom_emmission) then
        if ((size(sv%source) > 0) .and. (minval(sv%state(:,col_age)) == 0)) then
            !where the age of a tracer is 0, make the vertical position equal to the bathymetric value of the grid cell
            !where the tracer is located
            where (sv%state(:,col_age) == 0 .and. Globals%Sources%bottom_emission_depth(sv%source(:)) > 0) sv%state(:,3) = var_dt(:,col_bat) + Globals%Sources%bottom_emission_depth(sv%source(:))
        end if
    end if

    !marking tracers for deletion because they are in land
    if (Globals%simdefs%removelandtracer == 1) then
        where(int(sv%landintmask + Globals%mask%landval*0.05) == Globals%mask%landval) sv%active = .false.
    end if
                
    !marking tracers for deletion because they are old
    if (Globals%simdefs%tracerMaxAge > 0) then
        where(sv%state(:,col_age) >= Globals%simdefs%tracerMaxAge) sv%active = .false.
    end if
                
    !update resolution proxy
    col_res = Utils%find_str(var_name, Globals%Var%resolution,.true.)
    sv%resolution = var_dt(:,col_res)
    deallocate(var_name)
    deallocate(var_dt)
    
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
    integer :: nf, i, col_u, col_dwz, col_v, col_w, col_bat, part_idx
    real(prec), dimension(:,:), allocatable :: var_dt, var_hor_dt
    type(string), dimension(:), allocatable :: var_name, var_name_hor
    type(string), dimension(:), allocatable :: requiredVars, requiredHorVars
    real(prec), dimension(size(sv%state,1), size(sv%state,2)) :: LagrangianKinematic
    real(prec) :: VonKarman = 0.4
    real(prec) :: Hmin_Chezy = 0.1
    real(prec), dimension(size(sv%state,1)) :: chezyZ, dist2bottom
    real(8), dimension(size(sv%state,1)) :: aux_r8
    real(prec) :: threshold_bot_wat, landIntThreshold
    type(string) :: tag
    !-------------------------------------------------------------------------------------
    allocate(requiredVars(3))
    requiredVars(1) = Globals%Var%u
    requiredVars(2) = Globals%Var%v
    requiredVars(3) = Globals%Var%w
    
    allocate(requiredHorVars(3))
    requiredHorVars(1) = Globals%Var%u
    requiredHorVars(2) = Globals%Var%v
    requiredHorVars(3) = Globals%Var%w
    !requiredHorVars(4) = Globals%Var%ssh
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
    LagrangianKinematic = 0.0
    !Correct bottom values
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredHorVars, var_hor_dt, var_name_hor, reqVertInt = .false.)
    
    col_u = Utils%find_str(var_name_hor, Globals%Var%u, .true.)
    col_v = Utils%find_str(var_name_hor, Globals%Var%v, .true.)
    col_w = Utils%find_str(var_name_hor, Globals%Var%w, .false.)
    !col_ssh = Utils%find_str(var_name_hor, Globals%Var%ssh, .false.)
    col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
    col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    nf_u = Utils%find_str(var_name, Globals%Var%u, .true.)
    nf_v = Utils%find_str(var_name, Globals%Var%v, .true.)
    
    threshold_bot_wat = (Globals%Mask%waterVal + Globals%Mask%bedVal) * 0.5
    landIntThreshold = -0.98
    
    !if (col_ssh /= MV_INT) then
    !    !h = var_dt(:,col_bat)
    !    part_depth = sv%state(:,3)
    !else
    !    h = var_dt(:,col_bat) + var_hor_dt(:,col_ssh)
    !    part_depth = sv%state(:,3) + var_hor_dt(:,col_ssh)
    !end if
    
    !(Depth - Bathymetry)/(Bathymetry - (bathymetry - dwz) - dwz is a positive number and the rest are negative
    dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
    
    !Need to do this double check because landIntMask does not work properly when tracers are near a bottom wall
    !(interpolation gives over 1 but should still be bottom)
    where (dist2bottom < threshold_bot_wat)
        aux_r8 = max((sv%state(:,col_dwz)/2),Hmin_Chezy) / Globals%Constants%Rugosity
        chezyZ = (VonKarman / dlog(aux_r8))**2
        sv%state(:,4) = var_hor_dt(:,col_u) * chezyZ
        sv%state(:,5) = var_hor_dt(:,col_v) * chezyZ
    end where
    
    tag = 'particulate'
    part_idx = Utils%find_str(sv%varName, tag, .true.)
    
    where ((dist2bottom < landIntThreshold) .and. (sv%state(:,part_idx) == 1))
        !At the bottom and tracer is particulate
        LagrangianKinematic(:,1) = 0
        LagrangianKinematic(:,2) = 0
    elsewhere (dist2bottom < threshold_bot_wat)
        LagrangianKinematic(:,1) = Utils%m2geo(sv%state(:,4), sv%state(:,2), .false.)
        LagrangianKinematic(:,2) = Utils%m2geo(sv%state(:,5), sv%state(:,2), .true.)
    elsewhere
        LagrangianKinematic(:,1) = Utils%m2geo(var_dt(:, nf_u), sv%state(:,2), .false.)
        LagrangianKinematic(:,2) = Utils%m2geo(var_dt(:, nf_v), sv%state(:,2), .true.)
        sv%state(:,4) = var_dt(:,nf_u)
        sv%state(:,5) = var_dt(:,nf_v)
    end where

    nf = Utils%find_str(var_name, Globals%Var%w, .false.)
    if ((nf /= MV_INT) .and. (Globals%SimDefs%VerticalVelMethod == 1)) then
        !Make the vertical velocity 0 at the bottom.
        where ((dist2bottom < landIntThreshold) .and. (sv%state(:,part_idx) == 1))
            LagrangianKinematic(:,3) = 0
            sv%state(:,6) = 0
        elsewhere (dist2bottom < threshold_bot_wat)
            !Reduce velocity towards the bottom following a vertical logaritmic profile
            LagrangianKinematic(:,3) = var_hor_dt(:,col_w) * chezyZ
            sv%state(:,6) = LagrangianKinematic(:,3)
        elsewhere
            LagrangianKinematic(:,3) = var_dt(:, nf)
            sv%state(:,6) = var_dt(:, nf)
        end where
    else if ((nf /= MV_INT) .and. (Globals%SimDefs%VerticalVelMethod == 2)) then
        LagrangianKinematic(:,3) = VerticalMotion%Divergence(sv, bdata, time)
        sv%state(:,6) = LagrangianKinematic(:,3)
    else if ((nf == MV_INT) .or. (Globals%SimDefs%VerticalVelMethod == 3)) then
        LagrangianKinematic(:,3) = 0.0
        sv%state(:,6) = 0.0
    end if
    deallocate(var_dt)
    deallocate(var_hor_dt)
    deallocate(var_name_hor)
    deallocate(var_name)
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
    integer :: np, part_idx, col_dwz, col_bat
    real(prec), dimension(size(sv%state,1)) :: dist2bottom
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: DiffusionMixingLength
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v, rand_vel_w
    type(string) :: tag
    real(prec) :: landIntThreshold
    !Begin---------------------------------------------------------------------------
    landIntThreshold = -0.98
    tag = 'particulate'
    part_idx = Utils%find_str(sv%varName, tag, .true.)
    
    DiffusionMixingLength = 0.0
    if (Globals%Constants%DiffusionCoeff == 0.0) return
    !interpolate each background

    np = size(sv%active) !number of Tracers
    allocate(rand_vel_u(np), rand_vel_v(np), rand_vel_w(np))
    call random_number(rand_vel_u)
    call random_number(rand_vel_v)
    call random_number(rand_vel_w)
    
    !if we are still in the same path, use the same random velocity, do nothing
    !if we ran the path, new random velocities are generated and placed
    where ((sv%state(:,10) > 2.0*sv%resolution) .and. (sv%landIntMask < Globals%Mask%landVal))
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

    !update used mixing length
    DiffusionMixingLength(:,10) = sqrt(sv%state(:,4)*sv%state(:,4) + sv%state(:,5)*sv%state(:,5) + sv%state(:,6)*sv%state(:,6))
                
    !Deposited particles should not move
    if (any(sv%state(:,part_idx) == 1)) then
        col_dwz = Utils%find_str(sv%varname, Globals%Var%dwz, .true.)
        col_bat = Utils%find_str(sv%varname, Globals%Var%bathymetry, .true.)
        !(Depth - Bathymetry)/(Bathymetry - (bathymetry - dwz) - dwz is a positive number and the rest are negative
        dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
                    
        !if a single particle is particulate, check if they are at the bottom and don't move them if true.
        where ((dist2bottom < landIntThreshold) .and. (sv%state(:,part_idx) == 1))
            !update system positions
            DiffusionMixingLength(:,1) = 0
            DiffusionMixingLength(:,2) = 0
            DiffusionMixingLength(:,3) = 0
            DiffusionMixingLength(:,10) = 0
        end where
    end if
    deallocate(rand_vel_u, rand_vel_v, rand_vel_w)

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
    
    call KernelUtils%initialize() 
    
    end subroutine initKernel

    end module kernel_mod