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
	use, intrinsic :: ieee_arithmetic
	use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod
    use kernelUtils_mod
    
    use kernelLitter_mod
    use kernelVerticalMotion_mod
    use kernelColiform_mod
    use kernelDetritus_mod
    use kernelMOHIDWaterQuality_mod

    type :: kernel_class        !< Kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernel
    procedure :: run => runKernel
    procedure, private :: setCommonProcesses
    procedure, private :: interpolate_backgrounds
    procedure, private :: distance2bottom
    procedure, private :: LagrangianKinematic
    procedure, private :: DiffusionMixingLength
    procedure, private :: DiffusionIsotropic
    procedure, private :: StokesDrift
    procedure, private :: Windage
    procedure, private :: Beaching
    procedure, private :: FreeLitterAtBeaching
    procedure, private :: Aging
	procedure		   :: LagrangianVelModification																
    end type kernel_class

    type(kernelLitter_class) :: Litter       !< litter kernels
    type(kernelVerticalMotion_class) :: VerticalMotion   !< VerticalMotion kernels
    type(kernelColiform_class) :: Coliform !< coliform kernels
    type(kernelDetritus_class) :: Detritus !< coliform kernels
    type(kernelMOHIDWaterQuality_class) :: MOHIDWaterQuality   !< VerticalMotion kernels
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
    integer :: i
	logical, save :: printed = .false.
    !write(*,*)"Entrada Run Kernel tamanho bdata =", size(bdata)
    !do i = 1, size(bdata)
    !    write(*,*)"Tamanho background i = ", i, bdata(i)%fields%getSize()
    !enddo
    !running preparations for kernel lanch
    
    call self%setCommonProcesses(sv, bdata, time)
    !write(*,*)"Entrada interpolate_backgrounds"
    call self%interpolate_backgrounds(sv, bdata, time)
    !Computes distance to bottom for all tracers
    !write(*,*)"Entrada distance2bottom"
    call self%distance2bottom(sv)
    !running kernels for each type of tracer
    !write(*,*)"Entrada kernels"
	
	! Modify velocities based on Reichardt and log law for the inner turbulent layer near to the seabed.
	call self%LagrangianVelModification(sv, bdata, time)
	
    if (sv%ttype == Globals%Types%base) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv)
    else if (sv%ttype == Globals%Types%paper) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + Litter%DegradationLinear(sv) + VerticalMotion%Buoyancy(sv, bdata, time) + &
                    VerticalMotion%Resuspension(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + Litter%DegradationLinear(sv) + Litter%BioFouling(sv, dt) + VerticalMotion%Buoyancy(sv, bdata, time) + &
                    VerticalMotion%Resuspension(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%coliform) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + coliform%MortalityT90(sv, bdata, time) + coliform%Dilution(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%seed) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%DiffusionMixingLength(sv, bdata, time, dt) + self%Aging(sv) + &
                    VerticalMotion%Buoyancy(sv, bdata, time) + VerticalMotion%Resuspension(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%detritus) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%DiffusionMixingLength(sv, bdata, time, dt) + self%Aging(sv) + &
                    VerticalMotion%Buoyancy(sv, bdata, time) + VerticalMotion%Resuspension(sv, bdata, time, dt) + &
                    detritus%Degradation(sv, dt)
    else if (sv%ttype == Globals%Types%WaterQuality) then
        runKernel = self%LagrangianKinematic(sv, bdata, time) + self%StokesDrift(sv, bdata, time) + &
                    self%Windage(sv, bdata, time) + self%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%Aging(sv) + MOHIDWaterQuality%WQProcess(sv, bdata, time, dt) + MOHIDWaterQuality%Dilution(sv, bdata, time, dt)
    end if
    if (Globals%simDefs%FreeLitterAtBeaching == 1) then
        runKernel = self%FreeLitterAtBeaching(sv, bdata, time, runKernel, dt)
    else
        runKernel = self%Beaching(sv, runKernel)
    endif
    
    runKernel = VerticalMotion%CorrectVerticalBounds(sv, runKernel, bdata, dt)
    	
	
    end function runKernel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com
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
    integer :: i, j, col_age, col_bat, col_bat_sv, col_landintmask, col_res, col_ssh, col_DifVelStdr
	integer :: col_rugosityVar,col_rugosityVar_sv
	integer :: col_D50Var,col_D50Var_sv
	integer :: counterr
    real(prec) :: maxLevel(2)
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    type(string) :: tag
    logical bottom_emmission
    !-----------------------------------------------------------
    !write(*,*)"Entrada setCommonProcesses"
    allocate(requiredVars(6))
    requiredVars(1) = Globals%Var%landIntMask
    requiredVars(2) = Globals%Var%resolution
    requiredVars(3) = Globals%Var%bathymetry
    requiredVars(4) = Globals%Var%ssh
    requiredVars(5) = Globals%Var%rugosityVar
    requiredVars(6) = Globals%Var%D50Var
	
    !write(*,*)"Entrada setCommonProcesses interpolate"
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true.)
    
    !write(*,*)"Saida setCommonProcesses interpolate"
    bottom_emmission = .false.
    col_bat = Utils%find_str(var_name, Globals%Var%bathymetry, .false.)
    !Set tracers bathymetry
    col_bat_sv = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
    if (col_bat /= MV_INT) then
        sv%state(:,col_bat_sv) = var_dt(:,col_bat)
    else
        sv%state(:,col_bat_sv) = 0.0
    endif

    !set tracer bottom rugosity
    col_rugosityVar = Utils%find_str(var_name, Globals%Var%rugosityVar, .false.)
    col_rugosityVar_sv = Utils%find_str(sv%varName, Globals%Var%rugosityVar, .true.)
    sv%state(:,col_rugosityVar_sv) = max(var_dt(:,col_rugosityVar), 0.0001)

    !set tracer bottom D50
    col_D50Var = Utils%find_str(var_name, Globals%Var%D50Var, .false.)
    col_D50Var_sv = Utils%find_str(sv%varName, Globals%Var%D50Var, .true.)
    sv%state(:,col_D50Var_sv) = max(var_dt(:,col_D50Var), 0.0001)
	
!	counterr = 0
!	do i= 1, size(sv%state,1)
!		if (mod(counterr, 10) == 0) then
!			write(*,'(  A5, A12,, A12)') , " Id:", 'rugosityVar', 'D50Var'
!			write(*,*),' '
!		end if
!		counterr = counterr + 1
!		write(*,'( I5, F12.4, F12.4)') , i, sv%state(i,col_rugosityVar_sv), sv%state(i,col_D50Var_sv)
!D	end do	 	


    tag = 'age'
    col_age = Utils%find_str(sv%varName, tag, .true.)
    
    !Check for particles below sea bottom
!    if (col_bat /= MV_INT) then
!        where (sv%state(:,3) < var_dt(:,col_bat)) sv%state(:,3) = var_dt(:,col_bat)
!    endif
    
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
    
    !TODO : If the hdf5 does not have ssh, should use the verticalZ (meaning we must save an extra 2D var with the hdf original var (cell faces)
    if (Globals%simDefs%inputFromHDF5) then
        col_ssh = Utils%find_str(var_name, Globals%Var%ssh, .false.)
        if (col_ssh /= MV_INT) then
            where (sv%state(:,3) >  var_dt(:,col_ssh)) sv%state(:,3) = var_dt(:,col_ssh) - 0.00001
        else
            !ssh not found... assume 0.0 as the limit
            where (sv%state(:,3) >  0.0) sv%state(:,3) = - 0.00001
        endif
        
    else
        maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .false.)   
        if (maxLevel(2) /= MV) where (sv%state(:,3) > maxLevel(2)) sv%state(:,3) = maxLevel(2)-0.00001  
    endif
    
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
        where(int(abs(sv%landintmask) + Globals%mask%landval*0.05) == Globals%mask%landval) sv%active = .false.
    end if
                
    !marking tracers for deletion because they are old
    if (Globals%simdefs%tracerMaxAge > 0) then
        where(sv%state(:,col_age) >= Globals%simdefs%tracerMaxAge) sv%active = .false.
    end if
    
    !update resolution proxy
    col_res = Utils%find_str(var_name, Globals%Var%resolution,.true.)
    sv%resolution = var_dt(:,col_res)
    
    !Diffusion processes. Initialize diffusion velocity standard deviation to 0.
    if (Globals%SimDefs%DiffusionMethod == 2) then !SullivanAllen
        tag = 'VelStandardDeviation'
        col_DifVelStdr = Utils%find_str(sv%varName, tag, .true.)
        sv%state(:,col_DifVelStdr) = 0.0
    endif
        
    deallocate(var_name)
    deallocate(var_dt)
    !write(*,*)"Saida setCommonProcesses"
    end subroutine setCommonProcesses
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> interpolates background variables to tracers positions
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    subroutine interpolate_backgrounds(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: col_temp, col_sal, col_temp_sv, col_sal_sv, col_dwz, col_dwz_sv
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    !-----------------------------------------------------------
    
    if (sv%ttype == Globals%Types%base) then
        allocate(requiredVars(1))
        requiredVars(1) = Globals%Var%dwz
    else
        allocate(requiredVars(3))
        requiredVars(1) = Globals%Var%temp
        requiredVars(2) = Globals%Var%sal
        requiredVars(3) = Globals%Var%dwz
    endif

    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true.)
    
    !Set tracers dwz
    col_dwz = Utils%find_str(var_name, Globals%Var%dwz, .true.)
    col_dwz_sv = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    
    sv%state(:,col_dwz_sv) = var_dt(:,col_dwz)
    
    !Set tracers temperature
    !Not usable for dilution of temperature so will need to be changed in the future (for example save in a ambient_temp sv name)
    if (sv%ttype /= Globals%Types%base) then
        col_temp = Utils%find_str(var_name, Globals%Var%temp, .false.)
        col_temp_sv = Utils%find_str(sv%varName, Globals%Var%temp, .false.)
        col_sal = Utils%find_str(var_name, Globals%Var%sal, .false.)
        col_sal_sv = Utils%find_str(sv%varName, Globals%Var%sal, .false.)
        if (col_temp /= MV_INT .and. col_temp_sv /= MV_INT) then
           sv%state(:,col_temp_sv) = var_dt(:,col_temp) 
        endif
        
        if (col_sal /= MV_INT .and. col_sal_sv /= MV_INT) then
            sv%state(:,col_sal_sv) = var_dt(:,col_sal)
        endif
    endif
    
    deallocate(var_name)
    deallocate(var_dt)
    end subroutine interpolate_backgrounds
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> computes distance to bottom for all tracers (beware! does not yet consider water level)
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    subroutine distance2bottom(self, sv)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    integer :: col_dwz, col_bat, col_dist2bottom, i
    type(string) :: tag
    !-----------------------------------------------------------
    !Set tracers dwz
    col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
    tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)
    !Need to add water elevation into account because in low depth areas the result will be wrong
    sv%state(:,col_dist2bottom) = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
    
    end subroutine distance2bottom
    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com
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
    integer :: nf_w, part_idx, col_dist2bottom																												
    real(prec), dimension(size(sv%state,1), size(sv%state,2)) :: LagrangianKinematic
    real(prec), dimension(size(sv%state,1)) :: dist2bottom											  
    real(prec) :: threshold_bot_wat, landIntThreshold
	real(prec), dimension(size(sv%state,1)) :: Threshold_value
	real(prec), dimension(size(sv%state,1)) :: LandIntThreshold_value
    type(string) :: tag
    integer :: i,counterr
    !-------------------------------------------------------------------------------------
    !write(*,*)"Entrada kinematic"
    col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
	tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)
    dist2bottom = sv%state(:,col_dist2bottom)
	
    tag = 'particulate'
    part_idx = Utils%find_str(sv%varName, tag, .true.)

    ! threshold_bot_wat: equal to 0.0 will be correspond to the distance to the bottom equal to the last measurmet(given data).
    ! LandIntThresholdValue: will be correspond to the distance to the bottom which the volcity could be considered zero
	! LandIntThresholdValue gives the value the dist2bottom for this threshold.
	! u, v and w velocities for distance to the seabed(h(i) - bathymetry) < Rugosity reach to zero.    
		
	threshold_bot_wat = 0.0	
	! Threshold_value: distance from the bottom (seabed) in unit [meter].
	! It could be a constant * Globals%Constants%Rugosity
	Threshold_value = Globals%Constants%BedLoadThickness
	LandIntThreshold_value = VerticalMotion%LandIntThresholdValue(sv, bdata, time, Threshold_value)	
	
	LagrangianKinematic = 0.0
 
	where (dist2bottom < LandIntThreshold_value .and. (sv%state(:,part_idx) == 1))
		LagrangianKinematic(:,1) = 0.0
		LagrangianKinematic(:,2) = 0.0
		
	else where
		LagrangianKinematic(:,1) = Utils%m2geo(sv%state(:,4), sv%state(:,2), .false.)
		LagrangianKinematic(:,2) = Utils%m2geo(sv%state(:,5), sv%state(:,2), .true.)
	end where

	nf_w = Utils%find_str(sv%varName, Globals%Var%w, .false.)

!	New version: log law
    if ((nf_w /= MV_INT) .and. (Globals%SimDefs%VerticalVelMethod == 1)) then
        !Lagrangian dispacementin vertical direction
		where (dist2bottom < LandIntThreshold_value .and. (sv%state(:,part_idx) == 1))
			! It is assumed that there is not vertical Lagrangian dispacement for the particles lower than LandIntThreshold_value near to the seabed,
			LagrangianKinematic(:,3) = 0.0
	
		else where
		! Otherwise, it is assumed that the vertical Lagrangian dispacement for the particel is according to the velocity profile.
		! We should take into account based on LagrangianVelModification function, some partilces follows log-law.
			LagrangianKinematic(:,3) = sv%state(:,6)

		end where
		
    else if ((nf_w /= MV_INT) .and. (Globals%SimDefs%VerticalVelMethod == 2)) then
        LagrangianKinematic(:,3) = VerticalMotion%Divergence(sv, bdata, time)
        sv%state(:,6) = LagrangianKinematic(:,3)
    else if ((nf_w == MV_INT) .or. (Globals%SimDefs%VerticalVelMethod == 3)) then
        LagrangianKinematic(:,3) = 0.0
        sv%state(:,6) = 0.0
    end if
    
    !write(*,*)"Saida kinematic"

    end function LagrangianKinematic

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Computes the influence of wave velocity in tracer kinematics. computes wave velocity if not present in nc or hdf
    !> @param[in] self, sv, bdata, time
    function StokesDrift(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: col_vsdx, col_vsdy, col_hs, col_ts, col_wd, col_wl, col_bat, col_ssh, col_DifVelStdr
    integer :: bkg
    real(prec) :: waveCoeff
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: StokesDrift
    real(prec), dimension(size(sv%state,1)) :: depth, WaterDepth, WaveAmplitude
    real(prec), dimension(size(sv%state,1)) :: AngFrequency, WaveLength, WaveNumber, VelStokesDrift
    real(prec), dimension(size(sv%state,1)) :: C_Term, F_Aux, G_Aux
    real(prec) :: maxLevel(2)
    real(prec) :: Pi = 4*atan(1.0)
    type(string) :: tag
    !Begin--------------------------------------------------------------------
    allocate(requiredVars(7))
    requiredVars(1) = Globals%Var%vsdx
    requiredVars(2) = Globals%Var%vsdy
    requiredVars(3) = Globals%Var%hs !wave height
    requiredVars(4) = Globals%Var%ts !wave period
    requiredVars(5) = Globals%Var%wd !wave direction
    requiredVars(6) = Globals%Var%wl !wave lenght
    requiredVars(7) = Globals%Var%ssh !water level
    waveCoeff = 0.01
    StokesDrift = 0.0
    
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true., reqVertInt = .false.)

    col_vsdx = Utils%find_str(var_name, Globals%Var%vsdx, .false.)
    col_vsdy = Utils%find_str(var_name, Globals%Var%vsdy, .false.)
    col_hs = Utils%find_str(var_name, Globals%Var%hs, .false.)
    col_ts = Utils%find_str(var_name, Globals%Var%ts, .false.)
    col_wd = Utils%find_str(var_name, Globals%Var%wd, .false.)
    col_wl = Utils%find_str(var_name, Globals%Var%wl, .false.)
    col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
    col_ssh = Utils%find_str(var_name, Globals%Var%ssh, .false.)
    
    if (Globals%SimDefs%DiffusionMethod == 2) then !SullivanAllen
        tag = "VelStandardDeviation"
        col_DifVelStdr = Utils%find_str(sv%varName, tag, .true.)
    endif
    
    if (col_vsdx /= MV_INT .and. col_vsdy /= MV_INT) then
        !computing the depth weight
        depth = sv%state(:,3)
        where (depth>=0.0) depth = 0.0
        depth = exp(depth)
        !write dx/dt
        where(abs(sv%landIntMask) < Globals%Mask%landVal)
            StokesDrift(:,1) = Utils%m2geo(var_dt(:,col_vsdx), sv%state(:,2), .false.)*waveCoeff*depth
            StokesDrift(:,2) = Utils%m2geo(var_dt(:,col_vsdy), sv%state(:,2), .true.)*waveCoeff*depth
        endwhere

        if (Globals%SimDefs%DiffusionMethod == 2) then !SullivanAllen
            sv%state(:,col_DifVelStdr) = sv%state(:,col_DifVelStdr) + SQRT(var_dt(:,col_vsdx)**2 + var_dt(:,col_vsdy)**2)
        endif
        
    elseif (col_hs /= MV_INT .and. col_ts /= MV_INT .and. col_wd /= MV_INT) then
        !Could not find stokes velocity (try computing from other wave parameters)
        
        if (col_ssh /= MV_INT) then
            ! get water level and use it to compute particle depth
            !TODO - make sure the signs are correct
            depth = max(var_dt(:,col_ssh) - sv%state(:,3), 0.0)
        else
            maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .false.) 
            !use maxlevel
            if (maxLevel(2) /= MV) then
                !TODO: check if signs are correct
                depth = max(maxLevel(2) - sv%state(:,3), 0.0)
            else
                depth = 0.0
            endif
        endif
        
        !WavePeriodAux   = max(var_dt(:,col_ts), 0.01) !only needed if we want to add the uncertainty method of Rod
        AngFrequency    = 2.0 * Pi / var_dt(:,col_ts)

        !Depth               = CurrentPartic%Position%Z-             &
        !                        Me%EulerModel(emp)%SZZ(i, j, WS_KUB)
                                              
        !if (Depth < 0.)       Depth = 0. 
            
        WaterDepth     = var_dt(:,col_ssh) - sv%state(:,col_bat)
        WaveAmplitude  = var_dt(:,col_hs) / 2.
                        
        if (col_wl == MV_INT) then
            !Compute wave lenght from wave period
            where (var_dt(:,col_ts) < 1e-3) 
                WaveLength = 0.
                G_Aux = 0.0
                F_Aux = 0.0
            elsewhere
                G_Aux = ((2 * Pi / var_dt(:,col_ts))**2) * WaterDepth / 9.81
                F_Aux = G_Aux + (1 / (1. + 0.6522 * G_Aux + 0.4622 * (G_Aux**2) + &
                        0.0864 * (G_Aux**4) + 0.0675 * (G_Aux**5)))
                
                where (F_Aux > 0. .and. WaterDepth > 0.)
                    WaveLength = var_dt(:,col_ts) * sqrt(9.81 * WaterDepth / F_Aux)
                elsewhere
                    WaveLength = 0.
                endwhere
            endwhere
        else
            WaveLength = var_dt(:,col_wl)
        endif                      
        
        where (WaveLength > 1e-3)
            WaveNumber = max(2.0 * Pi / WaveLength, 0.0)
        elsewhere
            WaveNumber = 0.0
        endwhere
        
        
        where (WaveNumber == 0.0 .or. WaterDepth < 0.01)
            VelStokesDrift = 0.0
        elsewhere (WaterDepth > WaveLength / 2.0)
            !use longuetHiggins Deep
            !C_Term = - (WaveAmplitude**2 * AngFrequency * sinh(2.0 * WaveNumber * WaterDepth)) / &
            !                    ( 4 * WaterDepth * sinh(WaveNumber * WaterDepth)**2)
            VelStokesDrift  = WaveAmplitude**2 * AngFrequency * WaveNumber * exp(-2* WaveNumber * Depth)        

        elsewhere
            
            VelStokesDrift = WaveAmplitude**2 * AngFrequency * WaveNumber * (( cosh(2 * WaveNumber * (Depth - WaterDepth)) ) /         &
                             ( 2 * (sinh(WaveNumber * WaterDepth)* sinh(WaveNumber * WaterDepth)))) - (WaveAmplitude**2 * AngFrequency * sinh(2.0 * WaveNumber * WaterDepth)) / &
                             ( 4 * WaterDepth * sinh(WaveNumber * WaterDepth)**2)
        endwhere
        
        where (VelStokesDrift > 10.0) VelStokesDrift = 0.0
        
        StokesDrift(:,1) = Utils%m2geo(cos(var_dt(:,col_wd) * (Pi / 180.)) * VelStokesDrift, sv%state(:,2), .false.)
        StokesDrift(:,2) = Utils%m2geo(sin(var_dt(:,col_wd) * (Pi / 180.)) * VelStokesDrift, sv%state(:,2), .true.)
        
        if (Globals%SimDefs%DiffusionMethod == 2) then !SullivanAllen
            sv%state(:,col_DifVelStdr) = sv%state(:,col_DifVelStdr) + VelStokesDrift
        endif
        
    endif
    
    deallocate(var_name)
    deallocate(var_dt)
    

    end function StokesDrift

    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Computes the influence of wave velocity in tracer kinematics
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    !function StokesDrift(self, sv, bdata, time)
    !class(kernel_class), intent(inout) :: self
    !type(stateVector_class), intent(in) :: sv
    !type(background_class), dimension(:), intent(in) :: bdata
    !real(prec), intent(in) :: time
    !integer :: np, nf, bkg
    !real(prec) :: waveCoeff
    !real(prec), dimension(:,:), allocatable :: var_dt
    !type(string), dimension(:), allocatable :: var_name
    !type(string), dimension(:), allocatable :: requiredVars
    !real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: StokesDrift
    !real(prec), dimension(size(sv%state,1)) :: depth
    !allocate(requiredVars(2))
    !requiredVars(1) = Globals%Var%vsdx
    !requiredVars(2) = Globals%Var%vsdy
    !waveCoeff = 0.01
    !StokesDrift = 0.0
    !!interpolate each background
    !do bkg = 1, size(bdata)
    !    if (bdata(bkg)%initialized) then
    !        if(bdata(bkg)%hasVars(requiredVars)) then
    !            np = size(sv%active) !number of Tracers
    !            nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
    !            allocate(var_dt(np,nf))
    !            allocate(var_name(nf))
    !            !interpolating all of the data
    !            call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name)
    !            !computing the depth weight
    !            depth = sv%state(:,3)
    !            where (depth>=0.0) depth = 0.0
    !            depth = exp(depth)
    !            !write dx/dt
    !            nf = Utils%find_str(var_name, Globals%Var%vsdx, .true.)
    !            where(sv%landIntMask < Globals%Mask%landVal) StokesDrift(:,1) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .false.)*waveCoeff*depth
    !            nf = Utils%find_str(var_name, Globals%Var%vsdy, .true.)
    !            where(sv%landIntMask < Globals%Mask%landVal) StokesDrift(:,2) = Utils%m2geo(var_dt(:,nf), sv%state(:,2), .true.)*waveCoeff*depth
    !            deallocate(var_name)
    !            deallocate(var_dt)
    !        end if
    !    end if
    !end do
    !
    !end function StokesDrift
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC. Revision 30-08-2025 by Joao Sobrinho
    !> @brief
    !> Computes the influence of wind velocity in tracer kinematicstokes
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Windage(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: np, nf, nf2, nf3, bkg
    real(prec) :: windCoeff
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Windage
    real(prec), dimension(size(sv%state,1)) :: depth
    real(prec) :: maxLevel(2)
    integer                                 :: col_ssh, col_DifVelStdr, col_u10, col_v10
    type(string)                            :: tag
    !Begin-------------------------------------------------------------------
    
    allocate(requiredVars(3))
    requiredVars(1) = Globals%Var%u10
    requiredVars(2) = Globals%Var%v10
    requiredVars(3) = Globals%Var%ssh
    
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true., reqVertInt = .false.)
    
    windCoeff = Globals%Constants%WindDragCoeff
    Windage = 0.0
    
    col_ssh = Utils%find_str(var_name, Globals%Var%ssh, .false.)
    col_u10 = Utils%find_str(var_name, Globals%Var%u10, .false.)
    col_v10 = Utils%find_str(var_name, Globals%Var%v10, .false.)
    
    if (col_u10 == MV_INT .or. col_v10 == MV_INT) then
        !no wind available, skip windage
        return
    endif
    
    if (col_ssh /= MV_INT) then
        !water level exists, use it to get a more accurate solution
        depth = max(var_dt(:,col_ssh) - sv%state(:,3), 0.0)
                    
        !if tracer is less than 5 cm below ssh, consider full wind effect
        where (depth >= 0.05) depth = 0.0
    else
        ! use what is available (depth levels probably)
        maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .true.)
        if (maxLevel(2) /= MV) where (depth >= 0 .or. depth >= maxLevel(2) - 0.05) depth = 0.0
    endif
    
    !This wont do much... unless we actually start considering a depth profile from water level to the input wind data's altitude
    depth = exp(10.0*depth)
                
    !write dx/dt 
    where(abs(sv%landIntMask) < Globals%Mask%landVal)
        Windage(:,1) = Utils%m2geo(var_dt(:,col_u10), sv%state(:,2), .false.)*windCoeff * depth
        Windage(:,2) = Utils%m2geo(var_dt(:,col_v10), sv%state(:,2), .true.)*windCoeff * depth
    endwhere
    
    if (Globals%SimDefs%DiffusionMethod == 2) then !SullivanAllen
        tag = "VelStandardDeviation"
        col_DifVelStdr = Utils%find_str(sv%varName, tag, .true.)
        sv%state(:,col_DifVelStdr) = sv%state(:,col_DifVelStdr) + SQRT(Windage(:,1)**2 + Windage(:,2)**2) * Globals%Constants%VarVelHX
    endif
        
    deallocate(var_name)
    deallocate(var_dt)
    
    end function Windage


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC 
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.09.12 | Email:shabani.mohsen@outlook.com	
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
    real(prec), dimension(size(sv%state,1)) :: beachCoeffRand, beachCoeffRandC
    real(prec), dimension(size(sv%state,1)) :: beachWeight
    real(prec) :: lbound, ubound
    integer :: i
	
    beachCoeff = 1.0
	
	! Create  random values between [0,1) to compare with Beaching probability
    !call random_number(beachCoeffRand) !this is a uniform distribution generator
    !beachCoeffRand = max(0.0, beachCoeffRand - Globals%Constants%BeachingStopProb)  !clipping the last % to zero
    !where (beachCoeffRand /= 0.0 ) beachCoeffRand = beachCoeffRand*(1.0/maxval(beachCoeffRand))

	! Create  random values between [0,1) to compare with Beaching probability
    call random_number(beachCoeffRand) !this is a uniform distribution generator 
	beachCoeffRandC = max(0.0, Globals%Constants%BeachingStopProb - beachCoeffRand)  !Set random numbers greater than the probability to zero.(beaching will not happen)
    !where (beachCoeffRandC /= 0.0 ) beachCoeffRandC = beachCoeffRandC*(1.0/maxval(beachCoeffRandC))	!Normalize the values between [0,1]
    where (beachCoeffRandC /= 0.0 ) beachCoeffRandC = max(1.0,	beachCoeffRandC)!Normalize the values between [0,1]
	!beachCoeffRandC = 1.0
    Beaching = svDt

    if (Globals%Constants%BeachingStopProb /= 0.0) then !beaching is completely turned off if the stopping propability is zero

        !getting the bounds for the interpolation of the land interaction field that correspond to beaching
        !lbound = (Globals%Mask%beachVal + Globals%Mask%waterVal)*0.5
        !ubound = (Globals%Mask%beachVal + Globals%Mask%landVal)*0.5
        lbound =  Globals%Mask%landVal * (-1.0)
		ubound =  Globals%Mask%landVal * (+1.0)
        !beachWeight = 1 - 0.5*(sv%landIntMask - lbound)/(ubound-lbound) !linear distance weight for beaching
        !beachWeight = 1 - 0.9*(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound) !quadratic weight
        ! 2nd order weight
		!beachWeight = 1 - 1.0 *(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound) !quadratic weight
        ! 4th order weight [0.5, 2]
		!beachWeight = 1 - 1.0 *(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound) !4th order weight
		! 6th order weight axisymmetric respect to 0 in domain of [-2,2]
		beachWeight = 1 - 1.0 *(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound) *(sv%landIntMask - lbound)/(ubound-lbound)*(sv%landIntMask - lbound)/(ubound-lbound)!6th order weight
		
        !replacing 1.0 with a coefficient from beaching where needed
        where((sv%landIntMask <= ubound) .and. (sv%landIntMask >= lbound)) beachCoeff = (beachCoeffRandC) * beachWeight
		!where((sv%landIntMask > ubound)) beachCoeff= 0.0
		do i=1,3
            Beaching(:,i) = svDt(:,i)*beachCoeff !position derivative is affected
            sv%state(:,i+3) = sv%state(:,i+3)*beachCoeff !so are the velocities
        end do
    end if   

    end function Beaching
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Beaching Kernel, uses the already updated state vector and determines if
    !> and how beaching occurs using the formulation devised in Project FreeLitterAt. Affects the state vector and state vector derivative.
    !> @param[in] self, sv, svDt, dt
    !---------------------------------------------------------------------------
    function FreeLitterAtBeaching(self, sv, bdata, time, svDt, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)), intent(in) :: svDt
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: FreeLitterAtBeaching
    real(prec), intent(in) :: dt
    real(prec)  :: limLeft, limRight, limBottom, limTop, x, y, Bathymetry
    real(prec) :: WaterColumn, WaterLevel, rand1, Probability, Tbeach, Tunbeach, Hs, Tp, m, threshold
    type(vector), dimension(2) :: beachPolygonBBox !Vertices of the beching area polygon
    type(vector), dimension(:), allocatable :: beachPolygonVertices !Vertices of the beching area polygon
    integer :: col_bat, col_Hs, col_Tp, col_beachPeriod, col_beachedWaterLevel, col_ssh, col_beachAreaId
    type(string) :: tag, outext
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    integer :: i, np, j
    type(string), dimension(:), allocatable :: requiredVars
    type(string) :: temp_str
    !Begin------------------------------------------------------
    
    ! Making FreeLitterAtBeaching = svDt so it can be used to set x and y position derivatives to 0 if beached
    FreeLitterAtBeaching = svDt
    
    !bail early
    if (size(Globals%BeachingAreas%beachArea) < 1) then
        outext = '[Kernel::FreeLitterAtBeaching] At least 1 beaching area polygon is required for beaching to be computed'
        call Log%put(outext)
        return
    endif
       
    allocate(requiredVars(3))
    requiredVars(1) = Globals%Var%ssh
    requiredVars(2) = Globals%Var%hs
    requiredVars(3) = Globals%Var%ts
    
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true., reqVertInt = .false.)
    
    col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)       
    col_ssh = Utils%find_str(var_name, Globals%Var%ssh, .true.)
    
    tag = 'beachPeriod'
    col_beachPeriod = Utils%find_str(sv%varName, tag, .true.)
    
    tag = 'beachAreaId'
    col_beachAreaId = Utils%find_str(sv%varName, tag, .true.)
    
    tag = 'beachedWaterLevel'
    col_beachedWaterLevel = Utils%find_str(sv%varName, tag, .true.)
    
    if (any(Globals%BeachingAreas%beachArea(:)%par%runUpEffect == 1)) then
        col_Hs = Utils%find_str(var_name, Globals%Var%hs, .true.) ! Significant Wave Height      
        col_Tp = Utils%find_str(var_name, Globals%Var%ts, .true.) ! Wave Period     
    endif
    
    do i=1,size(Globals%BeachingAreas%beachArea)
        
        beachPolygonBBox = Geometry%getBoundingBox(Globals%BeachingAreas%beachArea(i)%par%geometry)
        allocate(beachPolygonVertices(size(Geometry%getPoints(Globals%BeachingAreas%beachArea(i)%par%geometry))))
        beachPolygonVertices = Geometry%getPoints(Globals%BeachingAreas%beachArea(i)%par%geometry)
        
        limLeft = beachPolygonBBox(1)%x !limit left
        limBottom = beachPolygonBBox(1)%y !limit bottom
        limRight = beachPolygonBBox(2)%x !limit right
        limTop = beachPolygonBBox(2)%y !limit top
        
        Tbeach      = Globals%BeachingAreas%beachArea(i)%par%beachTimeScale
        Tunbeach    = Globals%BeachingAreas%beachArea(i)%par%unbeachTimeScale

        do np=1, size(sv%state,1)
                
            if (sv%state(np,col_beachPeriod) > 0.0 .and. sv%state(np,col_beachAreaId) == Globals%BeachingAreas%beachArea(i)%par%id) then
                ! Already beached in this beaching area, so try and unbeach it. no need to check if it is inside the beaching area
                WaterLevel = var_dt(np,col_ssh)
                Bathymetry = sv%state(np,col_bat)
            
                WaterColumn =  WaterLevel - Bathymetry
                
                sv%state(np,col_beachPeriod) = sv%state(np,col_beachPeriod) + dt
                threshold = Globals%BeachingAreas%beachArea(i)%par%waterColumnThreshold
                
                if (WaterColumn > threshold .and. WaterLevel > Bathymetry) then
                
                    call random_number(rand1)
                    
                    Probability = 1 - exp(-dt/Tunbeach)

                    if (Probability > rand1) then
                            
                        if (Globals%BeachingAreas%beachArea(i)%par%runUpEffectUnbeach == 1) then
                            Hs         = var_dt(np,col_Hs)
                            Tp         = var_dt(np,col_Tp) 
                            m          = Globals%BeachingAreas%beachArea(i)%par%beachSlope

                            WaterLevel = var_dt(np,col_ssh) + WaveRunUpStockdon2006(Hs,Tp,m)   
                        endif

                        if (WaterLevel > sv%state(np,col_beachedWaterLevel)) then
                            sv%state(np,col_beachPeriod) = 0.0
                        endif
                    endif
                endif
            
            else !Not beached yet, try beaching it. Need to check if tracer is inside beaching area
                x = sv%state(np,1)
                y = sv%state(np,2)
                !Skip all tracers outside the limits of this beaching area
                if (Utils%isPointInsidePolygon(x, y, beachPolygonVertices)) then
                    WaterLevel = var_dt(np,col_ssh)
                    Bathymetry = sv%state(np,col_bat)
                    WaterColumn =  WaterLevel - Bathymetry
                    threshold = Globals%BeachingAreas%beachArea(i)%par%waterColumnThreshold
                    if (WaterColumn < threshold .and. WaterLevel > Bathymetry) then
                        
                        call random_number(rand1)
                
                        Probability = 1 - exp(-dt/Tbeach)
                    
                        if (Probability > rand1) then 
                            if (Globals%BeachingAreas%beachArea(i)%par%runUpEffect == 1) then
                                Hs         = var_dt(np,col_Hs)
                                Tp         = var_dt(np,col_Tp) 
                                m          = Globals%BeachingAreas%beachArea(i)%par%beachSlope

                                sv%state(np,col_beachedWaterLevel) = var_dt(np,col_ssh) + WaveRunUpStockdon2006(Hs,Tp,m)  
                            else
                                sv%state(np,col_beachedWaterLevel) = var_dt(np,col_ssh)
                            endif
                                    
                            sv%state(np,col_beachPeriod) = sv%state(np,col_beachPeriod) + dt
                        endif
                    endif
                endif
            endif
        enddo
    enddo
        
    where (sv%state(:,col_beachPeriod) > 0.0)
        FreeLitterAtBeaching(:,1) = 0.0 !Do not change positions
        FreeLitterAtBeaching(:,2) = 0.0 !Do not change positions
        FreeLitterAtBeaching(:,3) = 0.0 !Do not change positions
        sv%state(:,4) = 0.0 !nor velocities
        sv%state(:,5) = 0.0 !nor velocities
        sv%state(:,6) = 0.0 !nor velocities
        sv%state(:,7) = 0.0 !nor velocities
        sv%state(:,8) = 0.0 !nor velocities
        sv%state(:,9) = 0.0 !nor velocities
    endwhere
    
    
    deallocate(var_dt)
    deallocate(var_name)
    
    end function FreeLitterAtBeaching

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
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com	
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
    integer :: np, part_idx, col_dist2bottom
    real(prec), dimension(size(sv%state,1)) :: dist2bottom
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: DiffusionMixingLength
    real(prec), dimension(:), allocatable :: rand_vel_u, rand_vel_v, rand_vel_w
    type(string) :: tag
    real(prec) :: landIntThreshold
	real(prec), dimension(size(sv%state,1)) :: Threshold_value
	real(prec), dimension(size(sv%state,1)) :: LandIntThreshold_value
	
    !Begin---------------------------------------------------------------------------
    
    if (Globals%SimDefs%DiffusionMethod == 2) then !SullivanAllen
        DiffusionMixingLength = SullivanAllen(self, sv, bdata, time, dt)
        return
    endif
    
    landIntThreshold = -0.98
    tag = 'particulate'
    part_idx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)
    
    DiffusionMixingLength = 0.0
    if (Globals%Constants%DiffusionCoeff == 0.0) return
    !interpolate each background


	! Threshold_value: distance from the bottom (seabed) in unit [meter].
	! It could be a constant * Globals%Constants%Rugosity
	Threshold_value = Globals%Constants%BedLoadThickness 
	LandIntThreshold_value = VerticalMotion%LandIntThresholdValue(sv, bdata, time, Threshold_value)
	

    np = size(sv%active) !number of Tracers
    allocate(rand_vel_u(np), rand_vel_v(np), rand_vel_w(np))
    call random_number(rand_vel_u)
    call random_number(rand_vel_v)
    call random_number(rand_vel_w)
    
    !if we are still in the same path, use the same random velocity, do nothing
    !if we ran the path, new random velocities are generated and placed
    where ((sv%state(:,10) > 2.0*sv%resolution) .and. (abs(sv%landIntMask) < Globals%Mask%landVal))
        !DiffusionMixingLength(:,7) = (2.*rand_vel_u-1.)*sqrt(Globals%Constants%DiffusionCoeff*abs(sv%state(:,4))/dt)/dt
        !DiffusionMixingLength(:,8) = (2.*rand_vel_v-1.)*sqrt(Globals%Constants%DiffusionCoeff*abs(sv%state(:,5))/dt)/dt
        !DiffusionMixingLength(:,9) = (2.*rand_vel_w-1.)*sqrt(0.000001*Globals%Constants%DiffusionCoeff*abs(sv%state(:,6))/dt)/dt
        DiffusionMixingLength(:,7) = (2.*rand_vel_u-1.)*sqrt(Globals%Constants%DiffusionCoeff/dt)/dt
        DiffusionMixingLength(:,8) = (2.*rand_vel_v-1.)*sqrt(Globals%Constants%DiffusionCoeff/dt)/dt
        DiffusionMixingLength(:,9) = (2.*rand_vel_w-1.)*sqrt(0.000001*Globals%Constants%DiffusionCoeff/dt)/dt
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
        dist2bottom = sv%state(:,col_dist2bottom)
        !if a single particle is particulate, check if they are at the bottom and don't move them if true.
        where ((dist2bottom < LandIntThreshold_value) .and. (sv%state(:,part_idx) == 1))
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
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> mixing length diffusion kernel, computes random velocities at given
    !> instants to model diffusion processes. These are valid while the tracer
    !> travels a given mixing length, propotional to the resolution of the
    !> background (and its ability to resove motion scales)
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function SullivanAllen(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), intent(in) :: dt
    real(prec), dimension(size(sv%state,1)) :: VelModH, TlagrangeH, RAND, HD, UD, VD
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: SullivanAllen
    integer :: np, part_idx, col_DifVelStdr, col_TPathHor
    real(prec) :: Pi = 4*atan(1.0)
    type(string) :: tag
    !Begin----------------------------------------------------------------------
    SullivanAllen = 0.0
    
    VelModH  = sqrt(sv%state(:,4)**2 + sv%state(:,5)**2)
                    
    tag = 'VelStandardDeviation'
    col_DifVelStdr = Utils%find_str(sv%varName, tag, .true.)
    
    tag = 'TPathHor'
    col_TPathHor = Utils%find_str(sv%varName, tag, .true.)
    
	sv%state(:,col_DifVelStdr) = sv%state(:,col_DifVelStdr) + Globals%Constants%VarVelHX * VelModH + Globals%Constants%VarVelH 

    where (sv%state(:,col_DifVelStdr) > 0.0)
        TlagrangeH = Globals%Constants%MixingLength / sv%state(:,col_DifVelStdr)
    elsewhere
        TlagrangeH = 0.0     
    endwhere

    !Copied over from MOHID LagrangianGlobal module
    
    call random_number(RAND)
    
    where (sv%state(:,col_TPathHor) >= TlagrangeH)
                        
        ! First step - compute the modulus of turbulent vector
                        
        !SQRT(3.0)=1.732050808 
        HD                       = 1.732050808 * sv%state(:,col_DifVelStdr) * RAND

        ! Second step - Compute the modulus of each component of the turbulent vector

        !   From 0 to Pi/2 cos and sin have positive values
        UD                       = HD * cos(2 * Pi * RAND)
        VD                       = HD * sin(2 * Pi * RAND)

        !Third step - Compute the direction of the the turbulent vector taking in consideration the layers thickness gradients
        ! Spagnol et al. (Mar. Ecol. Prog. Ser., 235, 299-302, 2002).
                        
        !if (CurrentOrigin%Movement%TurbGradK) then
        !    ![m/s] = [m/s] + [m] * [m/s] * [1/m] 
        !    ! K = Turbulent Diffusion Coefficent = MixingLength * Stdev of turbulent velocity / 2. 
        !    UD = UD + MixingLength * 1.732050808 * sv%state(:,col_DifVelStdr) / 2. * GradDWx
        !    VD = VD + MixingLength * 1.732050808 * sv%state(:,col_DifVelStdr) / 2. * GradDWy
        !
        !endif                                    
                        
        sv%state(:,col_TPathHor) = dt
        !Set dif velocities
        sv%state(:,7)     = UD
        sv%state(:,8)     = VD
        
        SullivanAllen(:,1) = Utils%m2geo(UD, sv%state(:,2), .false.)
        SullivanAllen(:,2) = Utils%m2geo(VD, sv%state(:,2), .true.)
                    
    elsewhere (sv%state(:,col_DifVelStdr) > 0.0)
        !Change only the position. dif velocities stay the same
        SullivanAllen(:,1) = Utils%m2geo(sv%state(:,7), sv%state(:,2), .false.)
        SullivanAllen(:,2) = Utils%m2geo(sv%state(:,8), sv%state(:,2), .true.)
        sv%state(:,col_TPathHor) = sv%state(:,col_TPathHor) + dt
    elsewhere
        !remove diffusion when there is no source of velocity
        sv%state(:,7) = 0.0
        sv%state(:,8) = 0.0
    end where
    
    end function SullivanAllen
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
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.09.12 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Lagrangian Kernel, evaluate the velocities at given points
    !> using the interpolants and split the evaluation part from the solver module.
	!> Also, for the particles near to the seabed uses log-law.
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    subroutine LagrangianVelModification(self, sv, bdata, time)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer ::  nf_w, nf_u, nf_v, col_u, col_dwz, col_v, col_w, part_idx, col_dist2bottom, col_temp, col_sal
    integer ::  col_rugosityVar_sv, col_D50Var_sv
    real(prec), dimension(:,:), allocatable :: var_dt, var_hor_dt
    type(string), dimension(:), allocatable :: var_name, var_name_hor
    type(string), dimension(:), allocatable :: requiredVars, requiredHorVars
    real(prec) :: VonKarman = 0.4
    real(prec) :: Hmin_Chezy = 0.1
    real(prec) :: A_Reichardt =	11.0
	real(prec) :: B_Reichardt =	0.33
    real(prec) :: C_Reichardt =	7.8
    real(prec), dimension(size(sv%state,1)) :: chezyZ, dist2bottom
    real(8), dimension(size(sv%state,1)) :: aux_r8, aux_r9, aux_r10
	real(prec), dimension(size(sv%state,1)) :: Temperatur_list, Salinity_list, depth_list
    real(prec), dimension(size(sv%state,1)) :: U_asterisk, V_asterisk, W_asterisk, z_plus_U_astr, z_plus_V_astr, z_plus_W_astr
    real(prec), dimension(size(sv%state,1)) :: u_Reichardt, v_Reichardt, w_Reichardt, u_VonKarman, v_VonKarman, w_VonKarman	
	real(prec), dimension(size(sv%state,1)) :: fDensity, kVisco, kViscoRelation
    real(prec) :: threshold_bot_wat, landIntThreshold
	real(prec), dimension(size(sv%state,1)) :: Threshold_value
	real(prec), dimension(size(sv%state,1)) :: LandIntThreshold_value
    type(string) :: tag
    integer :: i
    !-------------------------------------------------------------------------------------
    !write(*,*)"Entrada kinematic"
    tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)
	col_bat  = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
	col_temp = Utils%find_str(sv%varname, Globals%Var%temp, .false.)
	col_sal  = Utils%find_str(sv%varname, Globals%Var%sal, .false.)	
 	col_rugosityVar_sv = Utils%find_str(sv%varName, Globals%Var%rugosityVar, .true.)
 	col_D50Var_sv = Utils%find_str(sv%varName, Globals%Var%D50Var, .true.)
	
    allocate(requiredVars(3))
    requiredVars(1) = Globals%Var%u
    requiredVars(2) = Globals%Var%v
    requiredVars(3) = Globals%Var%w
    
    allocate(requiredHorVars(3))
    requiredHorVars(1) = Globals%Var%u
    requiredHorVars(2) = Globals%Var%v
    requiredHorVars(3) = Globals%Var%w
    !write(*,*)"Entrada interpolacao kinematic 1"
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
    !write(*,*)"Saida interpolacao kinematic 1"

    !Correct bottom values
    !write(*,*)"Entrada interpolacao kinematic 2"
    call KernelUtils%getInterpolatedFields(sv, bdata, time, requiredHorVars, var_hor_dt, var_name_hor, reqVertInt = .false.)
    !write(*,*)"Entrada interpolacao kinematic 2"
    
    col_u = Utils%find_str(var_name_hor, Globals%Var%u, .true.)
    col_v = Utils%find_str(var_name_hor, Globals%Var%v, .true.)
    col_w = Utils%find_str(var_name_hor, Globals%Var%w, .false.)
    col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    nf_u = Utils%find_str(var_name, Globals%Var%u, .true.)
    nf_v = Utils%find_str(var_name, Globals%Var%v, .true.)
    nf_w = Utils%find_str(var_name, Globals%Var%w, .false.)
	
! New version by Mohsen	
    ! threshold_bot_wat: equal to 0.0 will be correspond to the distance to the bottom equal to the last measurmet(given data).
    ! LandIntThresholdValue: will be correspond to the distance to the bottom which the volcity could be considered zero
	! LandIntThresholdValue gives the value the dist2bottom for this threshold.
	! u, v and w velocities for distance to the seabed(h(i) - bathymetry) < Rugosity reach to zero.    
 
    threshold_bot_wat = 0.0

	! Threshold_value: distance from the bottom (seabed) in unit [meter].
	! It could be a constant * Globals%Constants%Rugosity
	Threshold_value = Globals%Constants%BedLoadThickness 
	LandIntThreshold_value = VerticalMotion%LandIntThresholdValue(sv, bdata, time, Threshold_value)
	
    dist2bottom = sv%state(:,col_dist2bottom)

!	Calculate the kinematic Visco and density of seawater
    if ((col_temp /= MV_INT) .and. (col_sal /= MV_INT)) then
	! If there are salt and temperature data along the time:
	!	Calculate the kinematic Visco and density of seawater according to 
	!	the given temperature and salinity of seawater during the time. 
		        
		Temperatur_list = sv%state(:,col_temp)
		Salinity_list	= sv%state(:,col_sal)
		depth_list		= sv%state(:,3)
		fDensity = VerticalMotion%seaWaterDensity(Salinity_list, Temperatur_list, depth_list)
		kVisco = VerticalMotion%absoluteSeaWaterViscosity(Salinity_list, Temperatur_list) / fDensity

		kViscoRelation = abs(1.-(kvisco/Globals%Constants%MeanKvisco))
		where(kViscoRelation >= 0.9)
			kVisco = Globals%Constants%MeanKvisco
			fDensity = Globals%Constants%MeanDensity
		endwhere
    else
		! If there is no salt and temperature:
		! Consider the kinematic Visco and density of sea water according to 
		! the given constant values in Globals%Constants

		kVisco = Globals%Constants%MeanKvisco
		fDensity = Globals%Constants%MeanDensity
    end if

	! ...
	var_dt = merge(var_dt, 0.0, ieee_is_finite(var_dt) .and. var_dt >= -1000.0)
	var_hor_dt = merge(var_dt, 0.0, ieee_is_finite(var_hor_dt) .and. var_hor_dt >= -1000.0)


!	New version: Method 1: less common method and needs modification
!!	Reichardt law for the inner turbulent boundary layer(the particles with the depth more than the last measurmet towards the seabed)
!!	To find the U_asterik, it need a iterative Newton method
!!	I used log law to find it, but it is not a good approximation !!!
!    where (dist2bottom < threshold_bot_wat)
!		! if the sv%state(:,col_dwz) = Rugosity, then the ln(1) = 0.0 and X_asterisk = X / 0.0 ???? 
!		! Later, we should think about it, but it is not common. Beacuse, Rugosity < 0.005 meter!
!		aux_r8 = max(sv%state(:,col_dwz),Globals%Constants%Rugosity) / Globals%Constants%Rugosity
!		U_asterisk = var_hor_dt(:,col_u) * VonKarman /  dlog(aux_r8)
!		V_asterisk = var_hor_dt(:,col_v) * VonKarman /  dlog(aux_r8)
!		W_asterisk = var_hor_dt(:,col_w) * VonKarman /  dlog(aux_r8)
!		z_plus_U_astr = (sv%state(:,3) - sv%state(:,col_bat)) * abs(U_asterisk) / kVisco
!		z_plus_V_astr = (sv%state(:,3) - sv%state(:,col_bat)) * abs(V_asterisk) / kVisco
!		z_plus_W_astr = (sv%state(:,3) - sv%state(:,col_bat)) * abs(W_asterisk) / kVisco			
!		u_Reichardt = (U_asterisk) * ( 	(1.0 / VonKarman) * dlog(1.0 + VonKarman * z_plus_U_astr) + &
!										 C_Reichardt * (1.0 - dexp(-z_plus_U_astr/A_Reichardt) 	- &
!														(-z_plus_U_astr/A_Reichardt) * dexp(-z_plus_U_astr/B_Reichardt)) )
!
!		v_Reichardt = (V_asterisk) * ( 	(1.0 / VonKarman) * dlog(1.0 + VonKarman * z_plus_V_astr) + &
!										 C_Reichardt * (1.0 - dexp(-z_plus_V_astr/A_Reichardt) 	- &
!														(-z_plus_V_astr/A_Reichardt) * dexp(-z_plus_V_astr/B_Reichardt)) )
!
!		w_Reichardt = (U_asterisk) * ( 	(1.0 / VonKarman) * dlog(1.0 + VonKarman * z_plus_W_astr) + &
!										 C_Reichardt * (1.0 - dexp(-z_plus_W_astr/A_Reichardt) 	- &
!														(-z_plus_W_astr/A_Reichardt) * dexp(-z_plus_W_astr/B_Reichardt)) )
!
!		sv%state(:,4) = u_Reichardt
!		sv%state(:,5) = v_Reichardt
!		sv%state(:,6) = w_Reichardt
!
!	else where
!	
!        sv%state(:,4) = var_dt(:,nf_u)
!        sv%state(:,5) = var_dt(:,nf_v)	
!		sv%state(:,6) = var_dt(:,nf_w)	
!		
!    end where	

!	New version: Method 2: more common and used method
!	Log law for the inner turbulent boundary layer(the particles with the depth more than the last measurmet towards the seabed)

	where (dist2bottom < threshold_bot_wat)

		aux_r8 = max(sv%state(:,col_dwz), 0.05) / sv%state(:,col_rugosityVar_sv)
		U_asterisk = var_hor_dt(:,col_u) * VonKarman / dlog(aux_r8)
		V_asterisk = var_hor_dt(:,col_v) * VonKarman / dlog(aux_r8)
		W_asterisk = var_hor_dt(:,col_w) * VonKarman / dlog(aux_r8)

		where (dist2bottom > LandIntThreshold_value)    !LandIntThreshold_value corresponded to bed-load/land thickness
			!--- Case 1: water column above bed-load/land interaction
			aux_r9 = max(sv%state(:,3) - sv%state(:,col_bat), sv%state(:,col_rugosityVar_sv)) / sv%state(:,col_rugosityVar_sv)
			sv%state(:,4) = (U_asterisk / VonKarman) * dlog(aux_r9)
			sv%state(:,5) = (V_asterisk / VonKarman) * dlog(aux_r9)
			sv%state(:,6) = (W_asterisk / VonKarman) * dlog(aux_r9)
			
			!z_plus_U_astr = (sv%state(:,3) - sv%state(:,col_bat)) * abs(U_asterisk) / kVisco
			!z_plus_V_astr = (sv%state(:,3) - sv%state(:,col_bat)) * abs(V_asterisk) / kVisco
			!z_plus_W_astr = (sv%state(:,3) - sv%state(:,col_bat)) * abs(W_asterisk) / kVisco			

		elsewhere
			!--- Case 2: inside bed-load/land interaction layer: use an average velocity of log-law for all particles in this layer
			aux_r10 = Globals%Constants%BedLoadThickness / sv%state(:,col_rugosityVar_sv)
			sv%state(:,4) = (U_asterisk / VonKarman) * (dlog(aux_r10) - 1.0 + 1.0/aux_r10)
			sv%state(:,5) = (V_asterisk / VonKarman) * (dlog(aux_r10) - 1.0 + 1.0/aux_r10)
			sv%state(:,6) = (W_asterisk / VonKarman) * (dlog(aux_r10) - 1.0 + 1.0/aux_r10)
		end where

	else where
		!--- Case 3: elsewhere, keep unmodified velocities
		sv%state(:,4) = var_dt(:,nf_u)
		sv%state(:,5) = var_dt(:,nf_v)
		sv%state(:,6) = var_dt(:,nf_w)
	end where
	
	deallocate(var_dt)
    deallocate(var_hor_dt)
    deallocate(var_name_hor)
    deallocate(var_name)
	
	end subroutine LagrangianVelModification

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
    !call MOHIDWaterQuality%initialize()
    
    call KernelUtils%initialize() 
    
    end subroutine initKernel

    end module kernel_mod