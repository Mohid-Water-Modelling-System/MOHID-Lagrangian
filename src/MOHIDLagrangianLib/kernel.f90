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

    use kernelSharedProcesses_mod    
    use kernelLitter_mod
    use kernelVerticalMotion_mod
    use kernelColiform_mod
    use kernelDetritus_mod
    ! use kernelMOHIDWaterQuality_mod
	use kernelFloater_mod

    type :: kernel_class        !< Kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
		type(kernelSharedProcesses_class) :: SharedProcesses   !< Shared Processes kernels
		type(kernelLitter_class) :: Litter       !< litter kernels
		type(kernelVerticalMotion_class) :: VerticalMotion   !< VerticalMotion kernels
		type(kernelColiform_class) :: Coliform !< coliform kernels
		type(kernelDetritus_class) :: Detritus !< coliform kernels
		! type(kernelMOHIDWaterQuality_class) :: MOHIDWaterQuality   !< VerticalMotion kernels
		type(kernelFloater_class) :: Floater   !< kernel Floater
		type(kernelUtils_class) :: KernelUtils   !< kernel utils
    contains
    procedure :: initialize => initKernel
    procedure :: run => runKernel
    procedure, private :: setCommonProcesses
    procedure, private :: interpolate_backgrounds
    procedure, private :: distance2bottom															
    end type kernel_class



	
 
    public :: kernel_class
    contains

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
	!> Modified @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Adaptation from runSolver (Ricardo) method that evaluates the specific
    !> kernel, according to the selected kernel
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function runKernel(self, sv, bdata, time, dt)
    class(kernel_class), intent(inout) :: self
    type(stateVector_class), intent(inout)	:: sv
    type(stateVector_class) 				:: sv_atDepth
	type(source_parameters_type)			:: floaterSrcPar
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: runKernel
    integer :: i, i_particle
	integer :: col_bathy	
	logical, save :: printed = .false.
    !write(*,*)"Entrada Run Kernel tamanho bdata =", size(bdata)
    !do i = 1, size(bdata)
    !    write(*,*)"Tamanho background i = ", i, bdata(i)%fields%getSize()
    !enddo
    !running preparations for kernel lanch
	
	logical, allocatable :: mask_atDepth(:)
    
    call self%setCommonProcesses(sv, bdata, time)
    !write(*,*)"Entrada interpolate_backgrounds"
    call self%interpolate_backgrounds(sv, bdata, time)
    !Computes distance to bottom for all tracers
    !write(*,*)"Entrada distance2bottom"
    call self%distance2bottom(sv)
    !running kernels for each type of tracer
    !write(*,*)"Entrada kernels"
	
	! Modify velocities based on Reichardt and log law for the inner turbulent layer near to the seabed.
	call self%SharedProcesses%LagrangianVelModification(sv, bdata, time)
	
    if (sv%ttype == Globals%Types%base) then
        runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    self%SharedProcesses%Windage(sv, bdata, time) + & 
					self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%SharedProcesses%Aging(sv)
    else if (sv%ttype == Globals%Types%paper) then
        runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    self%SharedProcesses%Windage(sv, bdata, time) + self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%SharedProcesses%Aging(sv) + self%Litter%DegradationLinear(sv) + self%VerticalMotion%Buoyancy(sv, bdata, time) + &
                    self%VerticalMotion%Resuspension(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%plastic) then
        runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    self%SharedProcesses%Windage(sv, bdata, time) + self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%SharedProcesses%Aging(sv) + self%Litter%DegradationLinear(sv) + self%Litter%BioFouling(sv, dt) + self%VerticalMotion%Buoyancy(sv, bdata, time) + &
                    self%VerticalMotion%Resuspension(sv, bdata, time, dt)
		
	else if (sv%ttype == Globals%Types%floater) then
		call self%Floater%findSourceParameters(sv, floaterSrcPar)
		call sv%copyState(sv_atDepth)
		
		!sv_atDepth%state(:,3) = sv%state(:,3)
		
		col_bathy = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
		sv_atDepth%state(:,3) = max(sv%state(:,3) - floaterSrcPar%P_D(:), min(sv%state(:,col_bathy) + 1.0, sv%state(:,3)))

		! sv_atDepth%state(:,3) = min( max(sv%state(:,3) - floaterSrcPar%P_D(:), sv%state(:,col_bathy) + 0.1_prec), 0.0_prec )

		! col_bathy = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
		! where (sv%state(:,3) - floaterSrcPar%P_D(:) > sv%state(:,col_bathy) + 0.01 ) 
			! sv_atDepth%state(:,3) = sv%state(:,3) - floaterSrcPar%P_D(:)
		! elsewhere
			! sv_atDepth%state(:,3) = sv%state(:,col_bathy) + 0.01
		! endwhere
		
		call self%setCommonProcesses(sv_atDepth, bdata, time)
		call self%interpolate_backgrounds(sv_atDepth, bdata, time) 
		call self%distance2bottom(sv_atDepth)	
		call self%SharedProcesses%LagrangianVelModification(sv_atDepth, bdata, time) 		
		runKernel = self%Floater%FloaterVelocity(sv, sv_atDepth, bdata, time, dt, floaterSrcPar)  + self%SharedProcesses%Aging(sv)
		call sv_atDepth%finalize()

	! else if (sv%ttype == Globals%Types%floater) then

		! call self%Floater%findSourceParameters(sv, floaterSrcPar)
		! call sv%copyState(sv_atDepth)

		! col_bathy = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)

		! allocate(mask_atDepth(size(sv%state,1)))

		! ! True where the requested depth is above the bathymetry
		! mask_atDepth(:) = sv%state(:,3) - floaterSrcPar%P_D(:) > sv%state(:,col_bathy)

		! ! Move only valid particles to the requested depth
		! where (mask_atDepth)
			! sv_atDepth%state(:,3) = sv%state(:,3) - floaterSrcPar%P_D(:)
		! end where

		! call self%setCommonProcesses(sv_atDepth, bdata, time)
		! call self%interpolate_backgrounds(sv_atDepth, bdata, time)
		! call self%distance2bottom(sv_atDepth)
		! call self%SharedProcesses%LagrangianVelModification(sv_atDepth, bdata, time)

		! ! Compute floater velocity
		! runKernel = self%Floater%FloaterVelocity(sv, sv_atDepth, bdata, time, dt, floaterSrcPar)

		! ! Set velocity to zero for invalid particles
		! where (spread(.not. mask_atDepth, dim=2, ncopies=size(runKernel,2)))
			! runKernel = 0.0_prec
		! end where

		! deallocate(mask_atDepth)

		! call sv_atDepth%finalize()	
			
		
    else if (sv%ttype == Globals%Types%coliform) then
        runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + &
                    self%SharedProcesses%Aging(sv) + self%Coliform%MortalityT90(sv, bdata, time) + self%Coliform%Dilution(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%seed) then
        runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + self%SharedProcesses%Aging(sv) + &
                    self%VerticalMotion%Buoyancy(sv, bdata, time) + self%VerticalMotion%Resuspension(sv, bdata, time, dt)
    else if (sv%ttype == Globals%Types%detritus) then
        runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + self%SharedProcesses%Aging(sv) + &
                    self%VerticalMotion%Buoyancy(sv, bdata, time) + self%VerticalMotion%Resuspension(sv, bdata, time, dt) + &
                    self%Detritus%Degradation(sv, dt)
    ! else if (sv%ttype == Globals%Types%WaterQuality) then
        ! runKernel = self%SharedProcesses%LagrangianKinematic(sv, bdata, time) + self%SharedProcesses%StokesDrift(sv, bdata, time) + &
                    ! self%SharedProcesses%Windage(sv, bdata, time) + self%SharedProcesses%DiffusionMixingLength(sv, bdata, time, dt) + &
                    ! self%SharedProcesses%Aging(sv) + self%MOHIDWaterQuality%WQProcess(sv, bdata, time, dt) + self%MOHIDWaterQuality%Dilution(sv, bdata, time, dt)
    end if
    if (Globals%simDefs%FreeLitterAtBeaching == 1) then
        runKernel = self%SharedProcesses%FreeLitterAtBeaching(sv, bdata, time, runKernel, dt)
    else
        runKernel = self%SharedProcesses%Beaching(sv, bdata, time, runKernel)
    endif
    
    runKernel = self%VerticalMotion%CorrectVerticalBounds(sv, runKernel, bdata, time, dt)
    	
	
    end function runKernel

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com
	!> Modified @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
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
    integer :: i, j, col_age, col_bat, col_bat_sv, col_landintmask, col_ssh, col_res, col_DifVelStdr
	integer :: col_rugosityVar,col_rugosityVar_sv
	integer :: col_D50Var,col_D50Var_sv
	integer :: counterr
	real(prec), dimension(2) :: maxLevel
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
    call self%KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true.)

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
    
	if (col_rugosityVar /= MV_INT) then
		sv%state(:,col_rugosityVar_sv) = max(var_dt(:,col_rugosityVar), 0.0001_prec)
	else
		sv%state(:,col_rugosityVar_sv) = 0.0001_prec
	end if

    !set tracer bottom D50
    col_D50Var = Utils%find_str(var_name, Globals%Var%D50Var, .false.)
    col_D50Var_sv = Utils%find_str(sv%varName, Globals%Var%D50Var, .true.)
	if (col_D50Var /= MV_INT) then
		sv%state(:,col_D50Var_sv) = max(var_dt(:,col_D50Var), 0.0001_prec)
	else
		sv%state(:,col_D50Var_sv) = 0.0001_prec
	end if
	
!	counterr = 0
!	do i= 1, size(sv%state,1)
!		if (mod(counterr, 10) == 0) then
!			write(*,'(  A5, A12, A12)') , " Id:", 'rugosityVar', 'D50Var'
!			write(*,*),' '
!		end if
!		counterr = counterr + 1
!		write(*,'( I5, F12.4, F12.4)') , i, sv%state(i,col_rugosityVar_sv), sv%state(i,col_D50Var_sv)
!	end do	 	

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

!	Mohsen: I think it is not needed. A change has been done in function CorrectVerticalBounds to consider this modification.
!    !TODO : If the hdf5 does not have ssh, should use the verticalZ (meaning we must save an extra 2D var with the hdf original var (cell faces)
!    if (Globals%simDefs%inputFromHDF5) then
!        col_ssh = Utils%find_str(var_name, Globals%Var%ssh, .false.)
!        if (col_ssh /= MV_INT) then
!            where (sv%state(:,3) >  var_dt(:,col_ssh)) sv%state(:,3) = var_dt(:,col_ssh) - 0.00001
!        else
!            !ssh not found... assume 0.0 as the limit
!            where (sv%state(:,3) >  0.0) sv%state(:,3) = - 0.00001
!        endif
!        
!    else
!        maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .false.)   
!        if (maxLevel(2) /= MV) where (sv%state(:,3) > maxLevel(2)) sv%state(:,3) = maxLevel(2)-0.00001  
!    endif

    !update land interaction status
    ! col_landintmask = Utils%find_str(var_name, Globals%Var%landIntMask)
    ! sv%landIntMask = var_dt(:,col_landintmask)
	col_landintmask = Utils%find_str(var_name, Globals%Var%landIntMask, .false.)
	if (col_landintmask /= MV_INT) then
		sv%landIntMask = var_dt(:,col_landintmask)
	else
		print *,"sv%landIntMask is failed"
		sv%landIntMask = 0.0_prec
	end if
	
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
	!> Modified @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
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
	real(prec), dimension(2) :: maxLevel
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

    call self%KernelUtils%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, justRequired = .true.)
    
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
	!> Modified @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
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
    !> @author Daniel Garaboa Paz - GFNL
	!> Modified @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernel(self)
    class(kernel_class), intent(inout) :: self
    type(string) :: interpName
    
    interpName = 'linear'
	
    call self%Interpolator%initialize(1,interpName)
    call self%KernelUtils%initialize() 
    call self%SharedProcesses%initialize() 
    call self%Litter%initialize()
    call self%VerticalMotion%initialize()
    call self%Coliform%initialize()
    call self%Detritus%initialize()
    ! call self%MOHIDWaterQuality%initialize()   
    call self%Floater%initialize()
	
    end subroutine initKernel

    end module kernel_mod