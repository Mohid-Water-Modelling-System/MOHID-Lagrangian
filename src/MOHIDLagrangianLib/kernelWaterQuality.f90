    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernelMOHIDWaterQuality
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for MOHIDWaterQuality associated processes.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach.
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
    !------------------------------------------------------------------------------

    module kernelMOHIDWaterQuality_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod
    use kernelUtils_mod
    use ModuleWaterQuality


    type :: kernelMOHIDWaterQuality_class        !< MOHIDWaterQuality kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelMOHIDWaterQuality
    procedure :: WQProcess
    procedure :: Dilution
    end type kernelMOHIDWaterQuality_class
    
    type(kernelUtils_class) :: KernelUtils_MOHIDWaterQuality   !< kernel utils

    public :: kernelMOHIDWaterQuality_class
    contains
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> changes concentration of water quality parameters using the MOHIDWaterQuality module
    !> @param[in] self, bdata, sv, time, dt
    !---------------------------------------------------------------------------
    function WQProcess(self, sv, bdata, time, dt)
    class(kernelMOHIDWaterQuality_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), intent(in) :: dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: WQProcess
    real(prec), dimension(13,size(sv%state,1)), target :: mass
    real(prec), dimension(:,:), pointer :: mass_pointer
    logical, dimension (size(sv%state,1)) :: computeFlag
    integer :: c_O2, c_NO3, c_NO2, c_PO4
    integer :: c_DON_NR, c_DOP_NR, c_DON_R, c_DOP_R
    integer :: c_Pytho, c_Zoo, c_Age, c_Temp, c_Sal, c_SWper, c_LextCoef, c_NH4, c_Rad
    integer :: c_PartOrgNit, c_PartOrgPho
    integer :: zoo, phyto, ammonia, nitrate, nitrite, dONRefractory, dONNonRefractory, partOrganicNitrogen
    integer :: dOPRefractory, dOPNonRefractory, partOrganicPhosphorus, inorganicPhosphorus, oxygen
    real(prec), dimension(2) :: maxLevel
    real(prec), dimension(size(sv%state,1)) :: depth
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    type(string) :: tag, outext
    integer :: MOHIDWaterQualityID = 0
    integer :: STAT_CALL = 0
    !begin---------------------------------------------------------------------------
    WQProcess = 0.0
    mass = 0.0
    c_Temp = Utils%find_str(sv%varName, Globals%Var%temp, .true.)
    c_Sal = Utils%find_str(sv%varName, Globals%Var%sal, .true.)
    tag = 'age'
    c_Age = Utils%find_str(sv%varName, tag, .true.)
    tag = 'dissolved_oxygen'
    c_O2 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'ammonia'
    c_NH4 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'nitrate'
    c_NO3 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'nitrite'
    c_NO2 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'inorganic_phosphorus'
    c_PO4 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DON_NonRefractory'
    c_DON_NR = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DOP_NonRefractory'
    c_DOP_NR = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DON_Refractory'
    c_DON_R = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DOP_Refractory'
    c_DOP_R = Utils%find_str(sv%varName, tag, .true.)
    tag = 'phytoplankton'
    c_Pytho = Utils%find_str(sv%varName, tag, .true.)
    tag = 'zooplankton'
    c_Zoo = Utils%find_str(sv%varName, tag, .true.)
    tag = 'partOrgNit'
    c_PartOrgNit = Utils%find_str(sv%varName, tag, .true.)
    tag = 'partOrgPho'
    c_PartOrgPho = Utils%find_str(sv%varName, tag, .true.)
    tag = 'sw_extinction_coef'
    c_LextCoef = Utils%find_str(sv%varName, tag, .true.)
    tag = 'sw_percentage'
    c_SWper = Utils%find_str(sv%varName, tag, .true.)
    
    !Creates flag to identify which tracers concentration will be updated after computation
    where (mod(sv%state(:, c_Age), Globals%SimDefs%WqDt) > 1E-3)
        computeFlag = .true.
    elsewhere
        computeFlag = .false.
    endwhere
    
    if (any(computeFlag)) then
        !Fills the 2D mass array WaterQuality module uses.
        !Gets the indexes from the WQM - for building mass matrix
        call GetWQPropIndex(0, Zoo = zoo, Phyto = phyto, Ammonia = ammonia, Nitrate = nitrate,                  &
                            Nitrite = nitrite, DissOrganicNitrogenRefractory = dONRefractory,                   &
                            DONNonRefractory = dONNonRefractory, PartOrganicNitrogen = partOrganicNitrogen,     &
                            DissOrganicPhosphorusRefractory = dOPRefractory,                                    &
                            DOPNonRefractory = dOPNonRefractory, PartOrganicPhosphorus = partOrganicPhosphorus, &
                            InorganicPhosphorus = inorganicPhosphorus, Oxygen = oxygen, STAT = STAT_CALL)
        if (STAT_CALL /= 0) then
            outext='Failed to get water property indexes from MOHID module water quality. stopping'
            call Log%put(outext)
            stop
        endif
        
        allocate(requiredVars(1))
        requiredVars(1) = Globals%Var%rad
        
        call KernelUtils_MOHIDWaterQuality%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
        !surface radiation
        c_rad = Utils%find_str(var_name, Globals%Var%rad, .true.)

        !Get radiation at vertical center of tracers
        var_dt(:, c_rad) = var_dt(:, c_rad) * sv%state(:,c_SWper)
        
        depth = sv%state(:,3)
        maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .false.)
        !compute distance to surface
        depth = maxLevel(2) - depth
        
        !everywhere else will be 0
        where (computeFlag == .true.)
            mass(zoo, :) = sv%state(:,c_Zoo)
            mass(phyto, :) = sv%state(:,c_Pytho)
            mass(ammonia, :) = sv%state(:,c_NH4)
            mass(nitrate, :) = sv%state(:,c_NO3)
            mass(nitrite, :) = sv%state(:,c_NO2)
            mass(dONRefractory, :) = sv%state(:,c_DON_R)
            mass(dONNonRefractory, :) = sv%state(:,c_DON_NR)
            mass(partOrganicNitrogen, :) = sv%state(:,c_PartOrgNit)
            mass(dOPRefractory, :) = sv%state(:,c_DOP_R)
            mass(dOPNonRefractory, :) = sv%state(:,c_DOP_NR)
            mass(partOrganicPhosphorus, :) = sv%state(:,c_PartOrgPho)
            mass(inorganicPhosphorus, :) = sv%state(:,c_PO4)
            mass(oxygen, :) = sv%state(:,c_O2)
        endwhere
        
        !Get a pointer to use in mohid water
        mass_pointer => mass
        !main call to water quality processes
        call WaterQuality(0, GetPointer(sv%state(:,c_Sal)), GetPointer(sv%state(:,c_Temp)),     &
                            GetPointer(var_dt(:, c_Rad)), GetPointer(sv%state(:,c_LextCoef)),   &
                            GetPointer(depth),                                                  &
                            mass_pointer,                                                       &
                            1, size(sv%state,1),                                                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= 0) then
            outext='Failed to iterate through module WaterQuality. stopping'
            call Log%put(outext)
            stop
        endif
    
        !Unpack new values (check for negative values?)
        !Only changes tracers that were computed
        where (computeFlag)
            WQProcess(:,c_Zoo) = (mass_pointer(Zoo, :) - sv%state(:,c_Zoo)) / dt
            WQProcess(:,c_Pytho) = (mass_pointer(Phyto, :) - sv%state(:,c_Phyto)) / dt
            WQProcess(:,c_NH4) = (mass_pointer(Ammonia, :) - sv%state(:,c_NH4)) / dt
            WQProcess(:,c_NO3) = (mass_pointer(Nitrate, :) - sv%state(:,c_NO3)) / dt
            WQProcess(:,c_NO2) = (mass_pointer(Nitrite, :) - sv%state(:,c_NO2)) / dt
            WQProcess(:,c_DON_R) = (mass_pointer(DONRefractory, :) - sv%state(:,c_DON_R)) / dt
            WQProcess(:,c_DON_NR) = (mass_pointer(DONNonRefractory, :) - sv%state(:,c_DON_NR)) / dt
            WQProcess(:,c_PartOrgNit) = (mass_pointer(PartOrganicNitrogen, :) - sv%state(:,c_PartOrgNit)) / dt
            WQProcess(:,c_DOP_R) = (mass_pointer(DOPRefractory, :) - sv%state(:,c_DOP_R)) / dt
            WQProcess(:,c_DOP_NR) = (mass_pointer(DOPNonRefractory, :) - sv%state(:,c_DOP_NR)) / dt
            WQProcess(:,c_PartOrgPho) = (mass_pointer(PartOrganicPhosphorus, :) - sv%state(:,c_PartOrgPho)) / dt
            WQProcess(:,c_PO4) = (mass_pointer(InorganicPhosphorus, :) - sv%state(:,c_PO4)) / dt
            WQProcess(:,c_O2) = (mass_pointer(Oxygen, :) - sv%state(:,c_O2)) / dt
        endwhere
        
    endif
    
    end function WQProcess
    
    !> @author Joao Sobrinho
    !> @brief
    !> Gets a pointer to a matrix, used in mohid waterquality
    !> @param[in] ptr
    !---------------------------------------------------------------------------
    function GetPointer(matrix) result(ptr)
    !Arguments-------------------------------------------------------------
    real(prec), dimension(:), intent(in), target :: matrix

    !Local-----------------------------------------------------------------
    real(prec), dimension(:), pointer             :: ptr

    ptr => matrix
    end function GetPointer

    !> @author Joao Sobrinho - +Atlantic
    !> @brief
    !> Computes the dilution of dissolved material in the water column by increasing its volume
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Dilution(self, sv, bdata, time, dt)
    class(kernelMOHIDWaterQuality_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: dt, time
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Dilution
    real(prec), dimension(size(sv%state,1)) :: velocity_mod, volume_new, volumeRatio, double_vol_time
    integer :: c_Volume, c_O2, c_NO3, c_NO2, c_PO4, c_DON_NR, c_DOP_NR, c_DON_R, c_DOP_R
    integer :: c_Pytho, c_Zoo, c_Age, c_NH4, c_PartOrgNit, c_PartOrgPho
    integer :: vol_Idx, nf_u, nf_v, nf_O2, nf_NH4, nf_NO3, nf_NO2, nf_PO4, nf_DON_NR, nf_DOP_NR, nf_DON_R, nf_DOP_R
    integer :: nf_Pytho, nf_Zoo, nf_PartOrgNit, nf_PartOrgPho
    type(string) :: tag
    
    !TODO : add temperature and salinity, and refactor routine to be added in kernel utils. This way dilution can be called
    ! for any kernel that requires it
    !TODO : create arrays for all these variables and create  routines to get an array of column IDs.
    Dilution = 0.0
    
    allocate(requiredVars(15))
    requiredVars(1) = Globals%Var%u
    requiredVars(2) = Globals%Var%v
    requiredVars(3) = Globals%Var%dissolved_oxygen
    requiredVars(4) = Globals%Var%nitrate
    requiredVars(5) = Globals%Var%nitrite
    requiredVars(6) = Globals%Var%ammonia
    requiredVars(7) = Globals%Var%partOrgNit
    requiredVars(8) = Globals%Var%partOrgPho
    requiredVars(9) = Globals%Var%DON_NonRefractory
    requiredVars(10) = Globals%Var%DOP_NonRefractory
    requiredVars(11) = Globals%Var%DON_Refractory
    requiredVars(12) = Globals%Var%DOP_Refractory
    requiredVars(13) = Globals%Var%phytoplankton
    requiredVars(14) = Globals%Var%zooplankton
    requiredVars(15) = Globals%Var%inorganic_phosphorus
    
    tag = 'age'
    c_Age = Utils%find_str(sv%varName, tag, .true.)
    tag = 'volume'
    vol_Idx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'dissolved_oxygen'
    c_O2 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'ammonia'
    c_NH4 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'nitrate'
    c_NO3 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'nitrite'
    c_NO2 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'inorganic_phosphorus'
    c_PO4 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DON_NonRefractory'
    c_DON_NR = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DOP_NonRefractory'
    c_DOP_NR = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DON_Refractory'
    c_DON_R = Utils%find_str(sv%varName, tag, .true.)
    tag = 'DOP_Refractory'
    c_DOP_R = Utils%find_str(sv%varName, tag, .true.)
    tag = 'phytoplankton'
    c_Pytho = Utils%find_str(sv%varName, tag, .true.)
    tag = 'zooplankton'
    c_Zoo = Utils%find_str(sv%varName, tag, .true.)
    tag = 'partOrgNit'
    c_PartOrgNit = Utils%find_str(sv%varName, tag, .true.)
    tag = 'partOrgPho'
    c_PartOrgPho = Utils%find_str(sv%varName, tag, .true.)
    
    call KernelUtils_MOHIDWaterQuality%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
    
    nf_u = Utils%find_str(var_name, Globals%Var%u, .true.)
    nf_v = Utils%find_str(var_name, Globals%Var%v, .true.)
    nf_O2 = Utils%find_str(var_name, Globals%Var%dissolved_oxygen, .true.)
    nf_NH4 = Utils%find_str(var_name, Globals%Var%ammonia, .true.)
    nf_NO3 = Utils%find_str(var_name, Globals%Var%nitrate, .true.)
    nf_NO2 = Utils%find_str(var_name, Globals%Var%nitrite, .true.)
    nf_PO4 = Utils%find_str(var_name, Globals%Var%inorganic_phosphorus, .true.)
    nf_DON_NR = Utils%find_str(var_name, Globals%Var%DON_NonRefractory, .true.)
    nf_DOP_NR = Utils%find_str(var_name, Globals%Var%DOP_NonRefractory, .true.)
    nf_DON_R = Utils%find_str(var_name, Globals%Var%DON_Refractory, .true.)
    nf_DOP_R = Utils%find_str(var_name, Globals%Var%DOP_Refractory, .true.)
    nf_Pytho = Utils%find_str(var_name, Globals%Var%phytoplankton, .true.)
    nf_Zoo = Utils%find_str(var_name, Globals%Var%zooplankton, .true.)
    nf_PartOrgNit = Utils%find_str(var_name, Globals%Var%partOrgNit, .true.)
    nf_PartOrgPho = Utils%find_str(var_name, Globals%Var%partOrgPho, .true.)
    
    ! compute velocity modulus.
    velocity_mod = sqrt(var_dt(:,nf_u)**2 + var_dt(:,nf_v)**2)
                
    !Time to double the volume (TVol200 in old lagrangian). Linear function with duplication in 30min for 2m/s and 2h for 0.1m/s.
    double_vol_time = -2840.0 * velocity_mod + 7485.0
    !Update volume
    volume_new = sv%state(:,vol_Idx) * exp((alog(2.)/double_vol_time) * dt)
                
    !Dilution in volume = Volume variation / dt
    Dilution(:, vol_Idx) = (volume_new - sv%state(:, vol_Idx)) / dt
    
    !Reducing concentration
    volumeRatio = (volume_new / sv%state(:, vol_Idx)) - 1.0
    
    !TODO : make a cicle for each variable (returned after implementing TODO in the begining of the routine)
    Dilution(:, c_O2) = - ((sv%state(:, c_O2) - var_dt(:, c_O2)) * volumeRatio) / dt
    Dilution(:, c_NH4) = - ((sv%state(:, c_NH4) - var_dt(:, c_NH4)) * volumeRatio) / dt
    Dilution(:, c_NO3) = - ((sv%state(:, c_NO3) - var_dt(:, c_NO3)) * volumeRatio) / dt
    Dilution(:, c_NO2) = - ((sv%state(:, c_NO2) - var_dt(:, c_NO2)) * volumeRatio) / dt
    Dilution(:, c_PO4) = - ((sv%state(:, c_PO4) - var_dt(:, c_PO4)) * volumeRatio) / dt
    Dilution(:, c_DON_NR) = - ((sv%state(:, c_DON_NR) - var_dt(:, c_DON_NR)) * volumeRatio) / dt
    Dilution(:, c_DOP_NR) = - ((sv%state(:, c_DOP_NR) - var_dt(:, c_DOP_NR)) * volumeRatio) / dt
    Dilution(:, c_DON_R) = - ((sv%state(:, c_DON_R) - var_dt(:, c_DON_R)) * volumeRatio) / dt
    Dilution(:, c_DOP_R) = - ((sv%state(:, c_DOP_R) - var_dt(:, c_DOP_R)) * volumeRatio) / dt
    Dilution(:, c_Pytho) = - ((sv%state(:, c_Pytho) - var_dt(:, c_Pytho)) * volumeRatio) / dt
    Dilution(:, c_Zoo) = - ((sv%state(:, c_Zoo) - var_dt(:, c_Zoo)) * volumeRatio) / dt
    Dilution(:, c_PartOrgPho) = - ((sv%state(:, c_PartOrgPho) - var_dt(:, c_PartOrgPho)) * volumeRatio) / dt
    Dilution(:, c_PartOrgPho) = - ((sv%state(:, c_PartOrgPho) - var_dt(:, c_PartOrgPho)) * volumeRatio) / dt
    
    ! Eliminate tracers older than 5 days (tracers will be too large)
    !TODO : implement an initial volume variable and use it to define when tracers are to be deleted. Copy from Detritus
    where(sv%state(:,c_Age) > 43200) sv%active = .false.
    
    end function Dilution
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelMOHIDWaterQuality(self)
    class(kernelMOHIDWaterQuality_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initKernelMOHIDWaterQuality

    end module kernelMOHIDWaterQuality_mod