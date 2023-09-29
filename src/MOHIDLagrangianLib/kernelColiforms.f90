    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernelColiform
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for coliform associated processes.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach.
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
    !------------------------------------------------------------------------------

    module kernelColiform_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod
    use kernelUtils_mod

    type :: kernelColiform_class        !< Coliform kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelColiform
    procedure :: MortalityT90
    procedure :: Dilution

    end type kernelColiform_class
    
    type(kernelUtils_class) :: KernelUtils_coliform   !< kernel utils

    public :: kernelColiform_class
    contains
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - +Atlantic
    !> @brief
    !> T_90 fecal coliforms decay kernel.
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function MortalityT90(self, sv, bdata, time)
    class(kernelColiform_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: MortalityT90
    real(prec), dimension(size(sv%state,1)) :: T90, Radiation
    type(string), dimension(:), allocatable :: requiredVars
    integer :: col_SWcoef, col_SWper, T90_method, col_conc, col_T90, col_T90_var, col_T90_varM
    integer :: col_rad, col_sal, col_temp
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string) :: tag
    !Begin------------------------------------------------------------------------
    MortalityT90 = 0.0
    tag = 'T90'
    col_T90 = Utils%find_str(sv%varName, tag, .true.)
    tag = 'T90_variable'
    col_T90_var = Utils%find_str(sv%varName, tag, .true.)
    tag = 'T90_method'
    col_T90_varM = Utils%find_str(sv%varName, tag, .true.)
    tag = 'sw_extinction_coef'
    col_SWcoef = Utils%find_str(sv%varName, tag, .true.)
    tag = 'sw_percentage'
    col_SWper = Utils%find_str(sv%varName, tag, .true.)
        
    if (all(sv%state(:,col_T90_varM) == 1)) then
        T90_method = 1
    else
        T90_method = 2
    end if
    if (all(sv%state(:,col_T90_var) == 1)) then
        allocate(requiredVars(1))
        requiredVars(1) = Globals%Var%rad
        
        call KernelUtils_coliform%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
        
        col_temp = Utils%find_str(sv%varname, Globals%Var%temp, .true.)
        col_sal = Utils%find_str(sv%varname, Globals%Var%sal, .true.)
        !surface radiation
        col_rad = Utils%find_str(var_name, Globals%Var%rad, .true.)
        
        !Computes radiation at vertical center of tracers
        call KernelUtils_coliform%getSWRadiation(sv, var_dt, bdata, col_rad, col_SWper, col_SWcoef, Radiation)
        
        !Compute T90
        if (T90_method == 1) then
            !Canteras
                
            T90 = 2.533 * (1.04**(sv%state(:,col_temp) - 20.)) * (1.012**sv%state(:,col_sal)) + (0.113 * Radiation)
            !T90 = 2.533 * (1.04**(var_dt(:,col_temp) - 20.)) * (1.012**var_dt(:,col_sal)) + (0.113 * Radiation)
                        
            ![]                        
            T90 = (2.303 / T90) * 24 * 3600                        
                        
            MortalityT90(:,col_conc)  =  - sv%state(:,col_conc) * (alog(10.0) / T90)
            !
            
        else if (T90_method == 2) then
            !Chapra not yet implemented
        end if
                    
        deallocate(var_name)
        deallocate(var_dt)
    else
        !Use constant T90
        MortalityT90(:,col_conc) = - (sv%state(:,col_conc) * sv%state(:,col_T90))
    endif
        
    !matar coliformes com menos de 100 UFC/100ml
    where(sv%state(:,col_conc) < 100) sv%active = .false.
    
    end function MortalityT90
    
    
    !> @author Joao Sobrinho - +Atlantic
    !> @brief
    !> Computes the dilution of dissolved material in the water column by increasing its volume
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Dilution(self, sv, bdata, time, dt)
    class(kernelColiform_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    integer :: np, nf, bkg, init_vol_IDx, vol_Idx, conc_idx, nf_u, nf_v
    real(prec), intent(in) :: dt, time
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Dilution
    real(prec), dimension(size(sv%state,1)) :: velocity_mod, volume_new, double_vol_time
    type(string) :: tag
    
    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%u
    requiredVars(2) = Globals%Var%v
    
    tag = 'initial_volume'
    init_vol_IDx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'volume'
    vol_Idx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'concentration'
    conc_idx = Utils%find_str(sv%varName, tag, .true.)
    Dilution = 0.0
    !interpolate each background
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(bdata(bkg)%hasVars(requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                nf_u = Utils%find_str(var_name, Globals%Var%u, .true.)
                nf_v = Utils%find_str(var_name, Globals%Var%v, .true.)
                ! compute velocity modulus.
                velocity_mod = sqrt(var_dt(:,nf_u)**2 + var_dt(:,nf_v)**2)
                
                !Time to double the volume (TVol200 in old lagrangian). Linear function with duplication in 30min for 2m/s and 2h for 0.1m/s.
                double_vol_time = -2840 * velocity_mod + 7485
                !Update volume
                volume_new = sv%state(:,vol_Idx) * exp((alog(2.)/double_vol_time) * dt)
                
                !Dilution in volume = Volume variation / dt
                Dilution(:, vol_Idx) = (volume_new - sv%state(:, vol_Idx)) / dt
                !Reducing concentration : dilution(dConcentration/dT) = concentration * (Initial_volume / New_Volume) / dt
                Dilution(:, conc_idx) = - (sv%state(:, conc_idx) - (sv%state(:, conc_idx) * sv%state(:, vol_Idx) / volume_new)) / dt
                deallocate(var_name)
                deallocate(var_dt) 
            end if
        end if
    end do
    
    ! Eliminate tracers with concentrations smaller than 15
    where(sv%state(:,conc_idx) < 15) sv%active = .false.
    
    end function Dilution
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelColiform(self)
    class(kernelColiform_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    call KernelUtils_coliform%initialize()
    end subroutine initKernelColiform

    end module kernelColiform_mod