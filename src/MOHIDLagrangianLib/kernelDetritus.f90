    !------------------------------------------------------------------------------
    !        IST/MARETEC/ColabAtlantic, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernelDetritus
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2019
    ! REVISION      : Sobrinho 0.1
    !> @author
    !> Joao Sobrinho Colab Atlantic
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for Detritus associated processes.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach.
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
    !------------------------------------------------------------------------------

    module kernelDetritus_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod
    use kernelUtils_mod

    type :: kernelDetritus_class        !< Detritus kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelDetritus
    procedure :: Degradation

    end type kernelDetritus_class
    
    type(kernelUtils_class) :: KernelUtils_detritus   !< kernel utils

    public :: kernelDetritus_class
    contains
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - +Atlantic
    !> @brief
    !> Detritus degradation rate kernel.
    !> @param[in] self, sv, bdata, dt
    !---------------------------------------------------------------------------
    function Degradation(self, sv, dt)
    class(kernelDetritus_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), intent(in) :: dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Degradation
    real(prec) :: ToptBMin, ToptBMax, TBacteriaMin, TBacteriaMax, BK1, BK2, BK3, BK4, MaxDegradationRate
    real(prec) :: max_age, threshold_bot_wat
    real(prec), dimension(size(sv%state,1)) :: s1, s2, ya, yb, xa, xb, limFactor, mass, volume_new, init_mass, temperature
    integer :: volume_col, radius_col, col_temp, initvol_col, density_col, age_col
    real(prec):: compute_rate = 7200
    type(string) :: tag
    !Begin------------------------------------------------------------------------
    Degradation = 0.0
    tag = 'volume'
    volume_col = Utils%find_str(sv%varName, tag, .true.)
    tag = 'density'
    density_col = Utils%find_str(sv%varName, tag, .true.)
    tag = 'radius'
    radius_col = Utils%find_str(sv%varName, tag, .true.)
    tag = 'initial_volume'
    initvol_col = Utils%find_str(sv%varName, tag, .true.)
    tag = 'age'
    age_col = Utils%find_str(sv%varName, tag, .true.)
    
    !Compute only when time of the oldest particle is a multiple of 7200s, to reduce simulation time
    max_age = maxval(sv%state(:,age_col))
    if (mod(max_age,compute_rate) < 1E-3) then
        !Method used in bacterial growth limitation by temperature in MOHIDWater
        ToptBMin = Globals%Constants%TOptBacteriaMin
        ToptBMax = Globals%Constants%TOptBacteriaMax
        TBacteriaMin = Globals%Constants%TBacteriaMin
        TBacteriaMax = Globals%Constants%TBacteriaMax
        BK1 = Globals%Constants%BK1
        BK2 = Globals%Constants%BK2
        BK3 = Globals%Constants%BK3
        BK4 = Globals%Constants%BK4
        MaxDegradationRate = Globals%Constants%MaxDegradationRate
        
        col_temp = Utils%find_str(sv%varName, Globals%Var%temp, .true.)
        threshold_bot_wat = (Globals%Mask%waterVal + Globals%Mask%bedVal) * 0.5
        temperature = sv%state(:,col_temp)
        
        mass = sv%state(:,density_col) * sv%state(:,volume_col)
        init_mass = sv%state(:,density_col) * sv%state(:,initvol_col)
    
        !Change radius keeping shperical form
        sv%state(:,radius_col) = (sv%state(:,volume_col)*(0.75)*(1.0/3.14159265))**(1.0/3.0)
        
        s1 = (1.0 / (ToptBMin - TBacteriaMin)) * dlog((BK2 * (1.0 - BK1)) / (BK1 * (1.0 - BK2)))
        s2 = (1.0 / (TBacteriaMax - ToptBMax)) * dlog((BK3 * (1.0 - BK4))  / (BK4 * (1.0 - BK3)))
    
        ya = exp(s1 * (temperature - TBacteriaMin))
        yb = exp(s2 * (TBacteriaMax - temperature))
    
        xa = (BK1 * ya) / (1 + BK1 * (ya - 1))
        xb = (BK4 * yb) / (1 + BK4 * (yb - 1))
    
        limFactor = xa * xb
        !Kg = Kg    -  Kg  *    []     *          1/s       * s
        mass = mass - init_mass * limFactor * MaxDegradationRate * compute_rate
    
        !keeping density untouched, but changing volume
        volume_new = max(mass / sv%state(:,density_col), 0.0)
    
        Degradation(:,volume_col) = - (sv%state(:,volume_col) - volume_new) / dt
        
        !where(sv%state(:,volume_col) < 0.02*sv%state(:,initvol_col)) sv%active = .false.
        where(mass < 0.002) sv%active = .false.
        
    end if
    
    end function Degradation
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelDetritus(self)
    class(kernelDetritus_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    call KernelUtils_detritus%initialize()
    end subroutine initKernelDetritus

    end module kernelDetritus_mod