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
    function Degradation(self, sv, bdata, time, dt)
    class(kernelDetritus_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Degradation
    real(prec) :: ToptBMin, ToptBMax, TBacteriaMin, TBacteriaMax, BK1, BK2, BK3, BK4, MaxDegradationRate, multiple, max_age
    real(prec), dimension(size(sv%state,1)) :: s1, s2, ya, yb, xa, xb, limFactor, mass, volume_new, init_mass
    type(string), dimension(:), allocatable :: requiredVars
    integer :: volume_col, radius_col, temp_col, initvol_col, density_col, age_col
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
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
    
    allocate(requiredVars(1))
    requiredVars(1) = Globals%Var%temp
        
    call KernelUtils_detritus%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name)
        
    temp_col = Utils%find_str(var_name, Globals%Var%temp, .true.)
    
    !Compute only when time of the oldest particle is a multiple of 3600s, to reduce simulation time
    multiple = 7200 / dt
    max_age = maxval(sv%state(:,age_col))
    if (mod(max_age,multiple) < 1E-3) then
        
        mass = sv%state(:,density_col) * sv%state(:,volume_col)
        init_mass = sv%state(:,density_col) * sv%state(:,initvol_col)
    
        !Change radius keeping shperical form
        sv%state(:,radius_col) = (sv%state(:,volume_col)*(0.75)/3.1415926)**(1/3)
    
        s1 = (1 / (ToptBMin - TBacteriaMin)) * dlog((BK2 * (1 - BK1)) / (BK1 * (1 - BK2)))
        s2 = (1 / (TBacteriaMax - ToptBMax)) * dlog((BK3 * (1 - BK4))  / (BK4 * (1 - BK3)))
    
        ya = exp(s1 * (var_dt(:,temp_col) - TBacteriaMin))
        yb = exp(s2 * (TBacteriaMax - var_dt(:,temp_col)))
    
        xa = (BK1 * ya) / (1 + BK1 * (ya - 1))
        xb = (BK4 * yb) / (1 + BK4 * (yb - 1))
    
        limFactor = xa * xb
        !Kg = Kg    -  Kg  *    []     *          1/s       * s
        mass = mass - init_mass * limFactor * MaxDegradationRate * dt
    
        !keeping density untouched, but changing volume
        volume_new = max(mass / sv%state(:,density_col), 0.0)
    
        Degradation(:,volume_col) = - (sv%state(:,volume_col) - volume_new) / dt        
    end if
    
    deallocate(var_name)
    deallocate(var_dt)
        
    !matar Detritus com menos de 100 UFC/100ml
    where(sv%state(:,volume_col) < 0.02*sv%state(:,initvol_col)) sv%active = .false.
    
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