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


    type :: kernelColiform_class        !< Coliform kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelColiform
    procedure :: MortalityT90
    !procedure :: Buoyancy
    end type kernelColiform_class

    public :: kernelColiform_class
    contains
    
    !---------------------------------------------------------------------------
    !> @author Joao Barros Sobrinho - +Atlantic
    !> @brief
    !> T_90 fecal coliforms decay kernel.
    !> @param[in] self, sv, bdata, dt
    !---------------------------------------------------------------------------
    function MortalityT90(self, sv, bdata, time)
    class(kernelColiform_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: MortalityT90
    real(prec), dimension(size(sv%state,1)) :: depth, Radiation_SW, T90, Radiation
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(2) :: maxLevel
    integer :: nf, SWcoefidx, SWperidx, np, bkg, temp, sal, rad, T90_method, conc_idx, T90_idx, T90_var_idx, T90_varM_idx
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string) :: tag
    logical :: T_S_Rad_Fields
    !Begin------------------------------------------------------------------------
    
    MortalityT90 = 0.0
    tag = 'concentration'
    conc_idx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'T90'
    T90_idx = Utils%find_str(sv%varName, tag, .true.)
    
    tag = 'T90_variable'
    T90_var_idx = Utils%find_str(sv%varName, tag, .true.)
    
    tag = 'T90_method'
    T90_varM_idx = Utils%find_str(sv%varName, tag, .true.)
    
    if (all(sv%state(:,T90_varM_idx) == 1)) then
        T90_method = 1
    else
        T90_method = 2
    end if
    
    if (all(sv%state(:,T90_var_idx) == .true.)) then
        allocate(requiredVars(3))
        requiredVars(1) = Globals%Var%temp
        requiredVars(2) = Globals%Var%sal
        requiredVars(3) = Globals%Var%rad
        T_S_Rad_Fields = .false.
        do bkg = 1, size(bdata)
            if (bdata(bkg)%initialized) then
                if (bdata(bkg)%hasVars(requiredVars)) then
                    np = size(sv%active) !number of Tracers
                    nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                    allocate(var_dt(np,nf))
                    allocate(var_name(nf))
                    !interpolating all of the data
                    call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                    tag = Globals%Var%temp
                    temp = Utils%find_str(var_name, tag, .true.)
                    tag = Globals%Var%sal
                    sal = Utils%find_str(var_name, tag, .true.)
                    !surface radiation
                    tag = Globals%Var%rad
                    rad = Utils%find_str(var_name, tag, .true.)
                    
                    !Get SW percentage
                    tag = 'sw_percentage'
                    SWperidx = Utils%find_str(sv%varName, tag, .true.)
                    !w/m2        =   w/m2          *     []
                    Radiation_SW = var_dt(:,rad) * sv%state(:,SWperidx)
                    
                    !compute light extintion
                    depth = sv%state(:,3)
                    
                    maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
                    
                    depth = maxLevel(2) - depth
                    
                    !Get SW extinction coef value
                    tag = 'sw_extinction_coef'
                    SWcoefidx = Utils%find_str(sv%varName, tag, .true.)
                    
                    !compute light exctintion in water column
                    Radiation = Radiation_SW * exp(-sv%state(:,SWcoefidx) * depth)
                    
                    !Compute T90
                    if (T90_method == 1) then
                        !Canteras
                
                        T90 = 2.533 * (1.04**(var_dt(:,temp) - 20.)) * (1.012**var_dt(:,sal)) + (0.113 * var_dt(:,rad))
                        
                        ![]                        
                        T90 = (2.303 / T90) * 24 * 3600
                        
                        !CurrentPartic%Concentration(nProp) = CurrentPartic%Concentration(nProp) / &
                        !                                                                 (1.0 + Me%DT_Partic *     &
                        !                                                                 (log(10.) / CurrentPartic%T90))
                        
                        
                        MortalityT90(:,conc_idx)  =  sv%state(:,conc_idx) * (alog(10.0) / T90)

                    else if (T90_method == 2) then
                        !Chapra not yet implemented
            
                    end if
                    
                    deallocate(var_name)
                    deallocate(var_dt)
                    T_S_Rad_Fields = .true.
                end if
            end if
        end do
        
        if (.not. T_S_Rad_Fields) then
            !Use constant T90
            MortalityT90(:,conc_idx) = - (sv%state(:,conc_idx) * sv%state(:,T90_idx))
        endif
        
    else
        MortalityT90(:,conc_idx) = - (sv%state(:,conc_idx) * sv%state(:,T90_idx))
    end if
    
    !matar coliformes com menos de 0.01 UFC/100ml
    where(sv%state(:,conc_idx) < 10) sv%active = .false.

    end function MortalityT90
    
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
    end subroutine initKernelColiform

    end module kernelColiform_mod