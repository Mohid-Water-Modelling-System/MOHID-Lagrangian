    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernelLitter
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : September 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for litter associated processes.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach.
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
    !------------------------------------------------------------------------------

    module kernelLitter_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod


    type :: kernelLitter_class        !< Litter kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelLitter
    procedure :: DegradationLinear
    !procedure :: Buoyancy
    end type kernelLitter_class

    public :: kernelLitter_class
    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Linear degradation kernel.
    !> @param[in] self, sv
    !---------------------------------------------------------------------------
    function DegradationLinear(self, sv)
    class(kernelLitter_class), intent(in) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: DegradationLinear
    integer :: nf, idx
    type(string) :: tag

    DegradationLinear = 0.0
    tag = 'condition'
    nf = Utils%find_str(sv%varName, tag, .true.)
    tag = 'degradation_rate'
    idx = Utils%find_str(sv%varName, tag, .true.)

    DegradationLinear(:,nf) = -sv%state(:,idx)
    where(sv%state(:,nf) < 0.0) sv%active = .false.

    end function DegradationLinear

    !!---------------------------------------------------------------------------
    !!> @author Ricardo Birjukovs Canelas - MARETEC
    !!> @brief
    !!> Computes the vertical velocity due to buoyancy of the litter tracers
    !!> @param[in] self, sv, bdata, time
    !!---------------------------------------------------------------------------
    !function Buoyancy(self, sv, bdata, time)
    !class(kernelLitter_class), intent(inout) :: self
    !type(stateVector_class), intent(in) :: sv
    !type(background_class), dimension(:), intent(in) :: bdata
    !real(prec), intent(in) :: time
    !integer :: np, nf, bkg, rIdx, rhoIdx
    !real(prec), dimension(:,:), allocatable :: var_dt
    !type(string), dimension(:), allocatable :: var_name
    !type(string), dimension(:), allocatable :: requiredVars
    !real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Buoyancy
    !real(prec), dimension(size(sv%state,1)) :: fDensity, kVisco
    !real(prec) :: P = 1013.
    !type(string) :: tag
    !logical :: ViscoDensFields = .false.
    !allocate(requiredVars(2))
    !requiredVars(1) = Globals%Var%temp
    !requiredVars(2) = Globals%Var%sal
    !
    !Buoyancy = 0.0
    !!interpolate each background
    !!Compute buoyancy using state equation for temperature and viscosity
    !do bkg = 1, size(bdata)
    !    if (bdata(bkg)%initialized) then
    !        if(bdata(bkg)%hasVars(requiredVars)) then
    !            np = size(sv%active) !number of Tracers
    !            nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
    !            allocate(var_dt(np,nf))
    !            allocate(var_name(nf))
    !            !interpolating all of the data
    !            call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
    !            fDensity = seaWaterDensity(var_dt(:,2), var_dt(:,1),P)
    !            kVisco = absoluteSeaWaterViscosity(var_dt(:,2), var_dt(:,1))
    !            tag = 'radius'
    !            rIdx = Utils%find_str(sv%varName, tag, .true.)
    !            tag = 'density'
    !            rhoIdx = Utils%find_str(sv%varName, tag, .true.)
    !            ! Three different models - using the number for reynolds regime 1<Re<100
    !            !Buoyancy(:,3) = ((fDensity-sv%state(:,rhoIdx))/abs(fDensity-sv%state(:,rhoIdx)))*sqrt(8.*sv%state(:,rIdx)*abs(9.81*(sv%state(:,rhoIdx) - fDensity))/(3.*1.4*fDensity))
    !            Buoyancy(:,3) = 1.82*((sv%state(:,rIdx)*9.81*(fDensity-sv%state(:,rhoIdx)))/fDensity)**(1./2.)
    !            !Buoyancy(:,3) = (2.0/9.0)*(sv%state(:,rhoIdx) - fDensity)*Globals%Constants%Gravity%z*sv%state(:,rIdx)*sv%state(:,rIdx)/kVisco
    !
    !            deallocate(var_dt)
    !            deallocate(var_name)
    !            ViscoDensFields = .true.
    !        endif
    !    endif
    !end do
    !
    !!I there is no salt and temperatue Compute buoyancy using constant density and temp
    !! and standardad terminal velocity
    !if (ViscoDensFields == .false.) then
    !    kVisco = 1.09E-3
    !    fDensity = 1023
    !    tag = 'radius'
    !    rIdx = Utils%find_str(sv%varName, tag, .true.)
    !    tag = 'density'
    !    rhoIdx = Utils%find_str(sv%varName, tag, .true.)
    !    Buoyancy(:,3) = (2.0/9.0)*(sv%state(:,rhoIdx) - fDensity)*Globals%Constants%Gravity%z*sv%state(:,rIdx)*sv%state(:,rIdx)/kVisco
    !endif
    !
    !end function Buoyancy

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelLitter(self)
    class(kernelLitter_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initKernelLitter

    end module kernelLitter_mod