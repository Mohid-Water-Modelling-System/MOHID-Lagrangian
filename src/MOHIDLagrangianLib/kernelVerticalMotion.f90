    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !        USC/GFNL, Group of NonLinear Physics, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : kernelVerticalMotion_mod
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : USC/MARETEC, Marine Modelling Group
    ! DATE          : November 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for buoyancy associated processes.
    !> This class has several methods, that should be designed on a one method - one
    !> process approach.
    !> The output of every kernel should be a 2D matrix, where a row represents the
    !> derivative of the state vector of a given tracer. n columns - n variables.
    !> This is the step were interpolation and physics actually happen.
    !------------------------------------------------------------------------------

    module kernelVerticalMotion_mod

    use common_modules
    use stateVector_mod
    use background_mod
    use interpolator_mod

    ! Density constants
    !  density of seawater at zero pressure constants
    real(prec),parameter :: a0 = 999.842594
    real(prec),parameter :: a1 =   6.793952e-2
    real(prec),parameter :: a2 =  -9.095290e-3
    real(prec),parameter :: a3 =   1.001685e-4
    real(prec),parameter :: a4 =  -1.120083e-6
    real(prec),parameter :: a5 =   6.536332e-9

    real(prec),parameter :: b0 =   8.24493e-1
    real(prec),parameter :: b1 =  -4.0899e-3
    real(prec),parameter :: b2 =   7.6438e-5
    real(prec),parameter :: b3 =  -8.2467e-7
    real(prec),parameter :: b4 =   5.3875e-9

    real(prec),parameter :: c0 =  -5.72466e-3
    real(prec),parameter :: c1 =   1.0227e-4
    real(prec),parameter :: c2 =  -1.6546e-6

    real(prec),parameter :: d0 =   4.8314e-4

    real(prec),parameter :: h0 =  3.239908
    real(prec),parameter :: h1 =  1.43713E-3
    real(prec),parameter :: h2 =  1.16092E-4
    real(prec),parameter :: h3 = -5.77905E-7
    real(prec),parameter :: k0 =  8.50935E-5
    real(prec),parameter :: k1 = -6.12293E-6
    real(prec),parameter :: k2 =  5.2787E-8
    real(prec),parameter :: e0 = 19652.21
    real(prec),parameter :: e1 = 148.4206
    real(prec),parameter :: e2 = -2.327105
    real(prec),parameter :: e3 =  1.360477E-2
    real(prec),parameter :: e4 = -5.155288E-5
    real(prec),parameter :: i0 =  2.2838E-3
    real(prec),parameter :: i1 = -1.0981E-5
    real(prec),parameter :: i2 = -1.6078E-6
    real(prec),parameter :: j0 =  1.91075E-4
    real(prec),parameter :: f0 = 54.6746
    real(prec),parameter :: f1 = -0.603459
    real(prec),parameter :: f2 =  1.09987e-2
    real(prec),parameter :: f3 = -6.1670e-5
    real(prec),parameter :: g0 =  7.944e-2
    real(prec),parameter :: g1 =  1.6483e-2
    real(prec),parameter :: g2 = -5.3009e-4
    real(prec),parameter :: m0 = -9.9348E-7
    real(prec),parameter :: m1 =  2.0816E-8
    real(prec),parameter :: m2 =  9.1697E-10

    ! viscosity Constants constants
    real(prec),parameter :: n1 = -3.79418
    real(prec),parameter :: n2 = 604.129
    real(prec),parameter :: n3 = 139.18
    real(prec),parameter :: o1 = 1.474E-3
    real(prec),parameter :: o2 = 1.5E-5
    real(prec),parameter :: o3 = -3.927E-8
    real(prec),parameter :: p1 = 1.0734E-5
    real(prec),parameter :: p2 = -8.5E-8
    real(prec),parameter :: p3 = 2.23E-10


    type :: kernelVerticalMotion_class        !< VerticalMotion kernel class
        type(interpolator_class) :: Interpolator !< The interpolator object for the kernel
    contains
    procedure :: initialize => initKernelVerticalMotion
    procedure :: Buoyancy
    procedure :: CorrectVerticalBounds
    end type kernelVerticalMotion_class

    public :: kernelVerticalMotion_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Computes the vertical velocity due to buoyancy of the tracers
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Buoyancy(self, sv, bdata, time)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: np, nf, bkg, rIdx, rhoIdx
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Buoyancy
    real(prec), dimension(size(sv%state,1)) :: fDensity, kVisco
    real(prec) :: P = 1013.
    type(string) :: tag
    logical :: ViscoDensFields = .false.
    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%temp
    requiredVars(2) = Globals%Var%sal
    
    Buoyancy = 0.0
    !interpolate each background
    !Compute buoyancy using state equation for temperature and viscosity
    do bkg = 1, size(bdata)
        if (bdata(bkg)%initialized) then
            if(bdata(bkg)%hasVars(requiredVars)) then
                np = size(sv%active) !number of Tracers
                nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                allocate(var_dt(np,nf))
                allocate(var_name(nf))
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                fDensity = seaWaterDensity(var_dt(:,2), var_dt(:,1),P)
                kVisco = absoluteSeaWaterViscosity(var_dt(:,2), var_dt(:,1))
                tag = 'radius'
                rIdx = Utils%find_str(sv%varName, tag, .true.)
                tag = 'density'
                rhoIdx = Utils%find_str(sv%varName, tag, .true.)
                ! Three different models - using the number for reynolds regime 1<Re<100
                !Buoyancy(:,3) = ((fDensity-sv%state(:,rhoIdx))/abs(fDensity-sv%state(:,rhoIdx)))*sqrt(8.*sv%state(:,rIdx)*abs(9.81*(sv%state(:,rhoIdx) - fDensity))/(3.*1.4*fDensity))
                Buoyancy(:,3) = 1.82*((sv%state(:,rIdx)*9.81*(fDensity-sv%state(:,rhoIdx)))/fDensity)**(1./2.)
                !Buoyancy(:,3) = (2.0/9.0)*(sv%state(:,rhoIdx) - fDensity)*Globals%Constants%Gravity%z*sv%state(:,rIdx)*sv%state(:,rIdx)/kVisco
                deallocate(var_dt)
                deallocate(var_name)
                ViscoDensFields = .true.
            endif
        endif
    end do

    !I there is no salt and temperatue Compute buoyancy using constant density and temp
    ! and standardad terminal velocity
    if (.not.ViscoDensFields) then
        kVisco = 1.09E-3
        fDensity = 1023
        tag = 'radius'
        rIdx = Utils%find_str(sv%varName, tag, .true.)
        tag = 'density'
        rhoIdx = Utils%find_str(sv%varName, tag, .true.)
        Buoyancy(:,3) = (2.0/9.0)*(sv%state(:,rhoIdx) - fDensity)*Globals%Constants%Gravity%z*sv%state(:,rIdx)*sv%state(:,rIdx)/kVisco
    endif
    
    end function Buoyancy

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Corrects vertical position of the tracers according to data limits
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function CorrectVerticalBounds(self, sv, svDt, bdata, time, dt)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: svDt
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: CorrectVerticalBounds
    integer :: np, nf, bkg
    real(prec) :: maxLevel(2)
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    
    CorrectVerticalBounds = svDt

    allocate(requiredVars(1))
    requiredVars(1) = Globals%Var%bathymetry

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
                if (maxLevel(2) /= MV) where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt > maxLevel(2)) CorrectVerticalBounds(:,3) = (maxLevel(2)-sv%state(:,3))/dt*0.99

                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)

                !update minium vertical position
                nf = Utils%find_str(var_name, Globals%Var%bathymetry)
                where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt < var_dt(:,nf)) CorrectVerticalBounds(:,3) = (var_dt(:,nf)-sv%state(:,3))/dt*0.99

                deallocate(var_dt)
                deallocate(var_name)
            end if
        end if
    end do

    end function CorrectVerticalBounds

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - GFNL
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
    subroutine initKernelVerticalMotion(self)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(string) :: interpName
    interpName = 'linear'
    call self%Interpolator%initialize(1,interpName)
    end subroutine initKernelVerticalMotion

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !orphan functions, enter at your own risk

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Method that returns density of seawater at zero pressure
    !---------------------------------------------------------------------------
    function seawaterDensityZeroPressure(S,T)
    real(prec), dimension(:),intent(in) :: S,T
    real(prec), dimension(size(S)) :: SMOW,RB,RC
    real(prec), dimension(size(S)) :: seawaterDensityZeroPressure
    ! --- Computations ---
    ! Density of pure water
    SMOW = a0 + (a1 + (a2 + (a3 + (a4 + a5*T)*T)*T)*T)*T
    ! More temperature polynomials
    RB = b0 + (b1 + (b2 + (b3 + b4*T)*T)*T)*T
    RC = c0 + (c1 + c2*T)*T
    seawaterDensityZeroPressure = SMOW + RB*S + RC*(S**1.5) + d0*S*S
    end function seawaterDensityZeroPressure

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Method that returns secant bulk modulus
    !---------------------------------------------------------------------------
    function secantBulkModulus(S, T, P)
    real(prec), dimension(:),intent(in) :: S,T
    real(prec),intent(in) :: P
    real(prec), dimension(size(S)):: AW,BW,KW,SR,A,K0,B
    real(prec),dimension(size(S)):: secantBulkModulus
    !--- Pure water terms ---
    AW = h0 + (h1 + (h2 + h3*T)*T)*T
    !--- seawater, P = 0 ---
    BW = k0 + (k1 + k2*T)*T
    KW = e0 + (e1 + (e2 + (e3 + e4*T)*T)*T)*T
    SR = S**0.5
    A  = AW + (i0 + (i1 + i2*T)*T + j0*SR)*S
    K0 = KW + (f0 + (f1 + (f2 + f3*T)*T)*T + (g0 + (g1 + g2*T)*T)*SR)*S
    ! --- General expression ---
    B = BW + (m0 + (m1 + m2*T)*T)*S
    secantBulkModulus = K0 + (A + B*P)*P
    end function secantBulkModulus

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Method that returns secant bulk modulus
    !---------------------------------------------------------------------------
    function seaWaterDensity(S, T, P)
    real(prec), dimension(:),intent(in) :: S,T
    real(prec),intent(in) :: P
    real(prec) :: Pbar
    real(prec),dimension(size(S)) :: seaWaterDensity
    !Compute density of seawater from salinity, temperature, and pressure
    !Usage: dens(S, T, [P])
    !Input:
    !    S = Salinity,     [PSS-78]
    !    T = Temperature,  [�C]
    !    P = Pressure,     [dbar = 10**4 Pa]
    !P is optional, with default value zero
    !Output:
    !    Density,          [kg/m**3]
    !Algorithm: UNESCO 1983
    !"""
    Pbar = 0.1*P !# Convert to bar
    seaWaterDensity= seawaterDensityZeroPressure(S,T)/(1 - Pbar/secantBulkModulus(S,T,Pbar))
    end function seaWaterDensity

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Method that returns secant bulk modulus
    !---------------------------------------------------------------------------
    function absoluteSeaWaterViscosity(S,T)
    real(prec), dimension(:),intent(in) :: S,T
    real(prec), dimension(size(S)):: mu_W,mu_R,A,B
    real(prec),dimension(size(S)):: absoluteSeaWaterViscosity
    ! The dynamic viscosity correlation of sea water is given by:
    ! [kg/m s]
    ! nu  = mu / density			[m²/s]
    ! El-Dessouky, Ettouny (2002): Fundamentals of Sea Water Desalination (Appendix A: Themodynamic Properties)

    mu_W =	exp(n1 + n2/(n3 +T))
    A = 	o1 + o2*T + o3*T*T
    B = 	p1 + p2*T + p3*T*T
    mu_R =	1 + A*S + B*S*S

    absoluteSeaWaterViscosity  = mu_w*mu_R*0.001
    end function absoluteSeaWaterViscosity

    end module kernelVerticalMotion_mod