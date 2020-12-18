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
    procedure :: Reynolds
    procedure :: DragCoefficient
    procedure :: SphericalShapeFactor
    procedure :: Divergence
    procedure :: Resuspension
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
    integer :: np, nf, bkg, rIdx, rhoIdx, areaIdx, volIdx, temp1, temp2
    real(prec), dimension(:,:), allocatable :: var_dt
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Buoyancy
    real(prec), dimension(size(sv%state,1)) :: fDensity, kVisco
    real(prec), dimension(size(sv%state,1)) :: signZ, shapeFactor, densityRelation, cd,Re
    real(prec), dimension(size(sv%state,1)) :: ReynoldsNumber, kViscoRelation
    real(prec), dimension(2) :: maxLevel
    real(prec) :: P = 1013.

    type(string) :: tag
    logical :: ViscoDensFields = .false.

    Buoyancy = 0.0
    
    allocate(requiredVars(2))
    requiredVars(1) = Globals%Var%temp
    requiredVars(2) = Globals%Var%sal

    !maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
    !if (maxLevel(2) /= MV) then
    !    return
    !end if
    
    !interpolate each background
    !Compute buoyancy using state equation for temperature and viscosity
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
                temp1 = Utils%find_str(var_name, tag, .true.)
                tag = Globals%Var%sal
                temp2 = Utils%find_str(var_name, tag, .true.)
                fDensity = seaWaterDensity(var_dt(:,temp2), var_dt(:,temp1),P)
                kVisco = absoluteSeaWaterViscosity(var_dt(:,temp2), var_dt(:,temp1))
                ! Viscosity on boundaries could be 0. Avoid overdumped values on viscosity
                ! Viscosity relation of 90 % related to mean water density are allowed.
                kViscoRelation = abs(1.-(kvisco/Globals%Constants%MeanKvisco))
                where(kViscoRelation >= 0.9)
                    kVisco = Globals%Constants%MeanKvisco
                    fDensity = Globals%Constants%MeanDensity
                endwhere
                ! Read variables
                tag = 'radius'
                rIdx = Utils%find_str(sv%varName, tag, .true.)
                tag = 'density'
                rhoIdx = Utils%find_str(sv%varName, tag, .true.)
                tag = 'volume'
                volIdx = Utils%find_str(sv%varName, tag, .true.)
                tag = 'area'
                areaIdx = Utils%find_str(sv%varName, tag, .true.)
                ! Get the direction of buoyancy ( + to the surface, - to the bottom)
                signZ = -1.0
                where(fDensity-sv%state(:,rhoIdx) >= 0)
                    signZ = 1.0
                endwhere
                ! Boundary density values could be 0. This avoid underdumped values on density. 
                ! Just density relation of 90 % related to mean water density are allowed.
                densityRelation = abs(1.- (sv%state(:,rhoIdx)/fDensity))
                where (densityRelation >= 0.9)
                    densityRelation = abs(1.- (sv%state(:,rhoIdx)/Globals%Constants%MeanDensity))
                endwhere    
                ! Get drag and shapefactor
                shapeFactor = self%SphericalShapeFactor(sv%state(:,areaIdx),sv%state(:,volIdx))
                reynoldsNumber = self%Reynolds(sv%state(:,3), kvisco, sv%state(:,rIdx)) 
                cd = self%dragCoefficient(shapeFactor, sv%state(:,rIdx), reynoldsNumber) 
                ! Get buoyancy
                where (reynoldsNumber /=0.)
                    Buoyancy(:,3) = signZ*sqrt((-2.*Globals%Constants%Gravity%z) * (shapeFactor/cd) * densityRelation)
                endwhere
                deallocate(var_name)
                deallocate(var_dt)
                ViscoDensFields = .true.
            end if
        end if
    end do
    
    !I there is no salt and temperatue Compute buoyancy using constant density and temp
    ! and standardad terminal velocity
    if (.not. ViscoDensFields) then

        ! Check if your domain data has depth level to move in vertical motion.
        requiredVars(1)=Globals%var%u
        requiredVars(2)=Globals%var%v
        do bkg = 1, size(bdata)
            if (bdata(bkg)%initialized) then
                if (bdata(bkg)%hasVars(requiredVars)) then
                    maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
                    if (maxLevel(2) == MV) then
                        return
                    end if
                end if
            end if
        end do

        ! interpolate each background
        ! Compute buoyancy using state equation for temperature and viscosity

        kVisco = Globals%Constants%MeanKvisco
        fDensity = Globals%Constants%MeanDensity
        tag = 'radius'
        rIdx = Utils%find_str(sv%varName, tag, .true.)
        tag = 'density'
        rhoIdx = Utils%find_str(sv%varName, tag, .true.)
        tag = 'volume'
        volIdx = Utils%find_str(sv%varName, tag, .true.)
        tag = 'area'
        areaIdx = Utils%find_str(sv%varName, tag, .true.)
        signZ = -1.0
        where(fDensity-sv%state(:,rhoIdx) >= 0)
            signZ = 1.0
        endwhere
       
        ! Boundary density values could be 0. This avoid underdumped values on density. 
        ! Just density relation of 90 % related to mean water density are allowed.
        densityRelation = abs(1.- (sv%state(:,rhoIdx)/fDensity))          
        ! Get drag and shapefactor
        shapeFactor = self%SphericalShapeFactor(sv%state(:,areaIdx),sv%state(:,volIdx))
        reynoldsNumber = self%Reynolds(sv%state(:,3), kvisco, sv%state(:,rIdx)) 
        cd = self%dragCoefficient(shapeFactor, sv%state(:,rIdx), reynoldsNumber) 
        ! Get buoyancy
        where (reynoldsNumber /=0.)
            Buoyancy(:,3) = signZ*sqrt((-2.*Globals%Constants%Gravity%z) * (shapeFactor/cd) * densityRelation)
        end where                
    end if

    end function Buoyancy


    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Get the particle reynolds number
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Reynolds(self, characteristicVelocity, kvisco, characteristicLength)
    class(kernelVerticalMotion_class), intent(inout) :: self
    real(prec), dimension(:) , intent(in) :: characteristicVelocity, kvisco, characteristicLength
    real(prec), dimension(size(characteristicVelocity)) :: Reynolds
    
    Reynolds = abs(characteristicVelocity)*characteristicLength/kvisco

    end function Reynolds

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Approach the particle drag coefficient based on shape factor
    !> @param[in] self, sv, shapeFactor, characteristicLength 
    !---------------------------------------------------------------------------
    function DragCoefficient(self, shapeFactor, charateristicLength, Reynolds)
    class(kernelVerticalMotion_class), intent(inout) :: self
    real(prec), dimension(:), intent(in):: shapeFactor, charateristicLength, reynolds
    real(prec), dimension(size(shapeFactor,1)) :: DragCoefficient

    dragCoefficient = (24.00/Reynolds)*(1. + 0.173*Reynolds**0.657) + 0.413/(1.+16300*Reynolds**(-1.09))

    ! Alternative implementation    
    !where (Reynolds <= 0.2)
    !    dragCoefficient = (24.00/Reynolds)
    !elsewhere((Reynolds < 0.2) * (Reynolds <= 1000))
    !    dragCoefficient = (24.00/Reynolds) + 3./sqrt(Reynolds) + 0.34
    !elsewhere((Reynolds <= 1000) * (Reynolds <= 100000))
    !    dragCoefficient = 0.44
    !elsewhere
    !    dragCoefficient = 0.2
    !endwhere

    end function DragCoefficient

        !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> check the buoyancy for abnormal values
    !> @param[in] self, sv, shapeFactor, characteristicLength 
    !---------------------------------------------------------------------------
    function checkBuoyancy(self, sv, values)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    real(prec), dimension(size(sv%state,1)) :: values
    logical :: checkBuoyancy

    if (any(abs(values) > 1)) then
        print*, '>> WARNING: Abnormal vertical velocity detected'
        checkBuoyancy = .false.
    else
        checkBuoyancy = .true.
    end if
        
    end function checkBuoyancy

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC 
    !> @brief
    !> Gets the spherical object shapefactor
    !> @param[in] self, sv, bdata, time
    !--------------------------------------------------------------------------- 
    function SphericalShapeFactor(self, area, volume)
    class(kernelVerticalMotion_class), intent(inout) :: self
    real(prec), dimension(:), intent(in):: area, volume
    real(prec), dimension(size(area)) :: SphericalShapeFactor
    real(prec) :: pi
    pi = 4.*atan(1.)
    SphericalShapeFactor = (6.*volume/pi)**(1./3.)/sqrt(4.*area/pi)
    end function SphericalShapeFactor

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
    
    ! if vertical dvdt is 0, don't correct
    if (all(CorrectVerticalBounds(:,3) == 0.)) then
        return
    end if

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
                if (maxLevel(2) /= MV) where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt >= maxLevel(2)) CorrectVerticalBounds(:,3) =((maxLevel(2)-sv%state(:,3))/dt)*0.99
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                !update minium vertical position
                nf = Utils%find_str(var_name, Globals%Var%bathymetry)
                where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt < var_dt(:,nf)) CorrectVerticalBounds(:,3) = ((var_dt(:,nf)-sv%state(:,3))/dt)*0.99
                deallocate(var_name)
                deallocate(var_dt)
            end if
        end if
    end do

    end function CorrectVerticalBounds


    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - MARETEC
    !> @brief
    !> Resuspend particles based on horizontal velocity module motion. 
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Resuspension(self, sv, bdata, time)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Resuspension
    integer :: np, nf, bkg, nf_u, nf_v, nf_lim
    real(prec) :: maxLevel(2), landIntThreshold
    real(prec), dimension(:,:), allocatable :: var_dt
    real(prec), dimension(size(sv%state,1)) :: velocity_mod
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    
    Resuspension = 0.0
    landIntThreshold = -0.75

    if (Globals%Constants%ResuspensionCoeff > 0.0) then

        allocate(requiredVars(3))
        requiredVars(1) = Globals%Var%u
        requiredVars(2) = Globals%Var%v
        requiredVars(3) = Globals%Var%landIntMask

        !interpolate each background
        do bkg = 1, size(bdata)
            if (bdata(bkg)%initialized) then
                if(bdata(bkg)%hasVars(requiredVars)) then
                    maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
                    if (maxLevel(2) == MV) then
                        return ! if it is 2 dimensional data, return
                    else if (maxLevel(2) /= MV) then
                        np = size(sv%active) !number of Tracers
                        nf = bdata(bkg)%fields%getSize() !number of fields to interpolate
                        allocate(var_dt(np,nf))
                        allocate(var_name(nf))
                        !interpolating all of the data
                        call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                        nf_u = Utils%find_str(var_name, Globals%Var%u, .true.)
                        nf_v = Utils%find_str(var_name, Globals%Var%v, .true.)
                        nf_lim = Utils%find_str(var_name, Globals%Var%landIntMask, .true.)
                        ! compute velocity modulus.
                        velocity_mod = sqrt(var_dt(:,nf_u)**2 + var_dt(:,nf_v)**2)
                        ! Resuspension apply if level data has more than 1 dimension.
                        ! Resuspension apply based on threshold.
                        ! Tracer jumps in the w equals to horizontal velocity.
                        where (var_dt(:,nf_lim) < landIntThreshold)  Resuspension(:,3) = Globals%Constants%ResuspensionCoeff*velocity_mod
                        deallocate(var_name)
                        deallocate(var_dt) 
                    end if
                end if
            end if
        end do
    end if
    
    end function Resuspension
    
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
    where (seaWaterDensity /= seaWaterDensity)
        seaWaterDensity = Globals%Constants%MeanDensity
    end where
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
    where (absoluteSeaWaterViscosity /= absoluteSeaWaterViscosity)
        absoluteSeaWaterViscosity = Globals%Constants%MeanKvisco
    end where
    end function absoluteSeaWaterViscosity


    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> WIP: Divergence Kernel, estimate the vertical velocity based on divergence.
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Divergence(self, sv, bdata, time)
        class(kernelVerticalMotion_class), intent(inout) :: self
        type(stateVector_class), intent(inout) :: sv
        type(stateVector_class) :: x0, x1, y0, y1
        type(background_class), dimension(:), intent(in) :: bdata
        real(prec), intent(in) :: time 
        real(prec) :: dr = 0.01
        integer :: np, nf, bkg, i
        real(prec), dimension(:,:), allocatable :: v1,v0,u1,u0, var_dt
        type(string), dimension(:), allocatable :: var_name
        type(string), dimension(:), allocatable :: requiredVars

        real(prec), dimension(size(sv%state,1)) :: Divergence, resolution
        real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: uv_x0, uv_x1, uv_y0, uv_y1
        real(prec), dimension(size(sv%state,1)) :: u_x0, u_x1, v_y0, v_y1, dx,dy
        
        allocate(requiredVars(3))
        requiredVars(1) = Globals%Var%u
        requiredVars(2) = Globals%Var%v
        requiredVars(3) = Globals%Var%resolution

        Divergence = 0.

        ! Creating points to compute derivative.
        call sv%copyState(x0)
        call sv%copyState(x1)
        call sv%copyState(y0)
        call sv%copyState(y1)

        ! interpolate each background
        do bkg = 1, size(bdata)
            if (bdata(bkg)%initialized) then
                if(bdata(bkg)%hasVars(requiredVars)) then
                    np = size(sv%active)             ! number of Tracers
                    nf = bdata(bkg)%fields%getSize() ! number of fields to interpolate
                    allocate(var_name(nf))
                    allocate(var_dt(np,nf))

                    call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                    
                    nf = Utils%find_str(var_name, Globals%Var%resolution, .true.)
                    resolution = var_dt(:,nf)
                    !if we are still in the same path, use the same random velocity, do nothing
                    !if we ran the path, new random velocities are generated and placed
  
                    dx = Utils%m2geo(resolution, sv%state(:,2), .false.)
                    dy = Utils%m2geo(resolution, sv%state(:,2), .true.)

                    x0%state(:,1) = x0%state(:,1) - dx
                    x1%state(:,1) = x1%state(:,1) + dx
                    y0%state(:,2) = y0%state(:,2) - dy
                    y1%state(:,2) = y1%state(:,2) + dy
            
                    call self%Interpolator%run(x0%state, bdata(bkg), time, uv_x0, var_name)
                    nf = Utils%find_str(var_name, Globals%Var%u, .true.)
                    u_x0 = Utils%m2geo(uv_x0(:,nf), y0%state(:,2), .false.)

                    call self%Interpolator%run(x1%state, bdata(bkg), time, uv_x1, var_name)
                    nf = Utils%find_str(var_name, Globals%Var%u, .true.)
                    u_x1 = Utils%m2geo(uv_x1(:,nf), y1%state(:,2), .false.)

                    call self%Interpolator%run(y0%state, bdata(bkg), time, uv_y0, var_name)
                    nf = Utils%find_str(var_name, Globals%Var%v, .true.)
                    v_y0 = Utils%m2geo(uv_y0(:,nf), y0%state(:,2), .true.)

                    call self%Interpolator%run(y1%state, bdata(bkg), time, uv_y1, var_name)
                    nf = Utils%find_str(var_name, Globals%Var%v, .true.)
                    v_y1 = Utils%m2geo(uv_y1(:,nf), y1%state(:,2), .true.)

                    deallocate(var_dt)
                    deallocate(var_name)
                end if
            end if
        end do
        
        Divergence = -((u_x1-u_x0)/(2.*dr) + (v_y1-v_y0)/(2.*dr))

    end function Divergence
     

    end module kernelVerticalMotion_mod