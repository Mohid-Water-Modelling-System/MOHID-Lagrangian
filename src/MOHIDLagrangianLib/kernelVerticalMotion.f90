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
    use kernelUtils_mod

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
    
    type(kernelUtils_class) :: KernelUtils_VerticalMotion   !< kernel utils

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
    integer :: rIdx, rhoIdx, areaIdx, volIdx
    integer :: col_temp, col_sal, col_dwz, col_bat
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Buoyancy
    real(prec), dimension(size(sv%state,1)) :: fDensity, kVisco, dist2bottom
    real(prec), dimension(size(sv%state,1)) :: signZ, shapeFactor, densityRelation, cd,Re
    real(prec), dimension(size(sv%state,1)) :: ReynoldsNumber, kViscoRelation
    real(prec), dimension(2) :: maxLevel
    real(prec) :: landIntThreshold = -0.98
    type(string) :: tag
    ! Begin -------------------------------------------------------------------------
    Buoyancy = 0.0
    
    col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
                
    dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
    col_temp = Utils%find_str(sv%varname, Globals%Var%temp, .false.)
    col_sal = Utils%find_str(sv%varname, Globals%Var%sal, .false.)
    tag = 'radius'
    rIdx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'density'
    rhoIdx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'volume'
    volIdx = Utils%find_str(sv%varName, tag, .true.)
    tag = 'area'
    areaIdx = Utils%find_str(sv%varName, tag, .true.)
    signZ = -1.0
    
    if ((col_temp /= MV) .and. (col_sal /= MV)) then
        
        fDensity = seaWaterDensity(sv%state(:,col_sal), sv%state(:,col_temp),sv%state(:,3))
        kVisco = absoluteSeaWaterViscosity(sv%state(:,col_sal), sv%state(:,col_temp))
        
        kViscoRelation = abs(1.-(kvisco/Globals%Constants%MeanKvisco))
        where(kViscoRelation >= 0.9)
            kVisco = Globals%Constants%MeanKvisco
            fDensity = Globals%Constants%MeanDensity
        endwhere
        
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
        reynoldsNumber = self%Reynolds(sv%state(:,6), kvisco, sv%state(:,rIdx)*2) 
        cd = self%dragCoefficient(shapeFactor, sv%state(:,rIdx), reynoldsNumber)
        ! Get buoyancy
        where ((reynoldsNumber /=0.) .and. (dist2bottom > landIntThreshold))
            Buoyancy(:,3) = signZ*sqrt((-2.*Globals%Constants%Gravity%z) * (shapeFactor/cd) * densityRelation)
        endwhere
        
    else
        !If there is no salt and temperatue Compute buoyancy using constant density and temp
        ! and standardad terminal velocity
        ! Compute buoyancy using state equation for temperature and viscosity

        kVisco = Globals%Constants%MeanKvisco
        fDensity = Globals%Constants%MeanDensity
        
        signZ = -1.0
        where(fDensity-sv%state(:,rhoIdx) >= 0)
            signZ = 1.0
        endwhere
       
        ! Boundary density values could be 0. This avoid underdumped values on density. 
        ! Just density relation of 90 % related to mean water density are allowed.
        densityRelation = abs(1.- (sv%state(:,rhoIdx)/fDensity))          
        ! Get drag and shapefactor
        shapeFactor = self%SphericalShapeFactor(sv%state(:,areaIdx),sv%state(:,volIdx))
        reynoldsNumber = self%Reynolds(sv%state(:,6), kvisco, sv%state(:,rIdx))
        cd = self%dragCoefficient(shapeFactor, sv%state(:,rIdx), reynoldsNumber) 
        ! Get buoyancy
        where ((reynoldsNumber /=0.) .and. (dist2bottom < landIntThreshold))
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
    integer :: id
    !To avoid reynolds with 0 value (which will give infinit value in DragCoefficient 
    Reynolds = max(abs(characteristicVelocity), 1.0E-8)*characteristicLength/kvisco
    !Reynolds = abs(characteristicVelocity)*characteristicLength/kvisco

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
    integer :: np, nf, col_bat, bkg
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
                !correcting for maximum admissible level in the background (only if option AddBottomCell is not used)
                ! this is to allow the particles to reach the bathymetry, instead of stoping at the layer centre depth.
                maxLevel = bdata(bkg)%getDimExtents(Globals%Var%level, .false.)
                if (maxLevel(2) /= MV) where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt >= maxLevel(2)) CorrectVerticalBounds(:,3) =((maxLevel(2)-sv%state(:,3))/dt)*0.99
                !interpolating all of the data
                call self%Interpolator%run(sv%state, bdata(bkg), time, var_dt, var_name, requiredVars)
                !update minium vertical position
                col_bat = Utils%find_str(var_name, Globals%Var%bathymetry)
                !Dont let particles drop below bathymetric value. in this case, vertical displacement is forced to make the particle reach the bathymetric value*0.99
                where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt < var_dt(:,col_bat)) CorrectVerticalBounds(:,3) = ((var_dt(:,col_bat)-sv%state(:,3))/dt)*0.99
                deallocate(var_name)
                deallocate(var_dt)
            end if
        end if
    end do

    end function CorrectVerticalBounds
    
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Resuspend particles based on shear erosion calculated from currents and waves. 
    !> @param[in] self, sv, bdata, time, dt
    !---------------------------------------------------------------------------
    function Resuspension(self, sv, bdata, time, dt)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Resuspension
    integer :: col_temp, col_sal, col_dwz, col_bat, col_hs, col_ts, col_wd
    real(prec) :: landIntThreshold
    real(prec), dimension(:,:), allocatable :: var_dt
    real(prec), dimension(size(sv%state,1)) :: velocity_mod, water_density, tension
    real(prec), dimension(size(sv%state,1)) :: dist2bottom, z0
    type(string), dimension(:), allocatable :: var_name
    type(string), dimension(:), allocatable :: requiredVars
    real(prec) :: P = 1013.
    real(prec) :: EP = 1/3
    real(prec) :: waterKinematicVisc = 1e-6
    real(prec) :: taum, taumax, cdr, cds, rec, rew, fw, fws, fwr, as, ar, t1, t2, t3, a1, a2
    real(prec) :: cdm, cdms, cdmaxs, cdmax, cdmr, cdmaxr, ubw, abw, waveLength, celerity
    real(prec) :: coefA, coefB, c0, omeg, wavn, cPhi, cWphi, dwz, bat
    real(prec) :: U, psi, dsilt, dsand, dgravel, kscr, kscmr, kscd, fcs, ffs, kscr_max, kscmr_max, kscd_max
    real(prec) :: ksc, ks, grainroughnessfactor, d50
    real(prec) :: abs_cos_angle, abs_sin_angle, ubw_aux, aux_1, aux_2, fws1, fwr1, dlog_t1, ts
    !Begin-----------------------------------------------------------------------------------
    Resuspension = 0
    landIntThreshold = -0.98
    z0 = Globals%Constants%Rugosity
    dsilt = 32e-6
    dsand = 62e-6
    dgravel = 0.002
    d50 = 0.0005
    grainroughnessfactor = 2.5
    
    if (Globals%Constants%ResuspensionCoeff > 0.0) then
        allocate(requiredVars(4))
        requiredVars(1) = Globals%Var%hs
        requiredVars(2) = Globals%Var%ts
        requiredVars(3) = Globals%Var%wd
        requiredVars(4) = Globals%Var%temp
        
        call KernelUtils_VerticalMotion%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, reqVertInt = .false.)
        
        col_temp = Utils%find_str(var_name, Globals%Var%temp, .true.)
        col_sal = Utils%find_str(sv%varName, Globals%Var%sal, .true.)
        col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
        col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
        col_hs = Utils%find_str(var_name, Globals%Var%hs, .false.)
        col_ts = Utils%find_str(var_name, Globals%Var%ts, .false.)
        col_wd = Utils%find_str(var_name, Globals%Var%wd, .false.)

        !Start computations --------------------------------------------------------------------
        velocity_mod = 0
        Tension = 0
        water_density = 0
        dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
        where (dist2bottom < landIntThreshold) water_density = seaWaterDensity(sv%state(:,col_sal), sv%state(:,col_temp),sv%state(:,3))
        !where (dist2bottom < landIntThreshold) water_density = seaWaterDensity(varVert_dt(:,col_sal), varVert_dt(:,col_temp),sv%state(:,3))
        if ((col_hs /= MV_INT) .and. (col_ts /= MV_INT)) then
            !Found wave fields to use
            do i=1, size(sv%state,1)
                if (dist2bottom(i) < landIntThreshold) then
                    !TAUM is used for determining the friction governing the current
                    taum=0.
                    !TAUMAX is used to determine the threshold of sediment motion
                    taumax=0.
                    ubw = 0.001
                    abw = 0.0011
                    bat = -sv%state(i,col_bat)
                        
                    if ((var_dt(i,col_hs) > 0.1) .and. (bat > 5.0) .and. (var_dt(i,col_ts) > 0.1)) then
            !----------------------------Start calculation of abw and ubw-------------------------------------------
                        ts = var_dt(i,col_ts)
                        hs = var_dt(i,col_hs)
                        !coefA = Gravity * wavePeriod / (2.*PI). 1.5613113 = gravity/2.*PI
                        coefA = 1.5613113 * ts
                        !coefB = 2*Gravity / wavePeriod
                        coefB = 6.28318 / ts
                        c0 = sqrt(9.81*bat)
                        !Velocity modulus without the logaritmic profile close to the bottom
                        velocity_mod(i) = sqrt(sv%state(i,4)**2+sv%state(i,5)**2)
                            
                        celerity = wave_celerity(c0, bat, coefA, coefB)
                        waveLength = celerity * ts
                        omeg = 6.28318 / ts
                        !wavn = 2*pi / wavelength. wavelength = celerity * waveperiod
                        wavn = 6.28318 / (celerity * ts)
                        !To avoid sinh results larger then 1e10
                        ubw= 0
                        if (wavn * bat < 23.7) ubw = 0.5 * omeg * hs / sinh(wavn * bat)
                        !To avoid overflow errors
                        if (ubw < 1e-6)  ubw = 0.    
                        abw = ubw / omeg
            !----------------------------End calculation of abw and ubw-------------------------------------------
                    end if
                        
            !---------------------------Compute rugosity----------------------------------------------------------
                    !Average velocity
                    U = sqrt(sv%state(i,4)**2.0 + sv%state(i,5)**2.0)/0.4* (dlog(bat/z0(i)) - 1.0 + z0(i)/bat)
                        
                    uwc2 = U**2.0 + ubw**2.0
                    !Using sand density instead of particle density,
                    !because it is assumed that the bottom is mostly filled by sand and not detritus
                    relativedensity = 2650.0/water_density(i)
                    !Mobility parameter
                    psi = uwc2/((relativedensity-1)*9.81*d50)

                    if(d50 < dsilt) then
                        kscr = 20.0 * dsilt
                        kscmr = 0.
                        kscd = 0.
                    else
                        fcs = 1.0
                        if(d50 > 0.25*dgravel) fcs = (0.25*dgravel/d50)**1.5

                        ffs = 1.0
                        if(d50 < 1.5*dsand) ffs = d50/(1.5*dsand)

                        !Ripples
                        kscr_max = 0.075
                        kscr = min(fcs*d50*(85.0-65.0*tanh(0.015*(psi-150.0))), kscr_max)

                        !Mega-ripples
                        kscmr_max = min(0.01*bat, 0.2)
                        kscmr = min(2.0e-6*ffs*bat*(1-exp(-0.05*psi))*(550.0-psi), kscmr_max)

                        !Dunes
                        kscd_max = min(0.04*bat, 0.4)
                        kscd = min(8.0e-6*ffs*bat*(1-exp(-0.02*psi))*(600.0-psi), kscd_max)
                    endif
                    !Current-related bed rougnhness
                    ksc = (kscr**2.0 + kscmr**2.0 + kscd**2.0)**0.5
                    !Bed rougnhness
                    ks =  grainroughnessfactor*d50 + ksc
                    !z0
                    z0(i) = Ks/30.0
            !---------------------------End Compute rugosity------------------------------------------------------
                        
                    dwz = sv%state(i,col_dwz)
                    aux = z0(i)*exp(1.001)
                        
                    if (dwz < aux) dwz = aux
                        
                    cdr=(0.40/(dlog(dwz/z0(i))-1.))**2.0
                    rec=velocity_mod(i)*dwz/waterKinematicVisc
                    !rec=velocity_mod(i)*bat/waterKinematicVisc - como esta no soulsbury
                    cds = 0.
                    if(velocity_mod(i) > 1.0e-3) cds=0.0001615*exp(6.0*rec**(-0.08))
                
                    cdmax = max(cdr,cds)
                    
                    if (ubw > 1e-3) then
                        !Current angle in cartesian convention (angle between the vector and positive x-axis)
                        !57.2958279 = 180/pi
                        cPhi = atan2(sv%state(i,5), sv%state(i,4)) * 57.2958279
                        !(0, 360)
                        if(cPhi < 0.) cPhi = cPhi + 360.0
                        cWphi = var_dt(i,col_wd) - cPhi !Wave - Current angle
                
            !---------------------------------Compute drag coefficient ------------------------------------------------
                        rew=ubw*abw/waterKinematicVisc
                        fws=0.0521*rew**(-0.187)
                        fwr=1.39*(abw/z0(i))**(-0.52)
                        !fwr=1.39*(abw/Globals%Constants%Rugosity)**(-0.52)
                        fw=max(fws,fwr)
                        if(velocity_mod(i) < 1.0e-3)then !wave-only flow
                            cdmax=fw
                        else !combined wave and current flow
                            !turbulent flow
                            !Rough-turbulent wave-plus-current shear-stress
                            abs_cos_angle = abs(cos(cWphi*0.01745328))
                            abs_sin_angle = abs(sin(cWphi*0.01745328))
                            ubw_aux = ubw/velocity_mod(i)
                            fwr1 = fwr*0.5 
                            aux_1 = ubw_aux*(fwr1)**0.5
                            fws1 = fws*0.5
                            aux_2 = ubw_aux*(fws1)**0.5
                            cds1 = cds**2.0
                            ar=0.24
                            t1=max(ar*(fwr1)**0.5*(abw/z0(i)),12.0)
                            t2=bat/(t1*z0(i))
                            t3=(cdr**2+(fwr1)**2*(ubw_aux)**4)**(0.25)
                            dlog_t1 = dlog(t1)
                            a1 = max(t3*(dlog(t2)-1)/(2*dlog_t1), 0.0)
                            a2 = max(0.40*t3/dlog_t1, 0.0)
                            cdmr=((a1**2+a2)**0.5-a1)**2.0
                            !cdmaxr=((cdmr+t3*ubw/velocity_mod(i)*(fwr/2)**0.5*abs(cos(cWphi*pi/180.)))**2+(t3*ubw/velocity_mod(i)*(fwr/2)**0.5*abs(sin(cWphi*pi/180.)))**2)**0.5
                            cdmaxr=((cdmr+t3*aux_1*abs_cos_angle)**2+(t3*aux_1*abs_sin_angle)**2)**0.5
                
                            !Smooth-turbulent wave-plus-current shear-stress
                            as=0.24
                            t1=9*as*rew*(fws1)**0.5*(cds1*(velocity_mod(i)/ubw)**4.0+(fws1)**2.0)**(0.25)
                            t2=(rec/rew)*(ubw_aux)*1.0/as*(2.0/fws)**0.5
                            t3=(cds1+(fws1)**2*(ubw_aux)**4)**(0.25)
                            dlog_t1 = dlog(t1)
                            a1 = max(t3*(dlog(t2)-1)/(2.0*dlog_t1), 0.0)
                            a2 = max(0.40*t3/dlog_t1, 0.0)
                            cdms=((a1**2.0+a2)**0.5-a1)**2.0
                            !cdmaxs=((cdmr+t3*ubw/velocity_mod(i)*(fws/2)**0.5*abs(cos(cWphi*pi/180.)))**2+(t3*ubw/velocity_mod(i)*(fws/2)**0.5*abs(sin(cWphi*pi/180.)))**2)**0.5
                            cdmaxs=((cdms+t3*aux_2*abs_cos_angle)**2+(t3*aux_2*abs_sin_angle)**2)**0.5
                            if(cdmaxr > cdmaxs)then !flow is rough turbulent
                                cdmax=cdmaxr
                            else !flow is smooth turbulent
                                cdmax=cdmaxs
                            endif
                        endif
            !------------------------------------End Compute drag coefficient ------------------------------------------------
                        if(velocity_mod(i) < 1.0e-3)then !wave-only flow
                            tension(i) = 0.5*water_density(i)*fw*ubw**2
                        else !combined wave and current flow
                            tension(i) = water_density(i)*cdmax*velocity_mod(i)**2
                        endif
                
                    else !Ubw==0.
                        tension(i)=water_density(i)*cdmax*velocity_mod(i)**2
                    endif
                end if
            end do
        else
            !Make calculations where tracer is very close to the bathymetric value
            !Using equations from the MOHIDWater interface_sediment_water module
            where (dist2bottom < landIntThreshold)
                velocity_mod = sqrt(sv%state(:,4)**2.0 + sv%state(:,5)**2.0)
                tension = velocity_mod * water_density
            end where
        end if
        where ((dist2bottom < landIntThreshold) .and. Tension>Globals%Constants%Critical_Shear_Erosion)
            !Tracer gets positive vertical velocity which corresponds to a percentage of the velocity modulus
            !Resuspension(:,3) = Globals%Constants%ResuspensionCoeff * velocity_mod
            !tracers gets brought up to 0.5m
            Resuspension(:,3) = 0.5/dt
        end where
        deallocate(var_name)
        deallocate(var_dt)
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
    integer :: id
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
    real(prec), dimension(:),intent(in) :: S,T,P
    !real(prec),intent(in) :: P
    real(prec), dimension(size(S)):: AW,BW,KW,SR,A,K0_local,B
    real(prec),dimension(size(S)):: secantBulkModulus
    integer :: id
    !--- Pure water terms ---
    AW = h0 + (h1 + (h2 + h3*T)*T)*T
    !--- seawater, P = 0 ---
    BW = k0 + (k1 + k2*T)*T
    KW = e0 + (e1 + (e2 + (e3 + e4*T)*T)*T)*T
    SR = S**0.5
    A  = AW + (i0 + (i1 + i2*T)*T + j0*SR)*S
    K0_local = KW + (f0 + (f1 + (f2 + f3*T)*T)*T + (g0 + (g1 + g2*T)*T)*SR)*S
    ! --- General expression ---
    B = BW + (m0 + (m1 + m2*T)*T)*S
    secantBulkModulus = K0_local + (A + B*P)*P
    end function secantBulkModulus

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Method that returns secant bulk modulus
    !---------------------------------------------------------------------------
    function seaWaterDensity(S, T, P)
    real(prec), dimension(:),intent(in) :: S,T,P
    real(prec),dimension(size(S)) :: seaWaterDensity, Pbar! Pbar enters as meters
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
    Pbar = -0.1*P !# Convert to bar
    seaWaterDensity= seawaterDensityZeroPressure(S,T)/(1 - Pbar/secantBulkModulus(S,T,Pbar))
    
    !where (seaWaterDensity /= seaWaterDensity)
    !    seaWaterDensity = Globals%Constants%MeanDensity
    !end where
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
    !where (absoluteSeaWaterViscosity /= absoluteSeaWaterViscosity)
    !    absoluteSeaWaterViscosity = Globals%Constants%MeanKvisco
    !end where
    end function absoluteSeaWaterViscosity

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - colab atlantic
    !> @brief
    !> Method that returns wave celerity
    !---------------------------------------------------------------------------
    function wave_celerity(c0, watercolumn, coefa, coefb)
    real(prec), intent(in) :: c0, watercolumn,coefa, coefb
    real(prec) :: wave_celerity
    real(prec) :: x1, x2, f2, f1, csi
    integer :: maxiterations = 100
    integer :: it
    real (prec) :: eps = 5.0e-5
    !Begin -----------------------------------------------------
    wave_celerity = c0
    x2 = wave_celerity + 0.01
    f2 = x2 - coefa * tanh(coefb * watercolumn / x2)
    do it = 1, maxiterations
        f1 = wave_celerity - coefa * tanh(coefb * watercolumn / wave_celerity)
        ! ---> convergence test
        if (abs(f1) < eps) return
        x1 = wave_celerity
        if (abs(f1) < abs(f2)) then
            csi= f1/f2
            wave_celerity = x1 + csi/(csi - 1.0)*(x2 - x1)
        else
            csi = f2/f1
            wave_celerity = x1 + (x2 - x1)/(1.0 - csi)
        endif
        f2 = f1
        x2 = x1
    end do
    end function wave_celerity
    
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
        !Begin -----------------------------------------------------------------------
        
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