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
                where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt < var_dt(:,nf)) CorrectVerticalBounds(:,3) = ((var_dt(:,nf)-sv%state(:,3))/dt)*0.99
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
    integer :: col_lim, col_temp, col_sal, col_dwz, col_bat, col_hs, col_ts, col_wd, col_age
    real(prec) :: landIntThreshold
    real(prec), dimension(:,:), allocatable :: var_dt, varVert_dt, var_dt2
    real(prec), dimension(size(sv%state,1)) :: velocity_mod, water_density, tension
    real(prec), dimension(size(sv%state,1)) :: dist2bottom
    type(string), dimension(:), allocatable :: var_name, var_name_vert, var_name_2
    type(string), dimension(:), allocatable :: requiredVars, requiredVerticalVars, requiredHorVars2
    real(prec) :: P = 1013.
    real(prec) :: EP = 1/3
    real(prec) :: waterKinematicVisc = 1e-6
    real(prec) :: taum, taumax, cdr, cds, rec, rew, fw, fws, fwr, as, ar, t1, t2, t3, a1, a2
    real(prec) :: cdm, cdms, cdmaxs, cdmax, cdmr, cdmaxr, ubw, abw, waveLength, celerity
    real(prec) :: coefA, coefB, c0, omeg, wavn, cPhi, cWphi, dwz, bat, z0
    real(prec) :: abs_cos_angle, abs_sin_angle, ubw_aux, aux_1, aux_2, fws1, fwr1, dlog_t1, teste1, teste2, ts
    type(string) :: tag
    logical :: escreve
    !Begin-----------------------------------------------------------------------------------
    Resuspension = 0
    landIntThreshold = -0.98
    escreve = .false.
    z0 = Globals%Constants%Rugosity
    
    if (Globals%Constants%ResuspensionCoeff > 0.0) then
        allocate(requiredVars(3))
        requiredVars(1) = Globals%Var%hs
        requiredVars(2) = Globals%Var%ts
        requiredVars(3) = Globals%Var%wd
        
        allocate(requiredVerticalVars(4))
        requiredVerticalVars(1) = Globals%Var%dwz
        requiredVerticalVars(2) = Globals%Var%bathymetry
        requiredVerticalVars(3) = Globals%Var%temp
        requiredVerticalVars(4) = Globals%Var%sal
        
        call KernelUtils_VerticalMotion%getInterpolatedFields(sv, bdata, time, requiredVars, var_dt, var_name, reqVertInt = .false.)
        call KernelUtils_VerticalMotion%getInterpolatedFields(sv, bdata, time, requiredVerticalVars, varVert_dt, var_name_vert, reqVertInt = .true.)
        col_temp = Utils%find_str(var_name_vert, Globals%Var%temp, .true.)
        col_sal = Utils%find_str(var_name_vert, Globals%Var%sal, .true.)
        col_dwz = Utils%find_str(var_name_vert, Globals%Var%dwz, .true.)
        col_bat = Utils%find_str(var_name_vert, Globals%Var%bathymetry, .true.)
        col_hs = Utils%find_str(var_name, Globals%Var%hs, .false.)
        col_ts = Utils%find_str(var_name, Globals%Var%ts, .false.)
        col_wd = Utils%find_str(var_name, Globals%Var%wd, .false.)
        
        if (col_hs /= MV_INT .and. col_ts /= MV_INT) then
            allocate(requiredHorVars2(2))
            requiredHorVars2(1) = Globals%Var%u
            requiredHorVars2(2) = Globals%Var%v
            call KernelUtils_VerticalMotion%getInterpolatedFields(sv, bdata, time, requiredHorVars2, var_dt2, var_name_2, reqVertInt = .false.)
            col_u = Utils%find_str(var_name_2, Globals%Var%u, .true.)
            col_v = Utils%find_str(var_name_2, Globals%Var%v, .true.)
            tag = 'age'
            col_age = Utils%find_str(sv%varName, tag, .true.)
        end if 

        !Start computations --------------------------------------------------------------------
        velocity_mod = 0
        Tension = 0
        water_density = 0
        dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - varVert_dt(:,col_bat)) / (varVert_dt(:,col_bat) - varVert_dt(:,col_dwz))

        where (dist2bottom < landIntThreshold) water_density = seaWaterDensity(varVert_dt(:,col_sal), varVert_dt(:,col_temp),P)
        if ((col_hs /= MV_INT) .and. (col_ts /= MV_INT)) then
            write(*,*)"Entrei na parte das ondas"
            !Found wave fields to use
            do i=1, size(sv%state,1)
                if (mod(sv%state(i,col_age),dt*30)==0) then
                    escreve = .true.
                    if (dist2bottom(i) < landIntThreshold) then
                        !TAUM is used for determining the friction governing the current
                        taum=0.
                        !TAUMAX is used to determine the threshold of sediment motion
                        taumax=0.
                        ubw = 0.001
                        abw = 0.0011
                        bat = -varVert_dt(i,col_bat)
                        
                        if ((var_dt(i,col_hs) > 0.1) .and. (bat > 5) .and. (var_dt(i,col_ts) > 0.1)) then
                !----------------------------Start calculation of abw and ubw-------------------------------------------
                            ts = var_dt(i,col_ts)
                            hs = var_dt(i,col_hs)
                            !coefA = Gravity * wavePeriod / (2.*PI). 1.5613113 = gravity/2.*PI
                            coefA = 1.5613113 * ts
                            !coefB = 2*Gravity / wavePeriod
                            coefB = 6.28318 / ts
                            c0 = sqrt(9.81*bat)
                            !Velocity modulus without the logaritmic profile close to the bottom
                            velocity_mod(i) = sqrt(var_dt2(i,col_u)**2+var_dt2(i,col_v)**2)
                            
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
                        
                        dwz = varVert_dt(i,col_dwz)
                        aux = z0*exp(1.001)
                        
                        if (dwz < aux) dwz = aux
                        
                        cdr=(0.40/(dlog(dwz/z0)-1.))**2
                        rec=velocity_mod(i)*dwz/waterKinematicVisc
                        !rec=velocity_mod(i)*bat/waterKinematicVisc - como esta no soulsbury
                        cds = 0.
                        if(velocity_mod(i) > 1e-3) cds=0.0001615*exp(6.*rec**(-0.08))
                
                        cdmax = max(cdr,cds)
                    
                        if (ubw > 1e-3) then
                            !Current angle in cartesian convention (angle between the vector and positive x-axis)
                            !57.2958279 = 180/pi
                            cPhi = atan2(sv%state(i,5), sv%state(i,4)) * 57.2958279
                            !(0, 360)
                            if(cPhi < 0.) cPhi = cPhi + 360
                            cWphi = var_dt(i,col_wd) - cPhi !Wave - Current angle
                
                !---------------------------------Compute drag coefficient ------------------------------------------------
                            rew=ubw*abw/waterKinematicVisc
                            fws=0.0521*rew**(-0.187)
                            fwr=1.39*(abw/z0)**(-0.52)
                            !fwr=1.39*(abw/Globals%Constants%Rugosity)**(-0.52)
                            fw=max(fws,fwr)
                            if(velocity_mod(i) < 1e-3)then !wave-only flow
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
                                cds1 = cds**2
                                ar=0.24
                                t1=max(ar*(fwr1)**0.5*(abw/z0),12.0)
                                t2=bat/(t1*z0)
                                t3=(cdr**2+(fwr1)**2*(ubw_aux)**4)**(0.25)
                                dlog_t1 = dlog(t1)
                                a1 = max(t3*(dlog(t2)-1)/(2*dlog_t1), 0.0)
                                a2 = max(0.40*t3/dlog_t1, 0.0)
                                cdmr=((a1**2+a2)**0.5-a1)**2
                                !cdmaxr=((cdmr+t3*ubw/velocity_mod(i)*(fwr/2)**0.5*abs(cos(cWphi*pi/180.)))**2+(t3*ubw/velocity_mod(i)*(fwr/2)**0.5*abs(sin(cWphi*pi/180.)))**2)**0.5
                                cdmaxr=((cdmr+t3*aux_1*abs_cos_angle)**2+(t3*aux_1*abs_sin_angle)**2)**0.5
                
                                !Smooth-turbulent wave-plus-current shear-stress
                                as=0.24
                                t1=9*as*rew*(fws1)**0.5*(cds1*(velocity_mod(i)/ubw)**4+(fws1)**2)**(0.25)
                                t2=(rec/rew)*(ubw_aux)*1/as*(2/fws)**0.5
                                t3=(cds1+(fws1)**2*(ubw_aux)**4)**(0.25)
                                dlog_t1 = dlog(t1)
                                a1 = max(t3*(dlog(t2)-1)/(2*dlog_t1), 0.0)
                                a2 = max(0.40*t3/dlog_t1, 0.0)
                                cdms=((a1**2+a2)**0.5-a1)**2
                                !cdmaxs=((cdmr+t3*ubw/velocity_mod(i)*(fws/2)**0.5*abs(cos(cWphi*pi/180.)))**2+(t3*ubw/velocity_mod(i)*(fws/2)**0.5*abs(sin(cWphi*pi/180.)))**2)**0.5
                                cdmaxs=((cdms+t3*aux_2*abs_cos_angle)**2+(t3*aux_2*abs_sin_angle)**2)**0.5
                                if(cdmaxr > cdmaxs)then !flow is rough turbulent
                                    cdmax=cdmaxr
                                else !flow is smooth turbulent
                                    cdmax=cdmaxs
                                endif
                            endif
                !------------------------------------End Compute drag coefficient ------------------------------------------------
                            if(velocity_mod(i) < 1e-3)then !wave-only flow
                                tension(i) = 0.5*water_density(i)*fw*ubw**2
                                if (i == 2) then
                                    write(*,*)"hs = ", hs
                                    write(*,*)"bat = ", bat
                                    write(*,*)"ts = ", ts
                                    write(*,*)"fw = ", fw
                                    write(*,*)"ubw = ", ubw
                                    write(*,*)"tension(i) so ondas = ", tension(i)
                                end if
                            else !combined wave and current flow
                                tension(i) = water_density(i)*cdmax*velocity_mod(i)**2
                                if (i == 2) then
                                    write(*,*)"ubw = ", ubw
                                    write(*,*)"cdmax = ", cdmax
                                    write(*,*)"fw = ", fw
                                    write(*,*)"water_density = ", water_density(i)
                                     write(*,*)"velocity_mod(i) = ", velocity_mod(i)
                                    write(*,*)"hs = ", hs
                                    write(*,*)"bat = ", bat
                                    write(*,*)"ts = ", ts
                                    write(*,*)"tension(i) correntes mais ondas = ", tension(i)
                                end if
                            endif
                
                        else !Ubw==0.
                            tension(i)=water_density(i)*cdmax*velocity_mod(i)**2
                            if (i == 2) then
                                write(*,*)"hs = ", hs
                                write(*,*)"bat = ", bat
                                write(*,*)"ts = ", ts
                                write(*,*)"tension(i) ubw = 0 = ", tension(i)
                            end if
                        endif
                    end if
                end if
            end do
            deallocate(var_name_2)
            deallocate(requiredHorVars2)
        else
            !Make calculations where tracer is very close to the bathymetric value
            !Using equations from the MOHIDWater interface_sediment_water module
            where (dist2bottom < landIntThreshold)
                velocity_mod = sqrt(sv%state(:,4)**2 + sv%state(:,5)**2)
                tension = velocity_mod * water_density
            end where
        end if
        
        where ((dist2bottom < landIntThreshold) .and. Tension>Globals%Constants%Critical_Shear_Erosion)
            !Tracer gets positive vertical velocity which corresponds to a percentage of the velocity module
            Resuspension(:,3) = Globals%Constants%ResuspensionCoeff * velocity_mod
        end where 
        
        deallocate(var_name)
        deallocate(var_name_vert)
        deallocate(var_dt)
        deallocate(varVert_dt)
        if ((size(sv%state,1) > 1) .and. (escreve)) then
            write(*,*)"dist2bottom = ", dist2bottom(2)
            write(*,*)"tension = ", tension(2)
            write(*,*)"velocity_mod = ", velocity_mod(2)
        end if
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