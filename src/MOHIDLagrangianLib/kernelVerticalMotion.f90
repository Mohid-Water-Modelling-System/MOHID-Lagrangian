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
	procedure :: LandIntThresholdValue								   
    procedure :: Buoyancy
    procedure :: CorrectVerticalBounds
    procedure :: Reynolds
    procedure :: DragCoefficient
    procedure :: SphericalShapeFactor
    procedure :: Divergence
    procedure :: Resuspension
	procedure, nopass :: seaWaterDensity
	procedure, nopass :: absoluteSeaWaterViscosity									 
    end type kernelVerticalMotion_class
    
    type(kernelUtils_class) :: KernelUtils_VerticalMotion   !< kernel utils

    public :: kernelVerticalMotion_class

    contains
	
	!---------------------------------------------------------------------------
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.09.12 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Evaluate the LandIntThresholdValue value  for a given Threshold_value (land threshold value for the dist2bottom function).
    !> Threshold_value is the regarding the position of the particle form seabed in [meter]
    !> @param[in] self, sv, bdata, time, Threshold_value
    !---------------------------------------------------------------------------
    function LandIntThresholdValue(self, sv, bdata, time, Threshold_value)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: part_idx, col_dist2bottom
    !integer :: nf_w, col_dist2bottom, part_idx
	
    real(prec), dimension(size(sv%state,1)) :: dist2bottom
    real(prec) :: threshold_bot_wat
	real(prec), dimension(size(sv%state,1)) :: Threshold_value	!distance from the bottom (seabed) in unit [meter]. It could be a constant * Globals%Constants%Rugosity
	real(prec), dimension(size(sv%state,1)) :: LandIntThresholdValue
    type(string) :: tag
	
    integer :: i,counterr
    !-------------------------------------------------------------------------------------
    !write(*,*)"Entrada kinematic"
    col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
	tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)
    dist2bottom = sv%state(:,col_dist2bottom)
	
    tag = 'particulate'
    part_idx = Utils%find_str(sv%varName, tag, .true.)
	
    ! threshold_bot_wat: equal to 0.0 will be correspond to the distance to the bottom equal to the last measurmet(given data).
    ! LandIntThresholdValue: will be correspond to the distance to the bottom which the volcity could be considered zero
	! LandIntThresholdValue gives the value the dist2bottom for this threshold.
	! u, v and w velocities for distance to the seabed(h(i) - bathymetry) < Rugosity reach to zero.    

	LandIntThresholdValue = -1 +   (Threshold_value)/(sv%state(:,col_dwz))
		
	end function  LandIntThresholdValue

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani CRETUS - GFNL- 2025.08.12 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the vertical velocity due to buoyancy of the tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> Based on papers by: 
	!> [1] Jalón-Rojas, I., Wang, X. H., & Fredj, E. (2019).3D numerical model to Track Marine Plastic Debris (TrackMPD)	Marine pollution bulletin, 141, 256-272.
	!> [2] Zhiyao, S., Tingting, W., Fumin, X., & Ruijie, L. (2008). A simple formula for predicting settling velocity of sediment particles. Water Science and Engineering, 1(1), 37-43.
	!> It should be noticed that in this paper the velocity is considered for a spherical particle with the diameter between 0.5 to 5 [millimeter]
    !---------------------------------------------------------------------------
    function Buoyancy(self, sv, bdata, time)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(in) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time
    integer :: rIdx, rhoIdx, areaIdx, volIdx
    integer :: col_temp, col_sal, col_dist2bottom, col_rugosityVar_sv, counterr
    !integer :: col_temp, col_sal, col_dwz, col_bat
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Buoyancy
    real(prec), dimension(size(sv%state,1)) :: fDensity, kVisco, dist2bottom
    real(prec), dimension(size(sv%state,1)) :: signZ, dimlessDiameter, densityRelation, cd,Re, settlingVelocity
    real(prec), dimension(size(sv%state,1)) :: ReynoldsNumber, kViscoRelation, A
    real(prec), dimension(2) :: maxLevel
    real(prec) :: landIntThreshold = -1.0
    type(string) :: tag
	! Logical mask_time for deactive calculation at initial time
	logical :: mask_time(size(sv%state,1))
	real(prec), dimension(size(sv%state,1)) :: Threshold_value
	real(prec), dimension(size(sv%state,1)) :: LandIntThreshold_value

    ! Begin -------------------------------------------------------------------------
    Buoyancy = 0.0
    tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)

    !col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    !col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)   
	
    !dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
    dist2bottom = sv%state(:,col_dist2bottom)
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
	tag = 'age'
    ageIdx = Utils%find_str(sv%varName, tag, .true.)
	col_rugosityVar_sv = Utils%find_str(sv%varName, Globals%Var%rugosityVar, .true.)
	
	! Threshold_value: distance from the bottom (seabed) in unit [meter].
	! It could be a constant * Globals%Constants%Rugosity
	Threshold_value = sv%state(:,col_rugosityVar_sv)
	LandIntThreshold_value =self%LandIntThresholdValue(sv, bdata, time, Threshold_value)
		
	! Initialize to .FALSE.
	mask_time = .false.
	! Apply condition only on the age column
	mask_time = (sv%state(:, ageIdx) /= 0.0)
  
    if ((col_temp /= MV_INT) .and. (col_sal /= MV_INT)) then
		        
		fDensity = seaWaterDensity(sv%state(:,col_sal), sv%state(:,col_temp),sv%state(:,3))
		kVisco = absoluteSeaWaterViscosity(sv%state(:,col_sal), sv%state(:,col_temp)) / fDensity
		
		kViscoRelation = abs(1.-(kvisco/Globals%Constants%MeanKvisco))
		where(kViscoRelation >= 0.9)
			kVisco = Globals%Constants%MeanKvisco
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

		! Get dimlessDiameter(dimensionless diameter) and settlingVelocity
		dimlessDiameter = (2.0 * sv%state(:,rIdx))**3.0 * ( -(Globals%Constants%Gravity%z) * densityRelation / kvisco**2.0 )
		settlingVelocity = (kvisco / (2.0 * sv%state(:,rIdx)) ) * dimlessDiameter * ( 38.1 + 0.93 * dimlessDiameter ** (4.0/7.0)) ** (-7.0/8.0)
		! Get buoyancy
		where ((dist2bottom > LandIntThreshold_value).and. mask_time)
			Buoyancy(:,3) = signZ*settlingVelocity
		endwhere
		!print*,'function:Buoyancy WiSalt in /MOHID-Lagrangian/src/MOHIDLagrangianLib'	
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
		
		! Get dimlessDiameter(dimensionless diameter) and settlingVelocity
		dimlessDiameter = (2.0 * sv%state(:,rIdx))**3.0 * ( -(Globals%Constants%Gravity%z) * densityRelation / kvisco**2.0 )
		settlingVelocity = (kvisco / (2.0 * sv%state(:,rIdx)) ) * dimlessDiameter * ( 38.1 + 0.93 * dimlessDiameter ** (4.0/7.0)) ** (-7.0/8.0)
!        ! Get buoyancy
		where ((dist2bottom > LandIntThreshold_value).and. mask_time)
			Buoyancy(:,3) = signZ*settlingVelocity
		end where
		!print*,'function:Buoyancy NoSalt in /MOHID-Lagrangian/src/MOHIDLagrangianLib'	
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
    !> https://www.researchgate.net/publication/244155878_Drag_Coefficient_and_Terminal_Velocity_of_Spherical_and_Non-Spherical_Particles
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
	!> @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Corrects vertical position of the tracers according to data limits
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function CorrectVerticalBounds(self, sv, svDt, bdata, dt)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: svDt
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: CorrectVerticalBounds
    integer :: col_bat, col_rugosityVar_sv
	real(prec), dimension(2) :: maxLevel
  
	logical, dimension(:), allocatable :: mask_0, mask_1, mask_2
	allocate(mask_0(size(sv%state, 1)))
	allocate(mask_1(size(sv%state, 1)))
	allocate(mask_2(size(sv%state, 1)))	
	
    CorrectVerticalBounds = svDt
	
    ! if vertical dvdt is 0, don't correct
    if (all(CorrectVerticalBounds(:,3) == 0.)) then
		!print*, 'All CorrectVerticalBounds(:,3) == 0 '
        return
    end if
    
    maxLevel = bdata(1)%getDimExtents(Globals%Var%level, .false.)
    if (maxLevel(2) /= MV) then 
		where (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt >= maxLevel(2))
			CorrectVerticalBounds(:,3) =((maxLevel(2)-sv%state(:,3))/dt)*0.9999
		end where
	end if
    
    col_bat = Utils%find_str(sv%varname, Globals%Var%bathymetry)
	col_rugosityVar_sv = Utils%find_str(sv%varName, Globals%Var%rugosityVar, .true.)	

	mask_1 = (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt < sv%state(:,col_bat) + sv%state(:,col_rugosityVar_sv)) .and. (sv%state(:,col_bat) <= -sv%state(:,col_rugosityVar_sv))
	mask_2 = (sv%state(:,3) + CorrectVerticalBounds(:,3)*dt < sv%state(:,col_bat) + sv%state(:,col_rugosityVar_sv)) .and. (sv%state(:,col_bat) >  -sv%state(:,col_rugosityVar_sv))
	mask_0 = (mask_1 .or. mask_2)
	
	! Apply only where `mask_1` is true: when z < bathymetry and bathymetry < 0
	where (mask_1)
		CorrectVerticalBounds(:,3) = ((sv%state(:,col_bat) + sv%state(:,col_rugosityVar_sv) - sv%state(:,3)) / dt) * 1.0
!		sv%state(:,4) = 0.0
!		sv%state(:,5) = 0.0
!		sv%state(:,6) = 0.0
	end where

	! Apply only where `mask_2` is true: when z < bathemetry and bathymetry is zero or less than the Rugosity!!
	where (mask_2)
		CorrectVerticalBounds(:,3) = (((0.0) - sv%state(:,3))/dt)*1.0
		sv%state(:,4) = 0.0
		sv%state(:,5) = 0.0
		sv%state(:,6) = 0.0
		sv%landIntMask = Globals%Mask%landVal
	end where
	
!	where (mask_0 .and. svDt(:,3) /= 0.)
!		CorrectVerticalBounds(:,1) = (CorrectVerticalBounds(:,3) / svDt(:,3)) *svDt(:,1) *0.999
!		CorrectVerticalBounds(:,2) = (CorrectVerticalBounds(:,3) / svDt(:,3)) *svDt(:,2) *0.999
!	end where	

!	where (mask_1 .and. sv%state(:,col_bat) >= Globals%Constants%BeachingLevel)
!	  sv%landIntMask = + Globals%Mask%landVal
!	elsewhere (mask_1)
!	  sv%landIntMask = - Globals%Mask%landVal
!	end where

	deallocate(mask_0)
	deallocate(mask_1)
	deallocate(mask_2)

    end function CorrectVerticalBounds
	
	!---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic 
	!> Modified @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Resuspend particles based on shear erosion calculated from currents and waves. 
    !> @param[in] self, sv, bdata, time, dt
	!---------------------------------------------------------------------------																				
    function Resuspension(self, sv, bdata, time, dt)
    class(kernelVerticalMotion_class), intent(inout) :: self
    type(stateVector_class), intent(inout) :: sv
    type(background_class), dimension(:), intent(in) :: bdata
    real(prec), intent(in) :: time, dt
    real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: Resuspension
	real(prec), dimension(size(sv%state,1)) :: ResuspensionRandmValue, ResuspensionProb_ResidenceTime
    integer :: col_age, col_temp, col_sal, col_dwz, col_bat, col_hs, col_ts, col_wd, col_dist2bottom, i
    integer :: col_rugosityVar_sv, col_D50Var_sv
	integer :: rhoIdx, rIdx
	integer :: col_u, col_v, col_w
    real(prec) :: landIntThreshold
	real(prec) :: tension_weight_coefficient
	real(prec), dimension(size(sv%state,1)) :: Shear2MShear1, Shear2PShear1, ResuspensionCrShear1, ResuspensionCrShear2
	real(prec), dimension(size(sv%state,1)) :: ResuspensionCrShield1, ResuspensionCrShield2
	real(prec), dimension(size(sv%state,1)) :: ResuspensionCrShieldSediment1, ResuspensionCrShieldSediment2
    real(prec), dimension(:,:), allocatable :: var_hor_dt
    real(prec), dimension(size(sv%state,1)) :: velocity_mod, velocity_mod_p2
    real(prec), dimension(size(sv%state,1)) :: water_density, kVisco, tension, tension_weight
    real(prec), dimension(size(sv%state,1)) :: dist2bottom, z0
    real(prec), dimension(size(sv%state,1)) :: densityRelation, densityRelationTracer, densityRelationSediment
    real(prec), dimension(size(sv%state,1)) :: kViscoRelation
    real(prec), dimension(size(sv%state,1)) :: dimlessDiameter, dimlessDiameterTracer, dimlessDiameterSediment
    type(string), dimension(:), allocatable :: var_hor_name
    type(string), dimension(:), allocatable :: requiredHorVars
    real(prec) :: P = 1013.
    real(prec) :: EP = 1/3
    real(prec) :: waterKinematicVisc = 1e-6
    real(prec) :: taum, taumax, cdr, cds, rec, rew, fw, fws, fwr, as, ar, t1, t2, t3, a1, a2
    real(prec) :: cds1, aux, relativedensity, uwc2, hs
    real(prec) :: cdm, cdms, cdmaxs, cdmax, cdmr, cdmaxr, ubw, abw, waveLength, celerity
    real(prec) :: coefA, coefB, c0, omeg, wavn, cPhi, cWphi, dwz, bat
    real(prec) :: U, psi, dsilt, dsand, dgravel, kscr, kscmr, kscd, fcs, ffs, kscr_max, kscmr_max, kscd_max
    real(prec) :: ksc, ks, grainroughnessfactor, d50
    real(prec) :: abs_cos_angle, abs_sin_angle, ubw_aux, aux_1, aux_2, fws1, fwr1, dlog_t1, ts
    type(string) :: tag
	real(prec), dimension(size(sv%state,1)) :: Threshold_value
	real(prec), dimension(size(sv%state,1)) :: LandIntThreshold_value
	real(prec), dimension(size(sv%state,1)) :: U_asterisk, V_asterisk, W_asterisk
	real(prec) :: VonKarman = 0.4
    real(prec) :: threshold_bot_wat	
	real(8), dimension(size(sv%state,1)) :: aux_r8, aux_r9
    real(prec) 		:: densitySediment	
	integer			:: counterr
    !Begin-----------------------------------------------------------------------------------
	
	densitySediment = Globals%Constants%D50Density
	
    tag = 'density'
    rhoIdx = Utils%find_str(sv%varName, tag, .true.)	

    tag = 'radius'
    rIdx = Utils%find_str(sv%varName, tag, .true.)
	
    tag = 'age'
    col_age = Utils%find_str(sv%varName, tag, .true.)	
	
    tag = 'dist2bottom'
    col_dist2bottom = Utils%find_str(sv%varName, tag, .true.)
    dist2bottom = sv%state(:,col_dist2bottom)
	col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
    threshold_bot_wat = (Globals%Mask%waterVal + Globals%Mask%bedVal) * 0.0

	col_rugosityVar_sv = Utils%find_str(sv%varName, Globals%Var%rugosityVar, .true.)	
	col_D50Var_sv = Utils%find_str(sv%varName, Globals%Var%D50Var, .true.)
	
	U_asterisk = 0.0
	V_asterisk = 0.0
	W_asterisk = 0.0
		
    Resuspension = 0

	! Threshold_value: distance from the bottom (seabed) in unit [meter]. It could be a constant * Globals%Constants%Rugosity
	Threshold_value = Globals%Constants%BedLoadThickness
	LandIntThreshold_value =self%LandIntThresholdValue(sv, bdata, time, Threshold_value)

    landIntThreshold = -0.98
    z0 = sv%state(:,col_rugosityVar_sv)
    dsilt = 32e-6
    dsand = 62e-6
    dgravel = 0.002
    d50 = 0.0005
    grainroughnessfactor = 2.5
    
    if (Globals%Constants%ResuspensionProb > 0.0) then

		! Create  random values between [0,1) to compare with Resuspension probability
		call random_number(ResuspensionRandmValue) 

        allocate(requiredHorVars(6))
        requiredHorVars(1) = Globals%Var%hs
        requiredHorVars(2) = Globals%Var%ts
        requiredHorVars(3) = Globals%Var%wd
		requiredHorVars(4) = Globals%Var%u
        requiredHorVars(5) = Globals%Var%v
        requiredHorVars(6) = Globals%Var%w
        
        call KernelUtils_VerticalMotion%getInterpolatedFields(sv, bdata, time, requiredHorVars, var_hor_dt, var_hor_name, reqVertInt = .false.)
        
		col_temp = Utils%find_str(sv%varName, Globals%Var%temp, .true.)
        col_sal = Utils%find_str(sv%varName, Globals%Var%sal, .true.)
        col_dwz = Utils%find_str(sv%varName, Globals%Var%dwz, .true.)
        col_bat = Utils%find_str(sv%varName, Globals%Var%bathymetry, .true.)
        col_hs = Utils%find_str(var_hor_name, Globals%Var%hs, .false.)
        col_ts = Utils%find_str(var_hor_name, Globals%Var%ts, .false.)
        col_wd = Utils%find_str(var_hor_name, Globals%Var%wd, .false.)
		col_u = Utils%find_str(var_hor_name, Globals%Var%u, .true.)
		col_v = Utils%find_str(var_hor_name, Globals%Var%v, .true.)
		col_w = Utils%find_str(var_hor_name, Globals%Var%w, .false.)

        !Start computations --------------------------------------------------------------------
        velocity_mod = 0
        velocity_mod_p2 = 0		
        Tension = 0
		tension_weight = 0
        water_density = 0
		kVisco = 0
		ResuspensionCrShear1 = 0
		ResuspensionCrShear2 = 0
		ResuspensionCrShield1 = 0
		ResuspensionCrShield2 = 0
		ResuspensionCrShieldSediment1 = 0
		ResuspensionCrShieldSediment2 = 0
		
        !dist2bottom = Globals%Mask%bedVal + (sv%state(:,3) - sv%state(:,col_bat)) / (sv%state(:,col_dwz))
        dist2bottom = sv%state(:,col_dist2bottom)
        !where (dist2bottom < landIntThreshold) water_density = seaWaterDensity(sv%state(:,col_sal), sv%state(:,col_temp),sv%state(:,3))
		if ((col_temp /= MV_INT) .and. (col_sal /= MV_INT)) then
			where ((dist2bottom < threshold_bot_wat))
				water_density = seaWaterDensity(sv%state(:,col_sal), sv%state(:,col_temp),sv%state(:,3))
				kVisco = absoluteSeaWaterViscosity(sv%state(:,col_sal), sv%state(:,col_temp)) / water_density
				
				kViscoRelation = abs(1.-(kvisco/Globals%Constants%MeanKvisco))
				where(kViscoRelation >= 0.9)
					kVisco = Globals%Constants%MeanKvisco
				endwhere

				! Boundary density values could be 0. This avoid underdumped values on density. 
				! Just density relation of 90 % related to mean water density are allowed.
				densityRelation = abs(1.- (sv%state(:,rhoIdx)/water_density))
				where (densityRelation >= 0.9)
					densityRelation = abs(1.- (sv%state(:,rhoIdx)/Globals%Constants%MeanDensity))
				endwhere
			end where
		else
			!If there is no salt and temperatue Compute buoyancy using constant density and temp
			where ((dist2bottom < threshold_bot_wat))
				kVisco = Globals%Constants%MeanKvisco
				water_density = Globals%Constants%MeanDensity
				! Boundary density values could be 0. This avoid underdumped values on density. 
				! Just density relation of 90 % related to mean water density are allowed.
				densityRelation = abs(1.- (sv%state(:,rhoIdx)/water_density))   
				densityRelationTracer 	= densityRelation 
				densityRelationSediment = abs(1.- (densitySediment /water_density))       				
			end where
		end if

		if ((col_hs /= MV_INT) .and. (col_ts /= MV_INT)) then
            !Found wave fields to use
            do i=1, size(sv%state,1)
                if (dist2bottom(i) < LandIntThreshold_value(i)) then
                    !TAUM is used for determining the friction governing the current
                    taum=0.
                    !TAUMAX is used to determine the threshold of sediment motion
                    taumax=0.
                    ubw = 0.001
                    abw = 0.0011
                    bat = -sv%state(i,col_bat)
                        
                    if ((var_hor_dt(i,col_hs) > 0.1) .and. (bat > 5.0) .and. (var_hor_dt(i,col_ts) > 0.1)) then
            !----------------------------Start calculation of abw and ubw-------------------------------------------
                        ts = var_hor_dt(i,col_ts)
                        hs = var_hor_dt(i,col_hs)
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
					! Hint: This average velocity needs modification!!!!!!!!
                    !Average velocity
                    !U = sqrt(sv%state(i,4)**2.0 + sv%state(i,5)**2.0)/0.4* (dlog(bat/z0(i)) - 1.0 + z0(i)/bat)
					!
					! In the BedLoadThickness zone: Corresponds to Case 2 of the velocity profile,
					! which uses an "average velocity" from the log-law for all particles in this layer.
					! Hence, we already have the average velocity (u,v,w) values in this zone at sv%sate(:,4:6) 
					U = sqrt(sv%state(i,4)**2.0 + sv%state(i,5)**2.0 + sv%state(i,6)**2.0)
					
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
                        cWphi = var_hor_dt(i,col_wd) - cPhi !Wave - Current angle
                
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
            !Using log law to calculate the (U,V,W)_asterisk, and later calculated wall tension 
		    
			where (dist2bottom < threshold_bot_wat) 
				aux_r8 = max(sv%state(:,col_dwz), 0.05) / sv%state(:,col_rugosityVar_sv)
				U_asterisk = var_hor_dt(:,col_u) * VonKarman /  dlog(aux_r8)
				V_asterisk = var_hor_dt(:,col_v) * VonKarman /  dlog(aux_r8)
				W_asterisk = var_hor_dt(:,col_w) * VonKarman /  dlog(aux_r8)
            end where
			
			! velocity_mod_p2 = velocity_mod ** 2.0
			velocity_mod_p2 = (U_asterisk ** 2.0 +  V_asterisk ** 2.0  + W_asterisk ** 2.0)
			tension =  water_density * velocity_mod_p2
				
        end if
		
		ResuspensionProb_ResidenceTime = Globals%Constants%ResuspensionProb * exp(- sv%state(:,col_age) / Globals%Constants%ResuspensionResidenceTime)
		
		! In the BedLoadThickness zone: Corresponds to Case 2 of the velocity profile, which uses an average velocity from the log-law for all particles in this layer.
		! It can be divided into 3 cases:
		!   Case 1: If 								tension < ResuspensionCriticalShear1, then particles do not move.
		!   Case 2: If ResuspensionCriticalShear1 < tension < ResuspensionCriticalShear2, then particles  move due to the drag force.
		!           In this case, we would need to solve the momentum equation, which is not cost-effective.
		!           Hence, we approximate the velocity with tension_weight * velocity in this zone, 
		!           where tension_weight is in the interval [0, 1].
		!   Case 3: If ResuspensionCriticalShear2 < tension								, then particles move according to the current velocities.
		!
		! To cover all cases in a single formulation, we can consider the following symmetrical sigmoidal:

		tension_weight_coefficient = +10.0

		if (Globals%SimDefs%ResuspensionCriticalShearMethod == 1) then
			
			where (dist2bottom < LandIntThreshold_value)
				ResuspensionCrShear1 = Globals%Constants%ResuspensionCriticalShear1
				ResuspensionCrShear2 = Globals%Constants%ResuspensionCriticalShear2
			end where
			
		else if (Globals%SimDefs%ResuspensionCriticalShearMethod == 2) then	

			where (dist2bottom < LandIntThreshold_value) dimlessDiameterTracer = (2.0 * sv%state(:,rIdx))**3.0 * ( -(Globals%Constants%Gravity%z) * densityRelationTracer / kvisco**2.0 )			

			where (dist2bottom < LandIntThreshold_value)
				ResuspensionCrShield1 = ((0.30 / (1.0 + 1.2 * dimlessDiameterTracer)) + 0.055 * (1.0 - exp(-0.02 * dimlessDiameterTracer))) 
				ResuspensionCrShield2 = ((0.30 / (1.0 + 1.0 * dimlessDiameterTracer)) + 0.100 * (1.0 - exp(-0.05 * dimlessDiameterTracer)))
			end where

			where (dist2bottom < LandIntThreshold_value)
				ResuspensionCrShear1 = max((ResuspensionCrShield1) * (-Globals%Constants%Gravity%z) * (sv%state(:,rhoIdx) - water_density) * (2.0 * sv%state(:,rIdx)), 0.0000010)
				ResuspensionCrShear2 = max((ResuspensionCrShield2) * (-Globals%Constants%Gravity%z) * (sv%state(:,rhoIdx) - water_density) * (2.0 * sv%state(:,rIdx)), 0.0000011)
			end where

		else if (Globals%SimDefs%ResuspensionCriticalShearMethod == 3) then	

			where (dist2bottom < LandIntThreshold_value) dimlessDiameterSediment = (2.0 * sv%state(:,col_D50Var_sv))**3.0 * ( -(Globals%Constants%Gravity%z) * densityRelationSediment / kvisco**2.0 )			

			where (dist2bottom < LandIntThreshold_value)
				ResuspensionCrShieldSediment1 = ((0.30 / (1.0 + 1.2 * dimlessDiameterSediment)) + 0.055 * (1.0 - exp(-0.02 * dimlessDiameterSediment))) 
				ResuspensionCrShieldSediment2 = ((0.30 / (1.0 + 1.0 * dimlessDiameterSediment)) + 0.100 * (1.0 - exp(-0.05 * dimlessDiameterSediment)))
			end where

			where (dist2bottom < LandIntThreshold_value)
				ResuspensionCrShield1 = 0.5588D0 * ResuspensionCrShieldSediment1 * ((2.D0 * sv%state(:,rIdx)) / sv%state(:,col_D50Var_sv) ) ** (-0.503D0)
				ResuspensionCrShield2 = 0.5588D0 * ResuspensionCrShieldSediment2 * ((2.D0 * sv%state(:,rIdx)) / sv%state(:,col_D50Var_sv) ) ** (-0.503D0)
			end where

			where (dist2bottom < LandIntThreshold_value)
				ResuspensionCrShear1 = max((ResuspensionCrShield1) * (-Globals%Constants%Gravity%z) * (sv%state(:,rhoIdx) - water_density) * (2.0 * sv%state(:,rIdx)), 0.0000010)
				ResuspensionCrShear2 = max((ResuspensionCrShield2) * (-Globals%Constants%Gravity%z) * (sv%state(:,rhoIdx) - water_density) * (2.0 * sv%state(:,rIdx)), 0.0000011)
			end where

		end if

		Shear2MShear1 = ResuspensionCrShear2 - ResuspensionCrShear1
		Shear2PShear1 = ResuspensionCrShear2 + ResuspensionCrShear1
		where ((dist2bottom < LandIntThreshold_value) .and. (ResuspensionRandmValue < ResuspensionProb_ResidenceTime) )
			tension_weight = 1.0 / ( 1.0 + exp( -(tension_weight_coefficient/Shear2MShear1) * (tension - 0.5 * Shear2PShear1) ) ) 
		end where

        where ((dist2bottom < LandIntThreshold_value) .and. (ResuspensionRandmValue < ResuspensionProb_ResidenceTime) )
			
			sv%state(:,4) = tension_weight(:) * sv%state(:,4)
			sv%state(:,5) = tension_weight(:) * sv%state(:,5)
			sv%state(:,6) = tension_weight(:) * sv%state(:,6)
						
			Resuspension(:,1) = Utils%m2geo(sv%state(:,4), sv%state(:,2), .false.)
			Resuspension(:,2) = Utils%m2geo(sv%state(:,5), sv%state(:,2), .true.)
			Resuspension(:,3) = sv%state(:,6) 

        end where
			
        deallocate(var_hor_name)
        deallocate(var_hor_dt)

	else
	
		where (dist2bottom < LandIntThreshold_value)
			sv%state(:,4) = 0.0
			sv%state(:,5) = 0.0
			sv%state(:,6) = 0.0
		end where

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