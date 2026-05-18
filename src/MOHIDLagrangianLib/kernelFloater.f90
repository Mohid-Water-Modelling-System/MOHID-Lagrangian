    !------------------------------------------------------------------------------
    !        Colab+Atlantic, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_floater
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : Colab+Atlantic, Marine Modelling Group
    ! DATE          : March 2026
    !> @author
    !> Mohsen Shabani Email:shabani.mohsen@outlook.com	
    !
    ! DESCRIPTION:
    !> Defines an abstract physics kernel class for floater object.
    !> 	
    !------------------------------------------------------------------------------
	module kernelFloater_mod

	use common_modules
	use stateVector_mod
	use background_mod
	use interpolator_mod
	use kernelUtils_mod
	use kernelNonlinearSolverScalar_mod
	use kernelNonlinearSolverVector_mod
	use kernelNonlinearSolverVectorCoupling_mod
	use sources_mod 
	use tracerFloater_mod!, only: floater_class  
	use kernelVerticalMotion_mod
	use kernelSharedProcesses_mod
	use, intrinsic :: ieee_arithmetic
	
	implicit none

	type :: source_parameters_type
		real(prec), allocatable :: Mass(:)
		integer,    allocatable :: S_N(:)
		real(prec), allocatable :: S_rd(:), S_CD(:), S_R(:)

		integer,    allocatable :: R_N(:)
		real(prec), allocatable :: R_rd(:), R_CD(:)
		real(prec), allocatable :: R_Lx(:), R_Ly(:), R_Lz(:)
		
		integer,    allocatable :: P_N(:)
		real(prec), allocatable :: P_rd(:), P_CD(:) , P_D(:)
		real(prec), allocatable :: P_Lx(:), P_Ly(:), P_Lz(:)
	end type source_parameters_type
	

	  
	type :: kernelFloater_class
		type(interpolator_class) 			:: Interpolator
		type(kernelUtils_class) 			:: KernelUtils   !< kernel utils	
		type(kernelSharedProcesses_class) 	:: SharedProcesses   !< Shared Processes kernels
		type(kernelVerticalMotion_class) 	:: VerticalMotion   !< VerticalMotion kernels
		type(kernelNonlinearSolverVectorCoupling_class) :: NonlinearSolverVectorCoupling
	
	contains
		procedure :: initialize => initKernelFloater
		procedure :: FloaterVelocity

		procedure :: DragCoefficient_K_wind_vector
		procedure :: DragCoefficient_K_current_vector
		procedure :: DragCoefficient_K_porous_vector
		procedure :: DragCoefficient_CD_Constant_scalar
		procedure :: DragCoefficient_CD_Sphere_scalar
		procedure :: DragCoefficient_CD_Rectangular_scalar
		procedure :: DragCoefficient_CD_Porous_scalar
		procedure :: DensityKViscosity_seawater
		procedure :: DensityKViscosity_air
		procedure :: Reynolds_X
		procedure :: findSourceParameters
			
		!procedure :: momentumBalance
	end type kernelFloater_class
	
	public :: kernelFloater_class, source_parameters_type


	contains

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  velocity of the floter tracers in seawater
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
	function FloaterVelocity(self, sv, sv_atDepth, bdata, time, dt, srcPar)
	class(kernelFloater_class), intent(inout) :: self
	type(stateVector_class), intent(inout) :: sv
	type(stateVector_class), intent(in) :: sv_atDepth
	type(background_class), dimension(:), intent(in) :: bdata
	real(prec), intent(in) :: time, dt
	type(source_parameters_type), intent(in) :: srcPar
	integer :: rIdx,idx
	integer :: col_bathy	
	integer :: col_ufloater, col_vfloater, col_wfloater
	real(prec), dimension(size(sv%state,1),size(sv%state,2)) :: FloaterVelocity
	real(prec), dimension(2) :: maxLevel
	real(prec) :: landIntThreshold = -1.0_prec
	type(string) :: tag
	type(string) :: outext
	integer :: i, i_particle
	real(prec) :: x0(2), root(2)
	logical :: ok_scalar
	logical, dimension(2) :: ok_vector
	real(prec), dimension(2) :: root_vector, x0_vector, x1_vector
	real(prec) :: u0, v0
	real(prec) :: uP_D, vP_D
	logical :: hasWind
	
	real(prec), allocatable :: KViscosity_seawater_vector_atSurfc(:)
	real(prec), allocatable :: Density_seawater_vector_atSurfc(:)
	real(prec), allocatable :: KViscosity_seawater_vector_atDepth(:)
	real(prec), allocatable :: Density_seawater_vector_atDepth(:)
	real(prec), allocatable :: KViscosity_air_vector_atSurfc(:)
	real(prec), allocatable :: Density_air_vector_atSurfc(:)

	real(prec), allocatable :: WindVel(:,:)
	real(prec), allocatable :: StokesVel(:,:)
	real(prec), allocatable :: vel_floater_new(:,:)


	allocate(WindVel(size(sv%state,1),2))
	allocate(StokesVel(size(sv%state,1),2))
	allocate(vel_floater_new(size(sv%state,1),3))

		vel_floater_new		= 0.0_prec 
		FloaterVelocity		= 0.0_prec
	
		! Find the values of source paramerters[call self%findSourceParameters(sv)]
		! It is in the self, so it should call one time and later it is  shared between alll function.
		!It moved to the main kerne, so it is not calculated two times
		!call self%findSourceParameters(sv)

		call self%DensityKViscosity_seawater(sv,         KViscosity_seawater_vector_atSurfc, Density_seawater_vector_atSurfc)
		call self%DensityKViscosity_seawater(sv_atDepth, KViscosity_seawater_vector_atDepth, Density_seawater_vector_atDepth)
		call self%DensityKViscosity_air(sv,              KViscosity_air_vector_atSurfc,       Density_air_vector_atSurfc)

		WindVel = self%SharedProcesses%WindVelocity(sv, bdata, time, hasWind)
		StokesVel = self%SharedProcesses%StokesDriftVelocity(sv, bdata, time)

		tag = 'bathymetry'
		col_bathy = Utils%find_str(sv%varName, tag, .true.)
		
		tag = 'ufloater'
		col_ufloater = Utils%find_str(sv%varName, tag, .true.)
		tag = 'vfloater'
		col_vfloater = Utils%find_str(sv%varName, tag, .true.)		
		tag = 'wfloater'
		col_wfloater = Utils%find_str(sv%varName, tag, .true.)	
				
		do i_particle = 1, size(sv%state,1)

			u0 = sv%state(i_particle,4)
			v0 = sv%state(i_particle,5)
			
			uP_D = sv_atDepth%state(i_particle,4)
			vP_D = sv_atDepth%state(i_particle,5)

			if (.not. ieee_is_finite(u0) .or. .not. ieee_is_finite(v0) .or. &
			    .not. ieee_is_finite(uP_D) .or. .not. ieee_is_finite(vP_D)) then
				outext = '[kernelFloater::FloaterVelocity] Non-finite current or porous velocity detected, stopping'
				call Log%put(outext)
				stop
			end if
			if (.not. all(ieee_is_finite(WindVel(i_particle,:)))) then
				outext = '[kernelFloater::FloaterVelocity] Non-finite wind velocity detected, stopping'
				call Log%put(outext)
				stop
			end if
			if (.not. all(ieee_is_finite(StokesVel(i_particle,:)))) then
				outext = '[kernelFloater::FloaterVelocity] Non-finite Stokes drift velocity detected; setting non-finite components to zero'
				call Log%put(outext)
				!stop
				! where (.not. ieee_is_finite(StokesVel(i_particle,:)))
					! StokesVel(i_particle,:) = 0.0
				! end where
				 StokesVel(i_particle,:) = 0.0
			end if
			if (.not. ieee_is_finite(srcPar%Mass(i_particle)) .or. srcPar%Mass(i_particle) <= 0.0_prec) then
				outext = '[kernelFloater::FloaterVelocity] Invalid floater mass detected, stopping'
				call Log%put(outext)
				stop
			end if
			if (.not. ieee_is_finite(srcPar%P_D(i_particle))) then
				outext = '[kernelFloater::FloaterVelocity] Non-finite porous depth detected, stopping'
				call Log%put(outext)
				stop
			end if
					
			! x0_vector = [u0, v0]
			x0_vector = [sv%state(i_particle, col_ufloater), sv%state(i_particle, col_vfloater)]
			x1_vector = x0_vector + [1.0e-1_prec, 1.0e-1_prec]

			! call self%NonlinearSolverVectorCoupling%NewtonSolver(momentumBalanceFloater, x0_vector, root_vector, ok_vector)

			! if (.not. all(ok_vector)) then
				! call self%NonlinearSolverVectorCoupling%SecantSolver(momentumBalanceFloater, x0_vector, x1_vector, root_vector, ok_vector)
			! end if

			! if (all(ok_vector) .and. -sv%state(i_particle, col_bathy) > srcPar%P_D(i_particle)) then
				! vel_floater_new(i_particle,1) = root_vector(1)
				! vel_floater_new(i_particle,2) = root_vector(2)
				! vel_floater_new(i_particle,3) = 0.0_prec
			! elseif (.not.all(ok_vector) .and. -sv%state(i_particle, col_bathy) > srcPar%P_D(i_particle)) then
				! print*, "NOT all(ok_vector):", all(ok_vector), "but is not beached"
				! vel_floater_new(i_particle,1) = x0_vector(1)
				! vel_floater_new(i_particle,2) = x0_vector(2)
				! vel_floater_new(i_particle,3) = 0.0_prec
			! else
				! print*, "(ok_vector):", all(ok_vector),i_particle, "bath:",-sv%state(i_particle, col_bathy),"P_D:",srcPar%P_D(i_particle)
				! vel_floater_new(i_particle,1) = 0.0_prec
				! vel_floater_new(i_particle,2) = 0.0_prec
				! vel_floater_new(i_particle,3) = 0.0_prec
			! end if

			if (-sv%state(i_particle, col_bathy) > srcPar%P_D(i_particle)) then
			
				call self%NonlinearSolverVectorCoupling%NewtonSolver(momentumBalanceFloater, x0_vector, root_vector, ok_vector)

				if (.not. all(ok_vector)) then
					call self%NonlinearSolverVectorCoupling%SecantSolver(momentumBalanceFloater, x0_vector, x1_vector, root_vector, ok_vector)
				end if

				if (all(ok_vector)) then
					vel_floater_new(i_particle,1) = root_vector(1)
					vel_floater_new(i_particle,2) = root_vector(2)
					vel_floater_new(i_particle,3) = 0.0_prec
				elseif (.not.all(ok_vector)) then
					print*, "NOT all(ok_vector)",i_particle
					vel_floater_new(i_particle,1) = 1.05 * x0_vector(1)
					vel_floater_new(i_particle,2) = 1.05 * x0_vector(2)
					vel_floater_new(i_particle,3) = 0.0_prec
				end if
				
			else
				! print*,"BEACHED --->",i_particle,"bath:",-sv%state(i_particle, col_bathy),"P_D:",srcPar%P_D(i_particle)
				vel_floater_new(i_particle,1) = 0.0_prec
				vel_floater_new(i_particle,2) = 0.0_prec
				vel_floater_new(i_particle,3) = 0.0_prec
			end if			
	
			! write(*,'(F12.4, I10, L10, F12.8, F12.8, F12.8)'), time, i_particle ,ok_vector(1),  sv%state(i_particle,4), WindVel(i_particle, 1), root_vector(1)
			! write(*,'(F12.4, I10, L10, F12.8, F12.8, F12.8)'), time, i_particle ,ok_vector(2),  sv%state(i_particle,5), WindVel(i_particle, 2), root_vector(2)
				
		end do
			
		sv%state(:, col_ufloater) = vel_floater_new(:,1)
		sv%state(:, col_vfloater) = vel_floater_new(:,2)	
		sv%state(:, col_wfloater) = vel_floater_new(:,3)	
		
		FloaterVelocity(:,1) = Utils%m2geo(sv%state(:, col_ufloater), sv%state(:,2), .false.)
		FloaterVelocity(:,2) = Utils%m2geo(sv%state(:, col_vfloater), sv%state(:,2), .true.)
		FloaterVelocity(:,3) = 0.0_prec


	contains

		!---------------------------------------------------------------------------
		!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
		!> @brief
		!> Computes the  momentum balance for a floter tracer
		!> @param[in] self, sv, bdata, time
		!--------------------------------------------------------------------------	  
		function momentumBalanceFloater(x) result(f)
			real(prec), intent(in) :: x(:)
			real(prec), dimension(size(x)) :: f
			real(prec), dimension(size(x)) :: rel_WindVel_vec, rel_WindVel_vec_abs
			real(prec), dimension(size(x)) :: rel_CurrentVel_vec, rel_CurrentVel_vec_abs
			real(prec), dimension(size(x)) :: rel_PorousVel_vec, rel_PorousVel_vec_abs
			real(prec), dimension(size(x)) :: KWind_vec, FWind_vec
			real(prec), dimension(size(x)) :: KCurrent_vec, FCurrent_vec 
			real(prec), dimension(size(x)) :: KPorous_vec, FPorous_vec			
			
			if (hasWind) then
				rel_WindVel_vec(1) = WindVel(i_particle,1) - x(1)
				rel_WindVel_vec(2) = WindVel(i_particle,2) - x(2)
				rel_WindVel_vec_abs(1) = abs(rel_WindVel_vec(1))
				rel_WindVel_vec_abs(2) = abs(rel_WindVel_vec(2))
				KWind_vec = self%DragCoefficient_K_wind_vector(i_particle, rel_WindVel_vec_abs, srcPar, KViscosity_air_vector_atSurfc, Density_air_vector_atSurfc)			
			else
				rel_WindVel_vec = 0.0_prec
				rel_WindVel_vec_abs = 0.0_prec
				KWind_vec = 0.0_prec
			end if 
						
			rel_CurrentVel_vec(1) = u0 + StokesVel(i_particle,1)  - x(1)
			rel_CurrentVel_vec(2) = v0 + StokesVel(i_particle,2) - x(2)
			rel_CurrentVel_vec_abs(1) = abs(rel_CurrentVel_vec(1))
			rel_CurrentVel_vec_abs(2) = abs(rel_CurrentVel_vec(2))
			KCurrent_vec = self%DragCoefficient_K_current_vector(i_particle, rel_CurrentVel_vec_abs, srcPar, KViscosity_seawater_vector_atSurfc, Density_seawater_vector_atSurfc)


			rel_PorousVel_vec(1) = uP_D - x(1)
			rel_PorousVel_vec(2) = vP_D - x(2)
			rel_PorousVel_vec_abs(1) = abs(rel_PorousVel_vec(1))
			rel_PorousVel_vec_abs(2) = abs(rel_PorousVel_vec(2))
			KPorous_vec = self%DragCoefficient_K_porous_vector(i_particle, rel_PorousVel_vec_abs, srcPar, KViscosity_seawater_vector_atDepth, Density_seawater_vector_atDepth)


			FWind_vec 		= KWind_vec * rel_WindVel_vec
			FCurrent_vec	= KCurrent_vec * rel_CurrentVel_vec
			FPorous_vec		= KPorous_vec * rel_PorousVel_vec
			! print*, "FPorous_vec:", FPorous_vec
			! print*, "FCurrent_vec:", FCurrent_vec			

			f(1) = x(1)  - sv%state(i_particle, col_ufloater) - (dt/srcPar%Mass(i_particle)) * (FWind_vec(1) + FCurrent_vec(1) + FPorous_vec(1) )
			f(2) = x(2)  - sv%state(i_particle, col_vfloater) - (dt/srcPar%Mass(i_particle)) * (FWind_vec(2) + FCurrent_vec(2) + FPorous_vec(2) )

		end function momentumBalanceFloater

	end function FloaterVelocity

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_K_wind_vector of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------
	function DragCoefficient_K_wind_vector(self, idx, vector, srcPar, KViscosity_air_vector_atSurfc, Density_air_vector_atSurfc) result(K_w_vector)
	class(kernelFloater_class), intent(in)		:: self
	integer, intent(in) 						:: idx
	real(prec), dimension(:), intent(in) 		:: vector
	type(source_parameters_type), intent(in)	:: srcPar
	real(prec), dimension(:), intent(in)		:: KViscosity_air_vector_atSurfc, Density_air_vector_atSurfc
	real(prec), dimension(size(vector)) 		:: K_w_vector
	real(prec)									:: CD_Constant_scalar, CD_Porous_scalar
	real(prec)									:: CD_Sphere_scalar,CD_Rectangular_scalar
	real(prec)									:: characteristicVelocity, Reynolds, D_eff_X


		characteristicVelocity = sqrt(vector(1) **2.0_prec +  vector(2) **2.0_prec)

		!for sphere: D_eff = S_R
		D_eff_X = srcPar%S_R(idx)
		Reynolds = self%Reynolds_X(characteristicVelocity, KViscosity_air_vector_atSurfc(idx), D_eff_X)
		CD_Sphere_scalar = self%DragCoefficient_CD_Sphere_scalar(idx, Reynolds, srcPar)

		!for Rectangular: D_eff = R_Lz
		D_eff_X = srcPar%R_Lz(idx)
		Reynolds = self%Reynolds_X(characteristicVelocity, KViscosity_air_vector_atSurfc(idx), D_eff_X)
		CD_Rectangular_scalar = self%DragCoefficient_CD_Rectangular_scalar(idx, Reynolds, srcPar)


		K_w_vector(1) = 0.5_prec * Density_air_vector_atSurfc(idx) * (srcPar%S_N(idx) * CD_Sphere_scalar * srcPar%S_rd(idx) * 3.14159_prec * (srcPar%S_R(idx)*srcPar%S_R(idx)) + &
																    srcPar%R_N(idx) * CD_Rectangular_scalar * srcPar%R_rd(idx) * (srcPar%R_Ly(idx)*srcPar%R_Lz(idx))) * vector(1)
		K_w_vector(2) = 0.5_prec * Density_air_vector_atSurfc(idx) * (srcPar%S_N(idx) * CD_Sphere_scalar * srcPar%S_rd(idx) * 3.14159_prec * (srcPar%S_R(idx)*srcPar%S_R(idx)) + &
																    srcPar%R_N(idx) * CD_Rectangular_scalar * srcPar%R_rd(idx) * (srcPar%R_Lx(idx)*srcPar%R_Lz(idx))) * vector(2)
	end function DragCoefficient_K_wind_vector

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_K_current_vector of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------
	function DragCoefficient_K_current_vector(self, idx, vector, srcPar, KViscosity_seawater_vector_atSurfc, Density_seawater_vector_atSurfc) result(K_c_vector)
	integer, intent(in)							:: idx
	class(kernelFloater_class), intent(in)		:: self
	real(prec), dimension(:), intent(in) 		:: vector
	type(source_parameters_type), intent(in)	:: srcPar
	real(prec), dimension(:), intent(in)		:: KViscosity_seawater_vector_atSurfc, Density_seawater_vector_atSurfc
	real(prec), dimension(size(vector))			:: K_c_vector
	real(prec)									:: CD_Constant_scalar, CD_Porous_scalar
	real(prec)									:: CD_Sphere_scalar,CD_Rectangular_scalar
	real(prec)									:: characteristicVelocity, Reynolds, D_eff_X	  

		characteristicVelocity = sqrt(vector(1) **2.0_prec +  vector(2) **2.0_prec)

		!for Shpere: D_eff = S_R
		D_eff_X	= srcPar%S_R(idx)
		Reynolds = self%Reynolds_X(characteristicVelocity, KViscosity_seawater_vector_atSurfc(idx), D_eff_X)
		CD_Sphere_scalar = self%DragCoefficient_CD_Sphere_scalar(idx, Reynolds, srcPar)
	
		! print*, "-----------------------------"
		! print*,"ReynoldsW:", Reynolds
		! print*, "-----------------------------"	
		! print*, "-----------------------------"
		! print*,"CD_Sphere_scalarW:", CD_Sphere_scalar
		! print*, "-----------------------------"	
	
		!for Rectangular: D_eff = R_Lz
		D_eff_X = srcPar%R_Lz(idx)
		Reynolds = self%Reynolds_X(characteristicVelocity, KViscosity_seawater_vector_atSurfc(idx), D_eff_X)
		CD_Rectangular_scalar = self%DragCoefficient_CD_Rectangular_scalar(idx, Reynolds, srcPar)	

		! print*, "-----------------------------"
		! print*,"ReynoldsC:", Reynolds
		! print*, "-----------------------------"	
		! print*, "-----------------------------"
		! print*,"CD_Rectangular_scalarC:", CD_Rectangular_scalar
		! print*, "-----------------------------"	
	
		K_c_vector(1) = 0.5_prec * Density_seawater_vector_atSurfc(idx) * (srcPar%S_N(idx) * CD_Sphere_scalar * (1.0_prec - srcPar%S_rd(idx)) * 3.14159_prec * (srcPar%S_R(idx)*srcPar%S_R(idx)) + &
																    srcPar%R_N(idx) * CD_Rectangular_scalar * (1.0_prec - srcPar%R_rd(idx)) * (srcPar%R_Ly(idx)*srcPar%R_Lz(idx))) * vector(1)
		K_c_vector(2) = 0.5_prec * Density_seawater_vector_atSurfc(idx) * (srcPar%S_N(idx) * CD_Sphere_scalar * (1.0_prec - srcPar%S_rd(idx)) * 3.14159_prec * (srcPar%S_R(idx)*srcPar%S_R(idx)) + &
																    srcPar%R_N(idx) * CD_Rectangular_scalar * (1.0_prec - srcPar%R_rd(idx)) * (srcPar%R_Lx(idx)*srcPar%R_Lz(idx))) * vector(2)

	
	end function DragCoefficient_K_current_vector


    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_K_porous_vector of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------
	function DragCoefficient_K_porous_vector(self, idx, vector, srcPar, KViscosity_seawater_vector_atDepth, Density_seawater_vector_atDepth) result(K_p_vector)
	integer, intent(in)							:: idx
	class(kernelFloater_class), intent(in)		:: self
	real(prec), dimension(:), intent(in) 		:: vector
	type(source_parameters_type), intent(in)	:: srcPar
	real(prec), dimension(:), intent(in)		:: KViscosity_seawater_vector_atDepth, Density_seawater_vector_atDepth
	real(prec), dimension(size(vector))			:: K_p_vector
	real(prec)									:: CD_Constant_scalar, CD_Porous_scalar
	real(prec)									:: CD_Sphere_scalar,CD_Rectangular_scalar
	real(prec)									:: characteristicVelocity ,Reynolds, D_eff_X	  

		characteristicVelocity = sqrt(vector(1) **2.0_prec +  vector(2) **2.0_prec)

		
		!for Rectangular: D_eff = P_Lz
		D_eff_X = 2.0_prec * srcPar%P_Lz(idx)
		Reynolds = self%Reynolds_X(characteristicVelocity, KViscosity_seawater_vector_atDepth(idx), D_eff_X)
		CD_Porous_scalar = self%DragCoefficient_CD_Porous_scalar(idx, Reynolds, srcPar)

		K_p_vector(1) = 0.5_prec * Density_seawater_vector_atDepth(idx) * (srcPar%P_N(idx) * CD_Porous_scalar * (1.0_prec - srcPar%P_rd(idx)) * (srcPar%P_Ly(idx)*srcPar%P_Lz(idx))) * vector(1)
		K_p_vector(2) = 0.5_prec * Density_seawater_vector_atDepth(idx) * (srcPar%P_N(idx) * CD_Porous_scalar * (1.0_prec - srcPar%P_rd(idx)) * (srcPar%P_Lx(idx)*srcPar%P_Lz(idx))) * vector(2)
	
	end function DragCoefficient_K_porous_vector
	
    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_CD_Sphere_scalar of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------	
	function DragCoefficient_CD_Sphere_scalar(self, idx, value_in, srcPar) result(CD_S_value)
	class(kernelFloater_class), intent(in)		:: self
	integer, intent(in)							:: idx
	real(prec), intent(in)						:: value_in
	type(source_parameters_type), intent(in)	:: srcPar
	real(prec) 									:: CD_S_value, value

		if (srcPar%S_CD(idx) == MV) then
			value = max(min(value_in, 1.0E6_prec), 1.0E-8_prec)
			CD_S_value = (24.0_prec/ value) + 2.6_prec * (value/5.0_prec) / (1.0_prec + (value/5.0_prec)**1.52_prec) + &
							0.411_prec * (value/(2.63E5_prec)) ** (-7.94_prec) / (1.0_prec + (value/(2.63E5_prec)) ** -8.0_prec) + &
							0.25_prec * (value/1.0E6_prec) / (1.0_prec + (value/1.0E6_prec)) 
		else
			CD_S_value = srcPar%S_CD(idx)
		end if	
		
	end function DragCoefficient_CD_Sphere_scalar

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_CD_Rectangular_scalar of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------
	
	function DragCoefficient_CD_Rectangular_scalar(self, idx, value_in, srcPar) result(CD_R_vlue)
	class(kernelFloater_class), intent(in)		:: self
	integer, intent(in)							:: idx
	real(prec), intent(in)						:: value_in
	type(source_parameters_type), intent(in)	:: srcPar
	real(prec)									:: CD_R_vlue, value

	
		if (srcPar%R_CD(idx) == MV) then
			value = max(min(value_in, 1.0E6_prec), 1.0E-8_prec)
			CD_R_vlue = 1.2_prec
		else
			CD_R_vlue = srcPar%R_CD(idx)
		end if
	end function DragCoefficient_CD_Rectangular_scalar

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_CD_Porous_scalar of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------
	
	function DragCoefficient_CD_Porous_scalar(self, idx, value_in, srcPar) result(CD_P_scalar)
	class(kernelFloater_class), intent(in)		:: self
	integer, intent(in)							:: idx
	real(prec), intent(in) 						:: value_in
	type(source_parameters_type), intent(in)	:: srcPar
	real(prec)									:: CD_P_scalar, value

		if (srcPar%P_CD(idx) == MV) then
			value = max(min(value_in, 1.0E6_prec), 1.0E-8_prec)
			CD_P_scalar = 1.5_prec	
		else
			CD_P_scalar = srcPar%P_CD(idx)
		end if

	end function DragCoefficient_CD_Porous_scalar

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DragCoefficient_CD_Constant_scalar of a tracers in seawater
    !> @param[in] self, sv, bdata, time
	!> vector is the magnitude of the reletive velocity :sqrt(vx**2 + vy**2)
    !---------------------------------------------------------------------------
	
	function DragCoefficient_CD_Constant_scalar(self, idx, value_in) result(CD_C_scalar)
	class(kernelFloater_class), intent(in)		:: self
	integer, intent(in)							:: idx
	real(prec) , intent(in)						:: value_in
	real(prec)									:: CD_C_scalar, value

		CD_C_scalar = 1.5_prec
	end function DragCoefficient_CD_Constant_scalar

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DensityKViscosity_seawater of a tracers in seawater
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    subroutine DensityKViscosity_seawater(self, sv, KViscosity_seawater_vector, Density_seawater_vector)
    class(kernelFloater_class), intent(in) 	:: self
    type(stateVector_class), intent(in) 		:: sv
	real(prec), allocatable, intent(inout) 	:: KViscosity_seawater_vector(:)
	real(prec), allocatable, intent(inout) 	:: Density_seawater_vector(:)
	real(prec), allocatable :: Temperatur_list(:), Salinity_list(:), depth_list(:)
	real(prec), allocatable :: fDensity(:), kVisco(:), kViscoRelation(:)
	integer :: col_temp, col_sal

	allocate(Temperatur_list(size(sv%state,1)))
	allocate(Salinity_list(size(sv%state,1)))	
	allocate(depth_list(size(sv%state,1)))	
	allocate(fDensity(size(sv%state,1)))
	allocate(kVisco(size(sv%state,1)))	
	allocate(kViscoRelation(size(sv%state,1)))	
    
		col_temp = Utils%find_str(sv%varname, Globals%Var%temp, .false.)
		col_sal = Utils%find_str(sv%varname, Globals%Var%sal, .false.)

	!	Calculate the kinematic Visco and density of seawater
		if ((col_temp /= MV_INT) .and. (col_sal /= MV_INT)) then
		! If there are salt and temperature data along the time:
		!	Calculate the kinematic Visco and density of seawater according to 
		!	the given temperature and salinity of seawater during the time. 
					
			Temperatur_list = sv%state(:,col_temp)
			Salinity_list	= sv%state(:,col_sal)
			depth_list		= sv%state(:,3)
			fDensity = self%VerticalMotion%seaWaterDensity(Salinity_list, Temperatur_list, depth_list)
			kVisco = self%VerticalMotion%absoluteSeaWaterViscosity(Salinity_list, Temperatur_list) / fDensity

			kViscoRelation = abs(1.0_prec - (kVisco/Globals%Constants%MeanKvisco))
			where(kViscoRelation >= 0.9_prec)
				kVisco = Globals%Constants%MeanKvisco
				fDensity = Globals%Constants%MeanDensity
			endwhere
		else
			! If there is no salt and temperature:
			! Consider the kinematic Visco and density of sea water according to 
			! the given constant values in Globals%Constants

			kVisco = Globals%Constants%MeanKvisco
			fDensity = Globals%Constants%MeanDensity
		end if

		call ensure_real_alloc(KViscosity_seawater_vector, size(sv%state,1))
		call ensure_real_alloc(Density_seawater_vector,    size(sv%state,1))

		KViscosity_seawater_vector = kVisco
		Density_seawater_vector = fDensity
	
    end subroutine DensityKViscosity_seawater
	
    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  DensityKViscosity_air of a tracers in seawater
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    subroutine DensityKViscosity_air(self, sv, KViscosity_air_vector, Density_air_vector)
    class(kernelFloater_class), intent(in) 	:: self
    type(stateVector_class), intent(in) 		:: sv
	real(prec), allocatable, intent(inout) 	:: KViscosity_air_vector(:)
	real(prec), allocatable, intent(inout) 	:: Density_air_vector(:)

		call ensure_real_alloc(KViscosity_air_vector, size(sv%state,1))
		call ensure_real_alloc(Density_air_vector,    size(sv%state,1))

		KViscosity_air_vector = 1.5E-5_prec
		Density_air_vector = 1.12_prec
	
    end subroutine DensityKViscosity_air

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Computes the  Reynolds_X of a tracers in seawater
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
    function Reynolds_X(self, characteristicVelocity, kvisco, characteristicLength) result(Reynolds)
    class(kernelFloater_class), intent(in) :: self
    real(prec), intent(in) :: characteristicVelocity, kvisco, characteristicLength
    real(prec)	:: Reynolds
    integer :: id
	
    !To avoid reynolds with 0 value (which will give infinit value in DragCoefficient 
    Reynolds = max(abs(characteristicVelocity), 1.0E-8_prec)*characteristicLength/kvisco
    !Reynolds = abs(characteristicVelocity)*characteristicLength/kvisco

    end function Reynolds_X

    !---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Find the source parameters and  put them as the properties of the corresponding tracer 
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------

	subroutine findSourceParameters(self, sv, srcPar)
		implicit none
		class(kernelFloater_class), intent(in) 	:: self
		type(stateVector_class), intent(in) 	:: sv
		type(source_parameters_type), intent(inout) :: srcPar
		integer 		:: i, n
		type(string) 	:: outext

		n = size(sv%state,1)
		!n = size(sv%trc)

		! if (size(sv%trc) < n) then
		   ! write(*,*) 'ERROR: size(sv%trc) < size(sv%state,1)'
		   ! stop
		! end if
		if (.not. allocated(sv%trc)) then
			outext = '[kernelFloater::findSourceParameters] State vector tracer pointer array is not allocated'
			call Log%put(outext)
			stop
		end if
		if (size(sv%trc) < n) then
			outext = '[kernelFloater::findSourceParameters] size(sv%trc) < size(sv%state,1), stopping'
			call Log%put(outext)
			stop
		end if
		
		call ensure_real_alloc(srcPar%Mass, n)

		call ensure_int_alloc (srcPar%S_N,  n)
		call ensure_real_alloc(srcPar%S_rd, n)
		call ensure_real_alloc(srcPar%S_CD, n)
		call ensure_real_alloc(srcPar%S_R,  n)

		call ensure_int_alloc (srcPar%R_N,  n)
		call ensure_real_alloc(srcPar%R_rd, n)
		call ensure_real_alloc(srcPar%R_CD, n)
		call ensure_real_alloc(srcPar%R_Lx, n)
		call ensure_real_alloc(srcPar%R_Ly, n)
		call ensure_real_alloc(srcPar%R_Lz, n)

		call ensure_int_alloc (srcPar%P_N,  n)
		call ensure_real_alloc(srcPar%P_rd, n)
		call ensure_real_alloc(srcPar%P_CD, n)
		call ensure_real_alloc(srcPar%P_Lx, n)
		call ensure_real_alloc(srcPar%P_Ly, n)
		call ensure_real_alloc(srcPar%P_Lz, n)
		call ensure_real_alloc(srcPar%P_D,  n)
		
		srcPar%Mass = 0.0_prec
		srcPar%S_N  = 0
		srcPar%S_rd = 0.0_prec
		srcPar%S_CD = 0.0_prec
		srcPar%S_R  = 0.0_prec

		srcPar%R_N  = 0
		srcPar%R_rd = 0.0_prec
		srcPar%R_CD = 0.0_prec
		srcPar%R_Lx = 0.0_prec
		srcPar%R_Ly = 0.0_prec
		srcPar%R_Lz = 0.0_prec

		srcPar%P_N  = 0
		srcPar%P_rd = 0.0_prec
		srcPar%P_CD = 0.0_prec
		srcPar%P_Lx = 0.0_prec
		srcPar%P_Ly = 0.0_prec
		srcPar%P_Lz = 0.0_prec
		srcPar%P_D  = 0.0_prec
		
		do i = 1, n
			!if (.not. associated(sv%trc(i)%ptr)) cycle
			if (.not. associated(sv%trc(i)%ptr)) then
				outext = '[kernelFloater::findSourceParameters] Tracer pointer not associated for floater state vector entry, stopping'
				call Log%put(outext)
				stop
			end if

			select type(aTracer => sv%trc(i)%ptr)
			type is (floater_class)
				srcPar%Mass(i) = aTracer%mpar%Mass
				srcPar%S_N(i)  = aTracer%mpar%Sphere_N
				srcPar%S_rd(i) = aTracer%mpar%Sphere_ratio_dry
				srcPar%S_CD(i) = aTracer%mpar%Sphere_CD
				srcPar%S_R(i)  = aTracer%mpar%Sphere_Radius

				srcPar%R_N(i)  = aTracer%mpar%Rectangular_N
				srcPar%R_rd(i) = aTracer%mpar%Rectangular_ratio_dry
				srcPar%R_CD(i) = aTracer%mpar%Rectangular_CD
				srcPar%R_Lx(i) = aTracer%mpar%Rectangular_Lx
				srcPar%R_Ly(i) = aTracer%mpar%Rectangular_Ly
				srcPar%R_Lz(i) = aTracer%mpar%Rectangular_Lz

				srcPar%P_N(i)  = aTracer%mpar%Porous_N
				srcPar%P_rd(i) = aTracer%mpar%Porous_ratio_dry
				srcPar%P_CD(i) = aTracer%mpar%Porous_CD
				srcPar%P_Lx(i) = aTracer%mpar%Porous_Lx
				srcPar%P_Ly(i) = aTracer%mpar%Porous_Ly
				srcPar%P_Lz(i) = aTracer%mpar%Porous_Lz
				srcPar%P_D(i)  = aTracer%mpar%Porous_Depth
			class default
				outext = '[kernelFloater::findSourceParameters] Non-floater tracer found in floater state vector, stopping'
				call Log%put(outext)
				stop
			end select
		end do
		
		if (any(srcPar%Mass <= 0.0_prec)) then
			outext = '[kernelFloater::findSourceParameters] Found floater with non-positive mass, stopping'
			call Log%put(outext)
			stop
		end if
		
	end subroutine findSourceParameters

	!---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Allocating the vectors
    !> @param[in] self, sv, bdata, time
    !---------------------------------------------------------------------------
	subroutine ensure_real_alloc(arr, n)
		real(prec), allocatable, intent(inout) :: arr(:)
		integer, intent(in) :: n

		if (.not. allocated(arr)) then
			allocate(arr(n))
		else if (size(arr) /= n) then
			deallocate(arr)
			allocate(arr(n))
		end if
	end subroutine ensure_real_alloc

	subroutine ensure_int_alloc(arr, n)
		integer, allocatable, intent(inout) :: arr(:)
		integer, intent(in) :: n

		if (.not. allocated(arr)) then
			allocate(arr(n))
		else if (size(arr) /= n) then
			deallocate(arr)
			allocate(arr(n))
		end if
	end subroutine ensure_int_alloc

	!---------------------------------------------------------------------------
	!> @author Mohsen Shabani - CoLab+Atlantic- 2026.05.01 | Email:shabani.mohsen@outlook.com	
    !> @brief
    !> Initializer method adpated from for kernel class. Sets the type of
    !> kernel and the interpolator to evaluate it.
    !---------------------------------------------------------------------------
	subroutine initKernelFloater(self)
	class(kernelFloater_class), intent(inout) :: self
	type(string) :: interpName

	interpName = 'linear'
	call self%Interpolator%initialize(1, interpName)
	
	call self%SharedProcesses%initialize()
	call self%KernelUtils%initialize()
	call self%VerticalMotion%initialize()
	
	end subroutine initKernelFloater

	end module kernelFloater_mod