
module tracer3D 

    use tracer_precision
    use tracer_interp 
    use bspline_module, only : bspline_3d    
    use nml 

    implicit none 

    type tracer_par_trans_class
        integer :: nt 

        real(prec), allocatable :: time(:)
        real(prec), allocatable :: H_min_dep(:)
        real(prec), allocatable :: dt_dep(:)
        integer,    allocatable :: n_max_dep(:) 
        real(prec), allocatable :: dt_write(:)
        
    end type 

    type tracer_par_class 
        integer :: n, n_active, n_max_dep, id_max 
        logical :: is_sigma                     ! Is the defined z-axis in sigma coords
        real(prec_time) :: dt, dt_dep, dt_write 
        real(prec) :: thk_min                   ! Minimum thickness of tracer (m)
        real(prec) :: H_min                     ! Minimum ice thickness to track (m)
        real(prec) :: depth_max                 ! Maximum depth of tracer (fraction)
        real(prec) :: U_max                     ! Maximum horizontal velocity of tracer to track (m/a)
        real(prec) :: U_max_dep                 ! Maximum horizontal velocity allowed for tracer deposition (m/a)
        real(prec) :: H_min_dep                 ! Minimum ice thickness for tracer deposition (m)
        real(prec) :: alpha                     ! Slope of probability function
        character(len=56) :: weight             ! Weighting function for generating prob. distribution
        logical    :: noise                     ! Add noise to gridded deposition location
        real(prec) :: dens_z_lim                ! Distance from surface to count density
        integer    :: dens_max                  ! Max allowed density of particles at surface
        character(len=56) :: interp_method  

        ! Transient parameters 
        character(len=512) :: par_trans_file 
        logical            :: use_par_trans
        type(tracer_par_trans_class) :: tpar 

    end type 

    type tracer_state_class
        real(prec_time) :: time, time_old
        real(prec_time) :: time_dep, time_write 
        real(prec_time) :: dt  
        integer, allocatable :: active(:), id(:)
        real(prec), allocatable :: x(:), y(:), z(:), sigma(:)
        real(prec), allocatable :: ux(:), uy(:), uz(:)
        real(prec), allocatable :: ax(:), ay(:), az(:)
        real(prec), allocatable :: dpth(:), z_srf(:)
        real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
        real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
        real(prec), allocatable :: H(:)

    end type 

    type tracer_stats_class
        ! All stats variable at precision of writing (prec_wrt), since
        ! this should not need high precision output 

        real(prec_wrt), allocatable :: x(:), y(:)
        real(prec_wrt), allocatable :: depth_norm(:)
        real(prec_wrt), allocatable :: age_iso(:) 

        real(prec_wrt), allocatable :: depth_iso(:,:,:)
        real(prec_wrt), allocatable :: depth_iso_err(:,:,:)
        real(prec_wrt), allocatable :: dep_z_iso(:,:,:)
        integer,    allocatable :: density_iso(:,:,:)
        
        real(prec_wrt), allocatable :: ice_age(:,:,:)
        real(prec_wrt), allocatable :: ice_age_err(:,:,:)
        integer,    allocatable :: density(:,:,:)
           
    end type

    type tracer_dep_class 
        ! Standard deposition information (time and place)
        real(prec), allocatable :: time(:) 
        real(prec), allocatable :: H(:) 
        real(prec), allocatable :: x(:), y(:), z(:)
        real(prec), allocatable :: lon(:), lat(:) 

        ! Additional tracer deposition information (climate, isotopes, etc)
        real(prec), allocatable :: t2m_ann(:), t2m_sum(:)      
        real(prec), allocatable :: pr_ann(:), pr_sum(:)     
        real(prec), allocatable :: t2m_prann(:) ! Precip-weighted temp
        real(prec), allocatable :: d18O_ann(:)
!         real(prec), allocatable :: dD(:)

    end type 

    type tracer_class 
        type(tracer_par_class)   :: par 
        type(tracer_state_class) :: now 
        type(tracer_dep_class)   :: dep
        type(tracer_stats_class) :: stats 

    end type 

    type(lin3_interp_par_type) :: par_lin 
    type(bspline_3d)           :: bspline3d_ux, bspline3d_uy, bspline3d_uz 

    private 

    ! For other tracer modules
    public :: tracer_par_class
    public :: tracer_state_class
    public :: tracer_dep_class
    public :: tracer_stats_class

    public :: tracer_reshape1D_vec
    public :: tracer_reshape2D_field 
    public :: tracer_reshape3D_field

    ! General public 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 
    public :: tracer_end 

contains 

    subroutine tracer_init(trc,filename,time,x,y,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:), y(:)
        logical,    intent(IN) :: is_sigma  
        real(prec_time), intent(IN) :: time 

        ! Local variables 
        integer :: i 

        ! Load the parameters
        call tracer_par_load(trc%par,filename,is_sigma)

        ! Update the transient parameters
        if (trc%par%use_par_trans) then 
            call tracer_par_update(trc%par,trc%par%tpar,time)
        end if 

        ! Allocate the state variables 
        call tracer_allocate(trc%now,trc%dep,n=trc%par%n)
        call tracer_allocate_stats(trc%stats,x,y)

        ! ===== Initialize stats depth axes ===============

        trc%stats%age_iso = [11.7,29.0,57.0,115.0,130.0]

        do i = 1, size(trc%stats%depth_norm)
            trc%stats%depth_norm(i) = 0.04*real(i)
        end do 

        ! =================================================

        ! Initialize state 
        trc%now%active    = 0 

        trc%now%id        = mv 
        trc%now%x         = mv 
        trc%now%y         = mv 
        trc%now%z         = mv 
        trc%now%sigma     = mv 
        trc%now%z_srf     = mv 
        trc%now%dpth      = mv 
        trc%now%ux        = mv 
        trc%now%uy        = mv 
        trc%now%uz        = mv 
        trc%now%ax        = mv 
        trc%now%ay        = mv 
        trc%now%az        = mv 
        trc%now%thk       = mv 
        trc%now%T         = mv 
        trc%now%H         = mv 

        trc%dep%time      = mv 
        trc%dep%H         = mv 
        trc%dep%x         = mv 
        trc%dep%y         = mv 
        trc%dep%z         = mv 

        trc%par%id_max    = 0 

        ! Initialize the time (to one older than now)
        trc%now%time   = time - 1000.0_dp
        trc%now%time_dep   = time - 1000.0 
        trc%now%time_write = time - 1000.0 

        ! Initialize random number generator 
        call random_seed() 

        return 

    end subroutine tracer_init

    subroutine tracer_update(trc,time,x,y,z,z_srf,H,ux,uy,uz,              &
                             lon,lat,t2m_ann,t2m_sum,pr_ann,pr_sum,d18O_ann, &
                             dep_now,stats_now,order,sigma_srf)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
        real(prec), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
        real(prec), intent(IN) :: lon(:,:), lat(:,:), t2m_ann(:,:), t2m_sum(:,:), pr_ann(:,:), pr_sum(:,:), d18O_ann(:,:) 
        logical, intent(IN) :: dep_now, stats_now  
        character(len=*), intent(IN), optional :: order 
        real(prec), intent(IN), optional :: sigma_srf     ! Value at surface by default (1 or 0?)

        ! Local variables  
        character(len=3) :: idx_order 
        integer    :: i, j, k, nx, ny, nz
        logical    :: rev_z 
        real(prec), allocatable :: x1(:), y1(:), z1(:)
        real(prec), allocatable :: zc(:)   ! Actual cartesian z-axis after applying sigma*H 
        real(prec), allocatable :: z_srf1(:,:), H1(:,:)
        real(prec), allocatable :: ux1(:,:,:), uy1(:,:,:), uz1(:,:,:)
        real(prec), allocatable :: usig1(:,:,:)
        real(prec), allocatable :: lon1(:,:), lat1(:,:), t2m_ann1(:,:), t2m_sum1(:,:), pr_ann1(:,:), pr_sum1(:,:), d18O_ann1(:,:)
        real(prec) :: ux0, uy0, uz0 
        real(prec) :: dt 

        ! Update the transient parameters
        if (trc%par%use_par_trans) then 
            call tracer_par_update(trc%par,trc%par%tpar,time)
        end if 

        ! Update current time and time step
        trc%now%time_old = trc%now%time 
        trc%now%time     = time 
        trc%now%dt       = real(dble(trc%now%time) - dble(trc%now%time_old),prec_time)

        ! Update record of last deposition time if dep_now
        if (dep_now) trc%now%time_dep = trc%now%time 

        ! Determine order of indices (default ijk)
        idx_order = "ijk"
        if (present(order)) idx_order = trim(order)

        ! Allocate helper z-axis variable 
        if (allocated(zc)) deallocate(zc)
        allocate(zc(size(z)))

        zc = z 

        if (trc%par%is_sigma) then 
            ! Ensure z-axis is properly bounded
            where (abs(zc) .lt. 1e-5) zc = 0.0 

            if (minval(zc) .lt. 0.0 .or. maxval(zc) .gt. 1.0) then 
                write(0,*) "tracer:: error: sigma axis not bounded between zero and one."
                write(0,*) "z = ", zc 
                stop 
            end if

        end if 

        ! Note: GRISLI (nx,ny,nz): sigma goes from 1 to 0, so sigma(1)=1 [surface], sigma(nz)=0 [base]
        !       SICO (nz,ny,nx): sigma(1) = 0, sigma(nz) = 1
        ! reshape routines ensure ascending z-axis (nx,ny,nz) with sigma(nz)=1 [surface]
        
        ! Correct the sigma values if necessary,
        ! so that sigma==0 [base]; sigma==1 [surface]
        if (trc%par%is_sigma .and. present(sigma_srf)) then 
            if (sigma_srf .eq. 0.0) then 
                ! Adjust sigma values 
                zc = 1.0 - zc 
            end if 
        end if 

        ! Also determine whether z-axis is initially ascending or descending 
        rev_z = (zc(1) .gt. zc(size(zc)))

        call tracer_reshape1D_vec(x, x1,rev=.FALSE.)
        call tracer_reshape1D_vec(y, y1,rev=.FALSE.)
        call tracer_reshape1D_vec(real(zc,kind=prec),z1,rev=rev_z)
        call tracer_reshape2D_field(idx_order,z_srf,z_srf1)
        call tracer_reshape2D_field(idx_order,H,H1)
        call tracer_reshape3D_field(idx_order,ux,ux1,rev_z=rev_z)
        call tracer_reshape3D_field(idx_order,uy,uy1,rev_z=rev_z)
        call tracer_reshape3D_field(idx_order,uz,uz1,rev_z=rev_z)
        
        if (dep_now) then 
            ! Also reshape deposition fields 
            call tracer_reshape2D_field(idx_order,lon,lon1)
            call tracer_reshape2D_field(idx_order,lat,lat1)
            call tracer_reshape2D_field(idx_order,t2m_ann,t2m_ann1)
            call tracer_reshape2D_field(idx_order,t2m_sum,t2m_sum1)
            call tracer_reshape2D_field(idx_order,pr_ann,pr_ann1)
            call tracer_reshape2D_field(idx_order,pr_sum,pr_sum1)
            call tracer_reshape2D_field(idx_order,d18O_ann,d18O_ann1)
            
        end if 

        ! Get axis sizes (if ny==2, this is a 2D profile domain) 
        nx = size(x1,1)
        ny = size(y1,1)
        nz = size(z1,1)

        if (trim(trc%par%interp_method) .eq. "spline") then

            ! Allocate z-velocity field in sigma coordinates 
            if (allocated(usig1)) deallocate(usig1)
            allocate(usig1(nx,ny,nz))

            usig1 = 0.0
            do k = 1, nz 
                where (H1 .gt. 0.0) usig1(:,:,k) = uz1(:,:,k) / H1
            end do 

            call interp_bspline3D_weights(bspline3d_ux,x1,y1,z1,ux1)
            call interp_bspline3D_weights(bspline3d_uy,x1,y1,z1,uy1)
            call interp_bspline3D_weights(bspline3d_uz,x1,y1,z1,usig1)
            
            write(*,*) "spline weights calculated."
        end if 

        ! Interpolate to the get the right elevation and other deposition quantities
        do i = 1, trc%par%n 

            if (trc%now%active(i) .eq. 2) then 

                ! Temporarily store velocity of this time step (for accelaration calculation)
                ux0 = trc%now%ux(i)
                uy0 = trc%now%uy(i)
                uz0 = trc%now%uz(i) 

                ! Linear interpolation used for surface position 
                par_lin = interp_bilinear_weights(x1,y1,xout=trc%now%x(i),yout=trc%now%y(i))
                trc%now%H(i)     = interp_bilinear(par_lin,H1)
                trc%now%z_srf(i) = interp_bilinear(par_lin,z_srf1)
                trc%now%z(i)     = trc%now%z_srf(i) - trc%now%dpth(i)

                ! Calculate zc-axis for the current point
                ! (z_bedrock + ice thickness)
                ! Note: equivalent to (z_srf - depth) = trc%now%z_srf(i) - (1.0-z1)*trc%now%H(i)
                zc = (trc%now%z_srf(i)-trc%now%H(i)) + z1*trc%now%H(i)

                if (trim(trc%par%interp_method) .eq. "linear") then 
                    ! Trilinear interpolation 

                    ! Note: currently we redundantly obtain bilinear (horizontal) weights, because
                    ! they are needed to calculate zc. In the future, this could be improved. 

                    par_lin = interp_trilinear_weights(x1,y1,zc,xout=trc%now%x(i),yout=trc%now%y(i),zout=trc%now%z(i))

                    trc%now%ux(i)  = interp_trilinear(par_lin,ux1)
                    trc%now%uy(i)  = interp_trilinear(par_lin,uy1)
                    trc%now%uz(i)  = interp_trilinear(par_lin,uz1)

                else
                    ! Spline interpolation 
                    trc%now%ux(i) = interp_bspline3D(bspline3d_ux,trc%now%x(i),trc%now%y(i),trc%now%z(i)/trc%now%H(i))
                    trc%now%uy(i) = interp_bspline3D(bspline3d_uy,trc%now%x(i),trc%now%y(i),trc%now%z(i)/trc%now%H(i))
                    trc%now%uz(i) = interp_bspline3D(bspline3d_uz,trc%now%x(i),trc%now%y(i),trc%now%z(i)/trc%now%H(i)) *trc%now%H(i)  ! sigma => m

                end if 

                ! Update acceleration term 
                trc%now%ax(i) = (trc%now%ux(i) - ux0) / trc%now%dt
                trc%now%ay(i) = (trc%now%uy(i) - uy0) / trc%now%dt
                trc%now%az(i) = (trc%now%uz(i) - uz0) / trc%now%dt

                ! Filler values of the tracer state, in the future these should
                ! equal the surface temperature and the accumulation rate at the time of
                ! deposition and be calculated otherwise
                trc%now%T(i)   = 260.0 
                trc%now%thk(i) = 0.3 

            end if 

        end do 

        ! Update the tracer thickness, then destroy points that are too thin 
        ! == TO DO == 

        ! Update the tracer positions 
        call calc_position(trc%now%x,trc%now%y,trc%now%z,trc%now%ux,trc%now%uy,trc%now%uz, &
                           trc%now%ax,trc%now%ay,trc%now%az,trc%now%dt,trc%now%active)

!         call calc_position(trc%now%x,trc%now%y,trc%now%dpth,trc%now%ux,trc%now%uy,-trc%now%uz,trc%now%dt,trc%now%active)
        trc%now%dpth = max(trc%now%z_srf - trc%now%z, 0.0) 

        ! Destroy points that moved outside the valid region 
        call tracer_deactivate(trc,x1,y1,maxval(H1))

        ! Activate new tracers if desired
        if (dep_now) call tracer_activate(trc%par,trc%now,x1,y1,H=H1,lat=lat1, &
                                ux_srf=ux1(:,:,nx),uy_srf=uy1(:,:,ny),nmax=trc%par%n_max_dep)

        ! Finish activation for necessary points 
        do i = 1, trc%par%n 

            if (trc%now%active(i) .eq. 1) then 
                ! Point became active now, further initializations needed below

                ! Point is at the surface, so only bilinear interpolation is needed
                par_lin = interp_bilinear_weights(x1,y1,xout=trc%now%x(i),yout=trc%now%y(i))

                ! Apply interpolation weights to variables
                trc%now%dpth(i)  = 0.01   ! Always deposit just below the surface (eg 1 cm) to avoid zero z-velocity
                trc%now%z_srf(i) = interp_bilinear(par_lin,z_srf1)
                trc%now%z(i)     = trc%now%z_srf(i)-trc%now%dpth(i)
                
                trc%now%H(i)     = interp_bilinear(par_lin,H1)
                trc%now%ux(i)    = interp_bilinear(par_lin,ux1(:,:,nz))
                trc%now%uy(i)    = interp_bilinear(par_lin,uy1(:,:,nz))
                trc%now%uz(i)    = interp_bilinear(par_lin,uz1(:,:,nz)) 
                trc%now%ax(i)    = 0.0 
                trc%now%ay(i)    = 0.0 
                trc%now%az(i)    = 0.0
                
                ! Initialize state variables
                trc%now%T(i)   = 260.0 
                trc%now%thk(i) = 0.3 

                ! Define deposition values 
                trc%dep%time(i) = trc%now%time 
                trc%dep%H(i)    = trc%now%H(i)
                trc%dep%x(i)    = trc%now%x(i)
                trc%dep%y(i)    = trc%now%y(i)
                trc%dep%z(i)    = trc%now%z(i) 

                trc%dep%lon(i)  = interp_bilinear(par_lin,lon1)
                trc%dep%lat(i)  = interp_bilinear(par_lin,lat1)

                trc%dep%t2m_ann(i)   = interp_bilinear(par_lin,t2m_ann1)
                trc%dep%t2m_sum(i)   = interp_bilinear(par_lin,t2m_sum1)
                trc%dep%pr_ann(i)    = interp_bilinear(par_lin,pr_ann1)
                trc%dep%pr_sum(i)    = interp_bilinear(par_lin,pr_sum1)
                trc%dep%t2m_prann(i) = MV
                trc%dep%d18O_ann(i)  = MV

                trc%now%active(i) = 2 

            end if 

        end do 


        ! == TO DO == 
        ! - Attach whatever information we want to trace (age, deposition elevation and location, climate, isotopes, etc)
        ! - Potentially attach this information via a separate subroutine, using a flag to see if
        !   it was just deposited, then in main program calling eg, tracer_add_dep_variable(trc,"T",T),
        !   where the argument "T" should match a list of available variables, and T should be the variable
        !   to be stored from the main program. 
        !   Downside to above approach, is re-calculating the par_lin object every time. 

        ! Update summary statistics 
        trc%par%n_active = count(trc%now%active.gt.0)

        ! Calculate some summary information on eulerian grid if desired 
        if (stats_now) call calc_tracer_stats(trc,x,y,z,z_srf,H)

        return 

    end subroutine tracer_update

    subroutine tracer_end(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        
        ! Allocate the state variables 
        call tracer_deallocate(trc%now,trc%dep)

        write(*,*) "tracer:: tracer object deallocated."
        
        return 

    end subroutine tracer_end

    ! ================================================
    !
    ! tracer management 
    !
    ! ================================================
    
    subroutine tracer_activate(par,now,x,y,H,lat,ux_srf,uy_srf,nmax)
        ! Use this to activate individual or multiple tracers (not more than nmax)
        ! Only determine x/y position here, later interpolate z_srf and deposition
        ! information 

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: x(:), y(:)
        real(prec), intent(IN) :: H(:,:), lat(:,:), ux_srf(:,:), uy_srf(:,:) 
        integer, intent(IN) :: nmax  

        integer :: ntot  
        real(prec) :: p(size(H,1),size(H,2)), p_init(size(H,1),size(H,2))
        integer :: i, j, k, ij(2)
        real(prec), allocatable :: jit(:,:), dens(:,:)
        real(prec) :: xmin, ymin, xmax, ymax 

        ! How many points can be activated?
        ntot = min(nmax,count(now%active == 0))


        if (ntot .gt. 0) then 
            ! Proceed with activation, since points are available 

            ! Determine initial desired distribution of points on low resolution grid
!             p_init = gen_distribution_thickness(H,H_min=par%H_min_dep,alpha=par%alpha,dist=par%weight)
            p_init = gen_distribution_vel(uv=sqrt(ux_srf**2+uy_srf**2),H=H,uv_max=par%U_max_dep,H_min=par%H_min_dep)

!             ! Additionally adjust distribution according to latitude 
!             where (lat .lt. 70.0) 
!                 p_init = 0.0 
!             end where 

            ! Normalize p_init, just in case
            p_init = p_init / sum(p_init)

            ! Generate random numbers to populate points 
            allocate(jit(2,ntot))

            if (par%noise) then 
                call random_number(jit)
                jit = (jit - 0.5)
                jit(1,:) = jit(1,:)*(x(2)-x(1)) 

                if (size(y,1) .gt. 2) then 
                    jit(2,:) = jit(2,:)*(y(2)-y(1)) 
                else   ! Profile
                    jit(2,:) = 0.0 
                end if 
            else 
                jit = 0.0 
            end if 

    !         write(*,*) "range jit: ", minval(jit), maxval(jit)
    !         write(*,*) "npts: ", count(now%active == 0)
    !         write(*,*) "ntot: ", ntot 
    !         stop 
        
            ! Calculate domain boundaries to be able to apply limits 
            xmin = minval(x) 
            xmax = maxval(x) 
            ymin = minval(y) 
            ymax = maxval(y) 

            if (maxval(p_init) .gt. 0.0) then 
                ! Activate points in locations with non-zero probability
                ! This if-statement ensures some valid points currently exist in the domain
                
                k = 0 
                p = p_init   ! Set probability distribution to initial distribution 

                do j = 1, par%n 

                    if (now%active(j)==0) then 

                        now%active(j) = 1
                        k = k + 1
                        par%id_max = par%id_max+1 
                        now%id(j)  = par%id_max 

                        ij = maxloc(p,mask=p.gt.0.0)
                        now%x(j) = x(ij(1)) + jit(1,k)
                        now%y(j) = y(ij(2)) + jit(2,k)
                        
                        if (now%x(j) .lt. xmin) now%x(j) = xmin 
                        if (now%x(j) .gt. xmax) now%x(j) = xmax 
                        if (now%y(j) .lt. ymin) now%y(j) = ymin 
                        if (now%y(j) .gt. ymax) now%y(j) = ymax 
                        
                        p(ij(1),ij(2)) = 0.0 
                         
                    end if 

                    ! Stop when all points have been allocated
                    if (k .ge. ntot) exit 

                    ! If there are no more points with non-zero probability, reset probability
                    if (maxval(p) .eq. 0.0) p = p_init 

                end do 
              
            end if

            ! Summary 
    !         write(*,*) "tracer_activate:: ", count(now%active == 0), count(now%active .eq. 1), count(now%active .eq. 2)

        end if 

        return 

    end subroutine tracer_activate 

    subroutine tracer_deactivate(trc,x,y,Hmax)
        ! Use this to deactivate individual or multiple tracers
        implicit none 

        type(tracer_class),   intent(INOUT) :: trc  
        real(prec), intent(IN) :: x(:), y(:) 
        real(prec), intent(IN) :: Hmax 

        ! Deactivate points where:
        !  - Thickness of ice sheet at point's location is below threshold
        !  - Point is above maximum ice thickness Hmax (interp error)
        !  - Point is past maximum depth into the ice sheet 
        !  - Velocity of point is higher than maximum threshold 
        !  - x/y position is out of boundaries of the domain 
        where (trc%now%active .gt. 0 .and. &
              ( trc%now%H .lt. trc%par%H_min                            .or. &
                trc%now%H .gt. Hmax                                     .or. &
                trc%now%dpth/trc%now%H .ge. trc%par%depth_max           .or. &
                sqrt(trc%now%ux**2 + trc%now%uy**2) .gt. trc%par%U_max  .or. &
                trc%now%x .lt. minval(x) .or. trc%now%x .gt. maxval(x)  .or. &
                trc%now%y .lt. minval(y) .or. trc%now%y .gt. maxval(y) ) ) 

            trc%now%active    = 0 

            trc%now%id        = mv 
            trc%now%x         = mv 
            trc%now%y         = mv 
            trc%now%z         = mv 
            trc%now%sigma     = mv 
            trc%now%z_srf     = mv 
            trc%now%dpth      = mv 
            trc%now%ux        = mv 
            trc%now%uy        = mv 
            trc%now%uz        = mv 
            trc%now%thk       = mv 
            trc%now%T         = mv 
            trc%now%H         = mv 

            trc%dep%time      = mv 
            trc%dep%H         = mv 
            trc%dep%x         = mv 
            trc%dep%y         = mv 
            trc%dep%z         = mv 

        end where 

        return 

    end subroutine tracer_deactivate 

    ! ================================================
    !
    ! tracer physics / stats
    !
    ! ================================================
    
    elemental subroutine calc_position(x,y,z,ux,uy,uz,ax,ay,az,dt,active)

        implicit none 

        real(prec),   intent(INOUT) :: x, y, z 
        real(prec),   intent(IN)    :: ux, uy, uz 
        real(prec),   intent(IN)    :: ax, ay, az 
        real(prec_time), intent(IN)    :: dt 
        integer,         intent(IN)    :: active 

        if (active .gt. 0) then 

            x = x + ux*dt + 0.5*ax*dt**2 
            y = y + uy*dt + 0.5*ay*dt**2 
            z = z + uz*dt + 0.5*az*dt**2
            
        end if 

        return 

    end subroutine calc_position

    function gen_distribution_vel(uv,H,uv_max,H_min) result(p)

        implicit none 

        real(prec), intent(IN) :: uv(:,:), H(:,:)
        real(prec), intent(IN) :: uv_max, H_min 
        real(prec) :: p(size(H,1),size(H,2))

        ! Local variables
        real(prec) :: p_sum 

        
        p = 0.0 
        where (uv .gt. 0.0 .and. uv .lt. uv_max .and. H .gt. H_min)
            p = 1.0 - uv/uv_max 
        end where 

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution_vel

    function gen_distribution_thickness(H,H_min,alpha,dist) result(p)

        implicit none 

        real(prec), intent(IN) :: H(:,:)
        real(prec), intent(IN) :: H_min, alpha 
        character(len=*), intent(IN) :: dist 
        real(prec) :: p(size(H,1),size(H,2))

        ! Local variables
        integer    :: k, ij(2)
        real(prec) :: p_sum 

        select case(trim(dist)) 

            case("linear")

                p = (alpha * max(H-H_min,0.0) / (maxval(H)-H_min))

            case("quadratic")

                p = (alpha * max(H-H_min,0.0) / (maxval(H)-H_min))**2

            case DEFAULT   ! "rand"

                ! Random even distribution (all points equally likely)
                call random_number(p)
                where (H .lt. H_min) p = 0.0 


        end select 

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution_thickness

    function gen_distribution_direction(x,y,u,v,theta_max) result(p)

        implicit none 

        real(prec), intent(IN) :: x(:), y(:), u(:,:), v(:,:) 
        real(prec), intent(IN) :: theta_max 
        real(prec) :: p(size(u,1),size(u,2))

        ! Local variables
        integer    :: k, ij(2)
        real(prec) :: p_sum 

        p = 1.0 

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution_direction

    elemental function calc_angle(x1,y1,x2,y2) result(theta)
        ! Given a vector, calculate the angle wrt unit circle 

        implicit none 

        real(prec), intent(IN) :: x1, y1, x2, y2 
        real(prec) :: theta 

        theta = atan2((y2-y1),(x2-x1))

        return 

    end function calc_angle 

    subroutine calc_tracer_stats(trc,x,y,z,z_srf,H)
        ! Convert tracer information to isochrone format matching
        ! Macgregor et al. (2015)
        ! Note: this should only be called at time t=0 ka BP, 
        ! since age is defined assuming that. 

        implicit none
        
        type(tracer_class), intent(INOUT) :: trc
        real(prec), intent(IN) :: x(:), y(:), z(:), z_srf(:,:), H(:,:)

        ! Local variables 
        integer :: i, j, k, q
        integer :: nx, ny, nz, nq
        real(prec), allocatable :: dx(:), dy(:) 
        real(prec) :: dt, dz, dz_now  
        real(prec) :: zc(size(z))
        integer       :: id(trc%par%n)
        real(prec)    :: dist(trc%par%n)
        integer, allocatable :: inds(:)
        integer :: n_ind 

        nx = size(x)
        ny = size(y)
        
        allocate(dx(nx+1),dy(ny+1))
        dx(2:nx) = (x(2:nx)-x(1:nx-1))/2.0
        dx(1)    = dx(2)
        dx(nx+1) = dx(nx)
        dy(2:ny) = (y(2:ny)-y(1:ny-1))/2.0
        dy(1)    = dy(2)
        dy(ny+1) = dy(ny)
        
        dt = 5.0 ! Isochrone uncertainty of ±5 ka  
        dz = (trc%stats%depth_norm(2) - trc%stats%depth_norm(1))/2.0   ! depth_norm is equally spaced

        ! Loop over grid and fill in information
        do j = 1, ny 
        do i = 1, nx 

            ! Calculate the isochrones
            nq = size(trc%stats%age_iso)
            do q = 1, nq 

                ! Filter for active particles within the grid box and age of interest
                call which (trc%now%active == 2 .and. &
                            trc%now%x .gt. x(i)-dx(i) .and. trc%now%x .le. x(i)+dx(i+1) .and. &
                            trc%now%y .gt. y(j)-dy(j) .and. trc%now%y .le. y(j)+dy(j+1) .and. &
                            (0.0 - trc%dep%time)*1e-3 .ge. trc%stats%age_iso(q)-dt .and. &
                            (0.0 - trc%dep%time)*1e-3 .le. trc%stats%age_iso(q)+dt, inds, n_ind) 

                ! Calculate range mean/sd depth for given age range
                if (n_ind .gt. 0) then
                    write(*,*) "isochrones: ", i, j, q, n_ind 
 
                    trc%stats%depth_iso(i,j,q)     = calc_mean(real(trc%now%dpth(inds),prec_wrt))
                    trc%stats%depth_iso_err(i,j,q) = calc_sd(real(trc%now%dpth(inds),prec_wrt), &
                                                             trc%stats%depth_iso(i,j,q))
                    trc%stats%density_iso(i,j,q)   = n_ind 

                    trc%stats%dep_z_iso(i,j,q)     = calc_mean(real(trc%dep%z(inds),prec_wrt))
                else 
                    trc%stats%depth_iso(i,j,q)     = MV 
                    trc%stats%depth_iso_err(i,j,q) = MV
                    trc%stats%density_iso(i,j,q)   = MV
                    
                    trc%stats%dep_z_iso(i,j,q)     = MV 

                end if 

            end do 
            
            ! Calculate the ages of each depth layer 
            nq = size(trc%stats%depth_norm)
            do q = 1, nq 

                ! Use dz_now to ensure that the first depth (0.04) includes all depths to the surface
                dz_now = dz 
                if (q .eq. 1) dz_now = dz*5.0 

                ! Filter for active particles within the grid box and age of interest
                call which (trc%now%active == 2 .and. &
                            trc%now%x .gt. x(i)-dx(i) .and. trc%now%x .le. x(i)+dx(i+1) .and. &
                            trc%now%y .gt. y(j)-dy(j) .and. trc%now%y .le. y(j)+dy(j+1) .and. &
                            trc%now%dpth/trc%now%H .gt. trc%stats%depth_norm(q)-dz_now  .and. &
                            trc%now%dpth/trc%now%H .le. trc%stats%depth_norm(q)+dz, inds, n_ind) 

                ! Calculate range mean/sd age for given depth range
                if (n_ind .gt. 0) then
                    write(*,*) "ice_ages: ", i, j, q, n_ind 
 
                    trc%stats%ice_age(i,j,q)     = calc_mean(real(0.0-trc%dep%time(inds),prec_wrt))*1e-3
                    trc%stats%ice_age_err(i,j,q) = calc_sd(real(0.0-trc%dep%time(inds),prec_wrt),trc%stats%ice_age(i,j,q))*1e-3
                    trc%stats%density(i,j,q)     = n_ind 
                else 
                    trc%stats%ice_age(i,j,q)     = MV 
                    trc%stats%ice_age_err(i,j,q) = MV 
                    trc%stats%density(i,j,q)     = MV 
                end if  

            end do 

        end do 
        end do  
        
        return

    end subroutine calc_tracer_stats

    function calc_mean(x) result(mean)

        implicit none 

        real(prec_wrt), intent(IN) :: x(:) 
        real(prec_wrt) :: mean 
        integer :: n 

        n = count(x.ne.MV)

        if (n .gt. 0) then 
            mean = sum(x,mask=x.ne.MV) / real(n)
        else 
            mean = MV 
        end if 
        
        return 

    end function calc_mean 

    function calc_sd(x,mean) result(stdev)

        implicit none 

        real(prec_wrt), intent(IN) :: x(:) 
        real(prec_wrt) :: mean 
        real(prec_wrt) :: stdev 
        integer :: n 

        n = count(x.ne.MV)

        if (n .gt. 0) then 
            stdev = sqrt( sum((x - mean)**2) / real(n) )
        else 
            stdev = MV 
        end if 

        return 

    end function calc_sd 

    ! ================================================
    !
    ! Initialization routines 
    !
    ! ================================================

    subroutine tracer_par_load(par,filename,is_sigma)

        implicit none 

        type(tracer_par_class), intent(OUT) :: par 
        character(len=*),       intent(IN)  :: filename 
        logical, intent(IN) :: is_sigma 

!         par%n         = 5000
!         par%n_max_dep = 100
        
!         par%thk_min   = 1e-2   ! m (minimum thickness of tracer 'layer')
!         par%H_min     = 1500.0 ! m 
!         par%depth_max = 0.99   ! fraction of thickness
!         par%U_max     = 200.0  ! m/a 

!         par%H_min_dep = 1000.0 ! m 
!         par%alpha     = 1.0 
!         par%weight    = "linear"
        
!         par%dens_z_lim = 50.0 ! m
!         par%dens_max   = 10   ! Number of points
        
        call nml_read(filename,"tracer_par","dt",            par%dt)
        call nml_read(filename,"tracer_par","n",             par%n)
        call nml_read(filename,"tracer_par","n_max_dep",     par%n_max_dep)
        call nml_read(filename,"tracer_par","dt_dep",        par%dt_dep)
        call nml_read(filename,"tracer_par","dt_write",      par%dt_write)
        call nml_read(filename,"tracer_par","thk_min",       par%thk_min)
        call nml_read(filename,"tracer_par","H_min",         par%H_min)
        call nml_read(filename,"tracer_par","depth_max",     par%depth_max)
        call nml_read(filename,"tracer_par","U_max",         par%U_max)
        call nml_read(filename,"tracer_par","U_max_dep",     par%U_max_dep)
        call nml_read(filename,"tracer_par","H_min_dep",     par%H_min_dep)
        call nml_read(filename,"tracer_par","alpha",         par%alpha)
        call nml_read(filename,"tracer_par","weight",        par%weight)
        call nml_read(filename,"tracer_par","noise",         par%noise)
        call nml_read(filename,"tracer_par","dens_z_lim",    par%dens_z_lim)
        call nml_read(filename,"tracer_par","dens_max",      par%dens_max)
        call nml_read(filename,"tracer_par","interp_method", par%interp_method)
        call nml_read(filename,"tracer_par","par_trans_file",par%par_trans_file)
    
        ! Define additional parameter values
        par%is_sigma  = is_sigma 
        par%n_active  = 0 

        par%use_par_trans = .FALSE.
        if (trim(par%par_trans_file) .ne. "None") then 
            par%use_par_trans = .TRUE. 

            call tracer_par_trans_load(par%tpar,par%par_trans_file)
        end if 

        ! Consistency checks 
        if (trim(par%interp_method) .ne. "linear" .and. &
            trim(par%interp_method) .ne. "spline" ) then 
            write(0,*) "tracer_init:: error: interp_method must be 'linear' &
            &or 'spline': "//trim(par%interp_method)
            stop 
        end if 


        return 

    end subroutine tracer_par_load
    
    subroutine tracer_par_update(par,tpar,time)
        ! Update transient parameter values for current time 

        implicit none 

        type(tracer_par_class), intent(INOUT) :: par 
        type(tracer_par_trans_class), intent(IN) :: tpar 
        real(prec_time) :: time 

        ! Local variables 
        integer :: i, n, k  

        n = size(tpar%time)

        ! Initially assume the first row is correct 
        k = 1
        
        ! Check to see if one of following rows is correct, update k 
        do i = 2, n 
            if (tpar%time(i) .gt. time) exit 
            k = k+1
        end do 

        par%H_min_dep = tpar%H_min_dep(k) 
        par%dt_dep    = tpar%dt_dep(k)
        par%n_max_dep = tpar%n_max_dep(k) 
        par%dt_write  = tpar%dt_write(k) 

        return 

    end subroutine tracer_par_update

    subroutine tracer_par_trans_load(tpar,filename)
        ! This subroutine will read a time series of
        ! several columns [time,par1,par2,...,parN] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        type(tracer_par_trans_class), intent(OUT) :: tpar 
        character(len=*), intent(IN) :: filename 

        ! Local variables 
        integer, parameter :: f = 191
        integer, parameter :: nmax = 10000

        integer :: i, stat, n 
        character(len=256) :: str, str1 
        real(4) :: x(nmax), y1(nmax), y2(nmax), y3(nmax), y4(nmax)

        ! Open file for reading 
        open(f,file=filename,status="old")

        ! Read the header in the first line: 
        read(f,*,IOSTAT=stat) str

        n = 0 

        do i = 1, nmax 
            read(f,'(a100)',IOSTAT=stat) str 

            ! Exit loop if the end-of-file is reached 
            if(IS_IOSTAT_END(stat)) exit 

            str1 = adjustl(trim(str))
!            str1=str
            if ( len(trim(str1)) .gt. 0 ) then 
                if ( .not. (str1(1:1) == "!" .or. &
                            str1(1:1) == "#") ) then 
                    n = n+1
                    read(str1,*) x(n), y1(n), y2(n), y3(n), y4(n)
                end if
            end if  
        end do 

        ! Close the file
        close(f) 

        if (n .eq. nmax) then 
            write(*,*) "tracer_par_trans_load:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data 
        call tracer_par_trans_allocate(tpar,n)

        tpar%time      =  x(1:n) 
        tpar%H_min_dep = y1(1:n) 
        tpar%dt_dep    = y2(1:n) 
        tpar%n_max_dep = y3(1:n) 
        tpar%dt_write  = y4(1:n) 

        write(*,*) "tracer_par_trans_load:: Time series read from file: "//trim(filename)
        write(*,"(a12,4a10)") "time", "H_min_dep", "dt_dep", "n_max_dep", "dt_write"
        do i = 1, n 
            write(*,"(g12.3,f10.1,f10.1,i10,f10.1)") tpar%time(i), tpar%H_min_dep(i), tpar%dt_dep(i), &
                       tpar%n_max_dep(i), tpar%dt_write(i) 
        end do

        return 

    end subroutine tracer_par_trans_load 

    subroutine tracer_par_trans_allocate(tpar,n)

        implicit none 

        type(tracer_par_trans_class), intent(INOUT) :: tpar 
        integer, intent(IN) :: n 

        ! Make sure all arrays are deallocated first 
        if (allocated(tpar%time))      deallocate(tpar%time)
        if (allocated(tpar%H_min_dep)) deallocate(tpar%H_min_dep)
        if (allocated(tpar%dt_dep))    deallocate(tpar%dt_dep)
        if (allocated(tpar%n_max_dep)) deallocate(tpar%n_max_dep)
        if (allocated(tpar%dt_write))  deallocate(tpar%dt_write)
        
        allocate(tpar%time(n))
        allocate(tpar%H_min_dep(n))
        allocate(tpar%dt_dep(n))
        allocate(tpar%n_max_dep(n))
        allocate(tpar%dt_write(n))
        
        return 

    end subroutine tracer_par_trans_allocate

    subroutine tracer_allocate(now,dep,n)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        integer, intent(IN) :: n

        ! Make object is deallocated
        call tracer_deallocate(now,dep)

        ! Allocate tracer 
        allocate(now%active(n))
        allocate(now%id(n))
        allocate(now%x(n),now%y(n),now%z(n),now%sigma(n))
        allocate(now%z_srf(n),now%dpth(n))
        allocate(now%ux(n),now%uy(n),now%uz(n))
        allocate(now%ax(n),now%ay(n),now%az(n))
        allocate(now%thk(n))
        allocate(now%T(n))
        allocate(now%H(n))

        ! Allocate deposition properties 
        
        allocate(dep%time(n), dep%H(n))
        allocate(dep%x(n), dep%y(n), dep%z(n), dep%lon(n), dep%lat(n))
        allocate(dep%t2m_ann(n), dep%t2m_sum(n), dep%pr_ann(n), dep%pr_sum(n),dep%t2m_prann(n))
        allocate(dep%d18O_ann(n))
        
        return

    end subroutine tracer_allocate

    subroutine tracer_deallocate(now,dep)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        
        ! Deallocate state objects
        if (allocated(now%active))    deallocate(now%active)
        if (allocated(now%x))         deallocate(now%x)
        if (allocated(now%y))         deallocate(now%y)
        if (allocated(now%z))         deallocate(now%z)
        if (allocated(now%sigma))     deallocate(now%sigma)
        if (allocated(now%z_srf))     deallocate(now%z_srf)
        if (allocated(now%dpth))      deallocate(now%dpth)
        if (allocated(now%ux))        deallocate(now%ux)
        if (allocated(now%uy))        deallocate(now%uy)
        if (allocated(now%uz))        deallocate(now%uz)
        if (allocated(now%ax))        deallocate(now%ax)
        if (allocated(now%ay))        deallocate(now%ay)
        if (allocated(now%az))        deallocate(now%az)
        if (allocated(now%thk))       deallocate(now%thk)
        if (allocated(now%T))         deallocate(now%T)
        if (allocated(now%H))         deallocate(now%H)

        ! Deallocate deposition objects
        if (allocated(dep%time))      deallocate(dep%time)
        if (allocated(dep%z))         deallocate(dep%z)
        if (allocated(dep%H))         deallocate(dep%H)
        if (allocated(dep%x))         deallocate(dep%x)
        if (allocated(dep%y))         deallocate(dep%y)
        if (allocated(dep%lon))       deallocate(dep%lon)
        if (allocated(dep%lat))       deallocate(dep%lat)
        if (allocated(dep%t2m_ann))   deallocate(dep%t2m_ann)
        if (allocated(dep%t2m_sum))   deallocate(dep%t2m_sum)
        if (allocated(dep%pr_ann))    deallocate(dep%pr_ann)
        if (allocated(dep%pr_sum))    deallocate(dep%pr_sum)
        if (allocated(dep%t2m_prann)) deallocate(dep%t2m_prann)
        if (allocated(dep%d18O_ann))  deallocate(dep%d18O_ann)
        
        return

    end subroutine tracer_deallocate

    subroutine tracer_allocate_stats(stats,x,y)

        implicit none 

        type(tracer_stats_class), intent(INOUT) :: stats 
        real(prec), intent(IN) :: x(:), y(:)
        
        ! Make surce object is deallocated
        call tracer_deallocate_stats(stats)

        ! Allocate tracer stats axes
        allocate(stats%x(size(x)))
        allocate(stats%y(size(y)))

        ! Allocate tracer stats objects
        allocate(stats%depth_norm(25))  ! To match Macgregor et al. (2015)
        allocate(stats%age_iso(5))      ! To match Macgregor et al. (2015)
        allocate(stats%depth_iso(size(x),size(y),size(stats%age_iso)))
        allocate(stats%depth_iso_err(size(x),size(y),size(stats%age_iso)))
        allocate(stats%dep_z_iso(size(x),size(y),size(stats%age_iso)))
        allocate(stats%density_iso(size(x),size(y),size(stats%age_iso)))
        allocate(stats%ice_age(size(x),size(y),size(stats%depth_norm)))
        allocate(stats%ice_age_err(size(x),size(y),size(stats%depth_norm)))
        allocate(stats%density(size(x),size(y),size(stats%depth_norm)))
        
        ! Also store axis information directly
        stats%x = x 
        stats%y = y 

        ! Initialize arrays to zeros 
        stats%depth_iso     = 0.0 
        stats%depth_iso_err = 0.0
        stats%dep_z_iso     = MV  
        stats%density_iso   = MV 
        stats%ice_age       = 0.0 
        stats%ice_age_err   = 0.0 
        stats%density       = MV 
        
        return

    end subroutine tracer_allocate_stats

    subroutine tracer_deallocate_stats(stats)

        implicit none 

        type(tracer_stats_class), intent(INOUT) :: stats 

        ! Deallocate stats objects
        if (allocated(stats%x))             deallocate(stats%x)
        if (allocated(stats%y))             deallocate(stats%y)
        if (allocated(stats%depth_norm))    deallocate(stats%depth_norm)
        if (allocated(stats%age_iso))       deallocate(stats%age_iso)
        if (allocated(stats%depth_iso))     deallocate(stats%depth_iso)
        if (allocated(stats%depth_iso_err)) deallocate(stats%depth_iso_err)
        if (allocated(stats%dep_z_iso))     deallocate(stats%dep_z_iso)
        if (allocated(stats%density_iso))   deallocate(stats%density_iso)
        if (allocated(stats%ice_age))       deallocate(stats%ice_age)
        if (allocated(stats%ice_age_err))   deallocate(stats%ice_age_err)
        if (allocated(stats%density))       deallocate(stats%density)
        
        return

    end subroutine tracer_deallocate_stats

    subroutine tracer_reshape1D_vec(var,var1,rev)

        implicit none 
     
        real(prec),    intent(IN) :: var(:)
        real(prec), intent(INOUT), allocatable :: var1(:)
        logical,    intent(IN) :: rev 

        integer :: i, nx

        nx = size(var,1)
        if (allocated(var1)) deallocate(var1)
        allocate(var1(nx))

        if (rev) then 
            do i = 1, nx
                var1(i) = var(nx-i+1)
            end do 
        else 
            var1 = var 
        end if 

        return 

    end subroutine tracer_reshape1D_vec

    subroutine tracer_reshape2D_field(idx_order,var,var1)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(prec),    intent(IN) :: var(:,:)
        real(prec), intent(INOUT), allocatable :: var1(:,:)
        integer :: i, j
        integer :: nx, ny

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

                nx = size(var,1)
                ny = size(var,2)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny))

                var1 = var 

            case("kji")
                ! z, y, x array order 

                nx = size(var,2)
                ny = size(var,1)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny))

                do i = 1, nx 
                do j = 1, ny 
                    var1(i,j)  = var(j,i)
                end do 
                end do 

            case DEFAULT 

                write(0,*) "tracer_reshape2D_field:: error: unrecognized array order: ",trim(idx_order)
                write(0,*) "    Possible choices are: ijk, kji"
                stop  

        end select 

        return 

    end subroutine tracer_reshape2D_field

    subroutine tracer_reshape3D_field(idx_order,var,var1,rev_z)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(prec),    intent(IN) :: var(:,:,:)
        real(prec), intent(INOUT), allocatable :: var1(:,:,:)
        logical,    intent(IN) :: rev_z   ! Reverse the z-axis? 
        integer :: i, j, k
        integer :: nx, ny, nz 

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

                nx = size(var,1)
                ny = size(var,2)
                nz = size(var,3)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny,nz))

                if (rev_z) then 
                    do i = 1, nx 
                    do j = 1, ny  
                    do k = 1, nz 
                        var1(i,j,k)  = var(i,j,nz-k+1) 
                    end do 
                    end do 
                    end do 
                else 
                    var1 = var 
                end if 

            case("kji")
                ! z, y, x array order 

                nx = size(var,3)
                ny = size(var,2)
                nz = size(var,1)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny,nz))

                if (rev_z) then 
                    do i = 1, nx 
                    do j = 1, ny 
                    do k = 1, nz 
                        var1(i,j,k)  = var(nz-k+1,j,i)
                    end do 
                    end do 
                    end do 
                else 
                    do i = 1, nx 
                    do j = 1, ny 
                    do k = 1, nz 
                        var1(i,j,k)  = var(k,j,i) 
                    end do 
                    end do 
                    end do 
                end if 

            case DEFAULT 

                write(0,*) "tracer_reshape3D_field:: error: unrecognized array order: ",trim(idx_order)
                write(0,*) "    Possible choices are: ijk, kji"
                stop  

        end select 

        return 

    end subroutine tracer_reshape3D_field

    subroutine which(x,ind,stat)
        ! Analagous to R::which function
        ! Returns indices that match condition x==.TRUE.

        implicit none 

        logical :: x(:)
        integer, allocatable :: tmp(:), ind(:)
        integer, optional :: stat  
        integer :: n, i  

        n = count(x)
        allocate(tmp(n))
        tmp = 0 

        n = 0
        do i = 1, size(x) 
            if (x(i)) then 
                n = n+1
                tmp(n) = i 
            end if
        end do 

        if (present(stat)) stat = n 

        if (allocated(ind)) deallocate(ind)

        if (n .eq. 0) then 
            allocate(ind(1))
            ind = -1 
        else
            allocate(ind(n))
            ind = tmp(1:n)
        end if 
        
        return 

    end subroutine which

end module tracer3D 


