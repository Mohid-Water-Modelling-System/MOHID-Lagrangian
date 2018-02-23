
module tracer_io 

    use ncio 

    use tracer_precision
    use tracer_interp 
    use tracer3D 

    implicit none 

    private 
    public :: tracer_write_init, tracer2D_write_init
    public :: tracer_write, tracer2D_write 
    public :: tracer_write_stats, tracer2D_write_stats
    public :: tracer_read 
    public :: tracer_align 
    public :: tracer_import_eulerian

contains 

    ! ================================================
    !
    ! I/O routines 
    !
    ! ================================================

    subroutine tracer_write_init(trc,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"pt",x=1,dx=1,nx=trc%par%n)
        call nc_write_dim(path_out,"time",x=real(mv,prec_wrt),unlimited=.TRUE.)

        return 

    end subroutine tracer_write_init 

    subroutine tracer_write(trc,time,fldr,filename,is2D)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time) :: time 
        character(len=*), intent(IN) :: fldr, filename 
        logical, intent(IN), optional :: is2D 

        ! Local variables 
        integer :: nt
        integer, allocatable :: dims(:)
        real(prec_wrt) :: time_in, mv_wrt   
        real(prec_wrt) :: tmp(size(trc%now%x))
        character(len=512) :: path_out 
        logical :: is_2D 

        trc%now%time_write = time 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

        ! Determine whether just writing a profile 
        is_2D = .FALSE. 
        if (present(is2D)) is_2D = is2D 

        ! Determine which timestep this is
        call nc_dims(path_out,"time",dims=dims)
        nt = dims(1)
        call nc_read(path_out,"time",time_in,start=[nt],count=[1])
        if (time_in .ne. MV .and. abs(time-time_in).gt.1e-2) nt = nt+1 

        call nc_write(path_out,"time",real(time,prec_wrt), dim1="time",start=[nt],count=[1],missing_value=mv_wrt)
        call nc_write(path_out,"n_active",trc%par%n_active,dim1="time",start=[nt],count=[1],missing_value=int(mv_wrt))
        
        tmp = trc%now%x
        where(trc%now%x .ne. mv_wrt) tmp = trc%now%x*1e-3
        call nc_write(path_out,"x",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="km")

        if (.not. is_2D) then 
            tmp = trc%now%y
            where(trc%now%y .ne. mv_wrt) tmp = trc%now%y*1e-3
            call nc_write(path_out,"y",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                            start=[1,nt],count=[trc%par%n ,1],units="km")
        end if 
        call nc_write(path_out,"z",real(trc%now%z,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"dpth",real(trc%now%dpth,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"z_srf",real(trc%now%z_srf,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"ux",real(trc%now%ux,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m/a")
        call nc_write(path_out,"uy",real(trc%now%uy,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m/a")
        call nc_write(path_out,"uz",real(trc%now%uz,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m/a")
        call nc_write(path_out,"thk",real(trc%now%thk,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"T",real(trc%now%T,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"H",real(trc%now%H,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")

        call nc_write(path_out,"id",trc%now%id,dim1="pt",dim2="time", missing_value=int(mv_wrt), &
                        start=[1,nt],count=[trc%par%n ,1])

        tmp = mv_wrt
        where(trc%dep%time .ne. mv_wrt) tmp = time-trc%dep%time
        call nc_write(path_out,"age",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="a")

        ! Write deposition information
        call nc_write(path_out,"dep_time",real(trc%dep%time,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="years")
        call nc_write(path_out,"dep_H",real(trc%dep%H,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        tmp = trc%dep%x
        where(trc%dep%x .ne. mv_wrt) tmp = trc%dep%x*1e-3
        call nc_write(path_out,"dep_x",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="km")
        
        if (.not. is_2D) then 
            tmp = trc%dep%y
            where(trc%dep%y .ne. mv_wrt) tmp = trc%dep%y*1e-3
            call nc_write(path_out,"dep_y",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                            start=[1,nt],count=[trc%par%n ,1],units="km")
        end if 
        call nc_write(path_out,"dep_z",real(trc%dep%z,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")

        return 

    end subroutine tracer_write 

    subroutine tracer_write_slice(trc,time,fldr,filename,is2D)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time) :: time 
        character(len=*), intent(IN) :: fldr, filename 
        logical, intent(IN), optional :: is2D 

        ! Local variables 
        real(prec_wrt) :: time_in, mv_wrt   
        real(prec_wrt), allocatable :: tmp(:)
        integer, allocatable :: inds(:)
        character(len=512) :: path_out 
        logical :: is_2D 

        trc%now%time_write = time 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

        ! Determine whether just writing a profile 
        is_2D = .FALSE. 
        if (present(is2D)) is_2D = is2D 

        if (trc%par%n_active .eq. 0) then 
            ! Do not write slice!!


        else 

            ! Create output file 
            call nc_create(path_out)
            call nc_write_dim(path_out,"pt",x=1,dx=1,nx=trc%par%n_active)
            call nc_write_dim(path_out,"time",x=real(time_in,prec_wrt),unlimited=.TRUE.)

            allocate(tmp(trc%par%n_active))

            ! Get indices of active particles 
            call which(trc%now%active .eq. 2, inds)

            tmp = trc%now%x(inds)*1e-3
            call nc_write(path_out,"x",tmp,dim1="pt", missing_value=mv_wrt,units="km")

            if (.not. is_2D) then 
                tmp = trc%now%y(inds)*1e-3
                call nc_write(path_out,"y",tmp,dim1="pt", missing_value=mv_wrt,units="km")
            end if 

            ! TO DO 
!             call nc_write(path_out,"z",real(trc%now%z,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")
!             call nc_write(path_out,"dpth",real(trc%now%dpth,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")
!             call nc_write(path_out,"z_srf",real(trc%now%z_srf,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")
!             call nc_write(path_out,"ux",real(trc%now%ux,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m/a")
!             call nc_write(path_out,"uy",real(trc%now%uy,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m/a")
!             call nc_write(path_out,"uz",real(trc%now%uz,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m/a")
!             call nc_write(path_out,"thk",real(trc%now%thk,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")
!             call nc_write(path_out,"T",real(trc%now%T,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1])
!             call nc_write(path_out,"H",real(trc%now%H,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")

!             call nc_write(path_out,"id",trc%now%id,dim1="pt",dim2="time", missing_value=int(mv_wrt), &
!                             start=[1,nt],count=[trc%par%n ,1])

!             tmp = mv_wrt
!             where(trc%dep%time .ne. mv_wrt) tmp = time-trc%dep%time
!             call nc_write(path_out,"age",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="a")

!             ! Write deposition information
!             call nc_write(path_out,"dep_time",real(trc%dep%time,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="years")
!             call nc_write(path_out,"dep_H",real(trc%dep%H,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")
!             tmp = trc%dep%x
!             where(trc%dep%x .ne. mv_wrt) tmp = trc%dep%x*1e-3
!             call nc_write(path_out,"dep_x",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="km")
            
!             if (.not. is_2D) then 
!                 tmp = trc%dep%y
!                 where(trc%dep%y .ne. mv_wrt) tmp = trc%dep%y*1e-3
!                 call nc_write(path_out,"dep_y",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
!                                 start=[1,nt],count=[trc%par%n ,1],units="km")
!             end if 
!             call nc_write(path_out,"dep_z",real(trc%dep%z,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
!                             start=[1,nt],count=[trc%par%n ,1],units="m")

        end if 

        return 

    end subroutine tracer_write_slice 

    subroutine tracer_write_diagnostic_init(trc,time,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec_time) :: time
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",        x=trc%stats%x*1e-3,     units="km")
        call nc_write_dim(path_out,"yc",        x=trc%stats%y*1e-3,     units="km")
        call nc_write_dim(path_out,"depth_norm",x=trc%stats%depth_norm, units="1")
        call nc_write_dim(path_out,"age_iso",   x=trc%stats%age_iso,    units="ka")
        call nc_write_dim(path_out,"time",      x=time,unlimited=.TRUE.,units="ka")
        
        return 

    end subroutine tracer_write_diagnostic_init 

    subroutine tracer_write_diagnostic_stats(trc,time,fldr,filename) !,z_srf,H)
        ! Write various meta-tracer information (ie, lagrangian => eulerian)
        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec_time) :: time
        character(len=*),   intent(IN) :: fldr, filename 
!         real(prec),         intent(IN) :: z_srf(:,:), H(:,:) 

        ! Local variables 
        character(len=512) :: path_out 
        real(prec_wrt) :: mv_wrt 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

!         call nc_write(path_out,"z_srf",z_srf,dim1="xc",dim2="yc",missing_value=mv_wrt, &
!                       units="m",long_name="Surface elevation")
!         call nc_write(path_out,"H",H,dim1="xc",dim2="yc",missing_value=mv_wrt, &
!                       units="m",long_name="Ice thickness")

        call nc_write(path_out,"ice_age",trc%stats%ice_age,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age")
        call nc_write(path_out,"ice_age_err",trc%stats%ice_age_err,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age - error")
        call nc_write(path_out,"density",trc%stats%density,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density")

        call nc_write(path_out,"depth_iso",trc%stats%depth_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone depth")
        call nc_write(path_out,"depth_iso_err",trc%stats%depth_iso_err,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone depth - error")
        call nc_write(path_out,"dep_z_iso",trc%stats%dep_z_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone deposition elevation")
        call nc_write(path_out,"density_iso",trc%stats%density_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density (for isochrones)")

        
        return 

    end subroutine tracer_write_diagnostic_stats
    
    subroutine tracer_write_stats(trc,time,fldr,filename) !,z_srf,H)
        ! Write various meta-tracer information (ie, lagrangian => eulerian)
        ! This output belongs to a specific time slice, usually at time = 0 ka BP. 

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec_time) :: time
        character(len=*),   intent(IN) :: fldr, filename 
!         real(prec),         intent(IN) :: z_srf(:,:), H(:,:) 

        ! Local variables 
        character(len=512) :: path_out 
        real(prec_wrt) :: mv_wrt 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",        x=trc%stats%x*1e-3,     units="km")
        call nc_write_dim(path_out,"yc",        x=trc%stats%y*1e-3,     units="km")
        call nc_write_dim(path_out,"depth_norm",x=trc%stats%depth_norm, units="1")
        call nc_write_dim(path_out,"age_iso",   x=trc%stats%age_iso,    units="ka")
        call nc_write_dim(path_out,"time",      x=time,unlimited=.TRUE.,units="ka")
        
!         call nc_write(path_out,"z_srf",z_srf,dim1="xc",dim2="yc",missing_value=mv_wrt, &
!                       units="m",long_name="Surface elevation")
!         call nc_write(path_out,"H",H,dim1="xc",dim2="yc",missing_value=mv_wrt, &
!                       units="m",long_name="Ice thickness")

        call nc_write(path_out,"ice_age",trc%stats%ice_age,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age")
        call nc_write(path_out,"ice_age_err",trc%stats%ice_age_err,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age - error")
        call nc_write(path_out,"density",trc%stats%density,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density")

        call nc_write(path_out,"depth_iso",trc%stats%depth_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone depth")
        call nc_write(path_out,"depth_iso_err",trc%stats%depth_iso_err,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone depth - error")
        call nc_write(path_out,"dep_z_iso",trc%stats%dep_z_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone deposition elevation")
        call nc_write(path_out,"density_iso",trc%stats%density_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density (for isochrones)")

        
        return 

    end subroutine tracer_write_stats
    
    subroutine tracer_read(trc,filename,time)

        implicit none 

        type(tracer_class), intent(OUT) :: trc
        character(len=*),   intent(IN)  :: filename 
        real(prec_time),    intent(IN)  :: time 

        ! Local variables 
        integer :: i, k 



        return 

    end subroutine tracer_read

    subroutine tracer_align(trc_new,trc_ref,trc,dxy_max,dz_max)
        ! Interpolate tracer ages from trc to those of trc_ref 
        ! Note: for this to work well, trc should be sufficiently high
        ! resolution to minimize interpolation errors 

        implicit none 

        type(tracer_class), intent(OUT) :: trc_new 
        type(tracer_class), intent(IN) :: trc_ref, trc  
        real(prec), intent(IN) :: dxy_max, dz_max  

        ! Local variables 
        integer :: i, k  
        real(prec) :: dist_xy(trc%par%n), dist_z(trc%par%n)

        ! Store reference tracer information in new object 
        trc_new = trc_ref 

        ! Make sure to set tagged info to missing, since
        ! this will not be valid for trc_new 
        trc_new%dep%time = MV 
        trc_new%dep%H    = MV 
        trc_new%dep%x    = MV 
        trc_new%dep%y    = MV 
        trc_new%dep%z    = MV 
        
        do i = 1, trc_new%par%n 

            if (trc_new%now%active(i) .eq. 2) then 
                ! Only treat active locations 

                dist_xy = MV 
                dist_z  = MV
                where (trc%now%active .eq. 2)
                    dist_xy = sqrt( (trc_new%now%x(i)-trc%now%x)**2 &
                                  + (trc_new%now%y(i)-trc%now%y)**2)
                    dist_z  = abs(trc_new%now%z(i)-trc%now%z)
                end where 

                k = minloc(dist_xy,mask=dist_xy.ne.MV.and.dist_z.le.dz_max,dim=1)

                if (dist_xy(k) .le. dxy_max) then 
                    trc_new%dep%time(i) = trc%dep%time(i)
                else 
                    trc_new%now%active(i) = 0
                    trc_new%now%H(i)      = MV 
                    trc_new%now%z_srf(i)  = MV
                    trc_new%now%x(i)      = MV 
                    trc_new%now%y(i)      = MV 
                    trc_new%now%z(i)      = MV 
         
                end if 

            end if 

        end do 

        return 

    end subroutine tracer_align

    subroutine tracer_import_eulerian(trc,time,x,y,z,age,z_srf,H,is_sigma,order,sigma_srf)
        ! Given a 3D field of Eulerian ages on a grid, 
        ! convert to tracer format. 

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:) 
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
        real(prec), intent(IN) :: age(:,:,:) 
        logical,    intent(IN) :: is_sigma 
        character(len=*), intent(IN), optional :: order 
        real(prec), intent(IN), optional :: sigma_srf     ! Value at surface by default (1 or 0?)

        ! Local variables  
        character(len=3) :: idx_order 
        integer :: nx, ny, nz 
        real(prec), allocatable :: x1(:), y1(:), z1(:)
        real(prec), allocatable :: z_srf1(:,:), H1(:,:)
        real(prec), allocatable :: age1(:,:,:) 
        real(prec) :: zc(size(z))
        logical :: rev_z 

        ! Determine order of indices (default ijk)
        idx_order = "ijk"
        if (present(order)) idx_order = trim(order)

        ! Correct the sigma values if necessary,
        ! so that sigma==0 [base]; sigma==1 [surface]
        zc = z 
        if (trc%par%is_sigma .and. present(sigma_srf)) then 
            if (sigma_srf .eq. 0.0) then 
                ! Adjust sigma values 
                zc = 1.0 - z 
            end if 
        end if 

        ! Also determine whether z-axis is initially ascending or descending 
        rev_z = (zc(1) .gt. zc(size(zc)))

        call tracer_reshape1D_vec(x, x1,rev=.FALSE.)
        call tracer_reshape1D_vec(y, y1,rev=.FALSE.)
        call tracer_reshape1D_vec(real(zc,kind=prec),z1,rev=rev_z)
        call tracer_reshape2D_field(idx_order,z_srf,z_srf1)
        call tracer_reshape2D_field(idx_order,H,H1)
        call tracer_reshape3D_field(idx_order,age,age1,rev_z=rev_z)
        
        ! Get axis sizes (if ny==2, this is a 2D profile domain) 
        nx = size(x1,1)
        ny = size(y1,1)
        nz = size(z1,1)





        return 

    end subroutine tracer_import_eulerian


    ! ======================================================================
    !
    ! 2D writing interface 
    !
    ! ======================================================================

    subroutine tracer2D_write_init(trc,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 

        call tracer_write_init(trc,fldr,filename)

        return 

    end subroutine tracer2D_write_init 

    subroutine tracer2D_write(trc,time,fldr,filename)
        ! Wrapper to calling normal tracer_write routine
        
        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time) :: time 
        character(len=*), intent(IN) :: fldr, filename 

        call tracer_write(trc,time,fldr,filename,is2D=.TRUE.)

        return 

    end subroutine tracer2D_write 

    subroutine tracer2D_write_stats(trc,time,fldr,filename) !,z_srf,H)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec_time) :: time
        character(len=*), intent(IN)   :: fldr, filename 
!         real(prec),         intent(IN) :: z_srf(:), H(:) 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 
        real(prec_wrt) :: mv_wrt 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",        x=trc%stats%x*1e-3,     units="km")
        call nc_write_dim(path_out,"depth_norm",x=trc%stats%depth_norm, units="1")
        call nc_write_dim(path_out,"age_iso",   x=trc%stats%age_iso,    units="ka")
        call nc_write_dim(path_out,"time",      x=time,unlimited=.TRUE.,units="ka")
        
!         call nc_write(path_out,"z_srf",z_srf,dim1="xc",missing_value=mv_wrt, &
!                       units="m",long_name="Surface elevation")
!         call nc_write(path_out,"H",H,dim1="xc",missing_value=mv_wrt, &
!                       units="m",long_name="Ice thickness")

        call nc_write(path_out,"ice_age",trc%stats%ice_age(:,1,:),dim1="xc",dim2="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age")
        call nc_write(path_out,"ice_age_err",trc%stats%ice_age_err(:,1,:),dim1="xc",dim2="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age - error")

        call nc_write(path_out,"depth_iso",trc%stats%depth_iso(:,1,:),dim1="xc",dim2="age_iso",missing_value=mv_wrt, &
                      units="ka",long_name="Isochrone depth")
        call nc_write(path_out,"depth_iso_err",trc%stats%depth_iso_err(:,1,:),dim1="xc",dim2="age_iso",missing_value=mv_wrt, &
                      units="ka",long_name="Isochrone depth - error")
        
        call nc_write(path_out,"density",trc%stats%density(:,1,:),dim1="xc",dim2="depth_norm",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density")
        
        return 

    end subroutine tracer2D_write_stats

end module tracer_io
