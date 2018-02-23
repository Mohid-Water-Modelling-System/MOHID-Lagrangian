
module tracer2D 
    ! Module to wrap a 2D (profile) version of the tracer model
    ! Makes calls to main tracer code by reshaping profile into
    ! 3D array with y-dimension thickness of 1. 

    use tracer_precision
    use tracer3D
    use ncio    

    implicit none 

    private 
    public :: tracer2D_init 
    public :: tracer2D_update 
    public :: tracer2D_end 

contains


    subroutine tracer2D_init(trc,filename,time,x,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:)
        logical,    intent(IN) :: is_sigma 
        real(prec_time) :: time 

        real(prec) :: y(5) 

        ! Define the ghost y-dimension
        y(1:5) = [0.0,0.25,0.50,0.75,1.0] 

        ! Call 3D tracer_init
        call tracer_init(trc,filename,time,x,y,is_sigma)

        return 

    end subroutine tracer2D_init


    subroutine tracer2D_update(trc,time,x,z,z_srf,H,ux,uz, &
                               lon,lat,t2m_ann,t2m_sum,pr_ann,pr_sum,d18O_ann, &
                               dep_now,stats_now)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), z(:)
        real(prec), intent(IN) :: z_srf(:), H(:)
        real(prec), intent(IN) :: ux(:,:), uz(:,:)
        real(prec), intent(IN) :: lon(:), lat(:), t2m_ann(:), t2m_sum(:), pr_ann(:), pr_sum(:), d18O_ann(:)
        logical,    intent(IN) :: dep_now, stats_now 

        ! Local variables
        real(prec) :: y(2) 
        real(prec), allocatable :: z_srf_2D(:,:), H_2D(:,:)
        real(prec), allocatable :: ux_3D(:,:,:), uy_3D(:,:,:), uz_3D(:,:,:)
        real(prec), allocatable :: lon_2D(:,:), lat_2D(:,:), t2m_ann_2D(:,:), t2m_sum_2D(:,:) 
        real(prec), allocatable :: pr_ann_2D(:,:), pr_sum_2D(:,:), d18O_ann_2D(:,:)
        
        integer :: j, ny 

        ny = size(y,1)

        ! Define ghost dimension and data 
        allocate(z_srf_2D(size(x,1),ny))
        allocate(H_2D(size(x,1),ny))
        allocate(ux_3D(size(ux,1),ny,size(ux,2)))
        allocate(uy_3D(size(ux,1),ny,size(ux,2)))
        allocate(uz_3D(size(ux,1),ny,size(ux,2)))
        allocate(lon_2D(size(x,1),ny))
        allocate(lat_2D(size(x,1),ny))
        allocate(t2m_ann_2D(size(x,1),ny))
        allocate(t2m_sum_2D(size(x,1),ny))
        allocate(pr_ann_2D(size(x,1),ny))
        allocate(pr_sum_2D(size(x,1),ny))
        allocate(d18O_ann_2D(size(x,1),ny))
        
        ! Set y-dimension to one value of zero
        y(1:2) = [0.0,1.0] 

        ! Reshape input data with a ghost y-dimension of length two
        do j = 1, size(y)

            z_srf_2D(:,j)    = z_srf 
            H_2D(:,j)        = H 

            ux_3D(:,j,:)     = ux 
            uy_3D            = 0.0 
            uz_3D(:,j,:)     = uz 

            lon_2D(:,j)      = lon 
            lat_2D(:,j)      = lat 
            t2m_ann_2D(:,j)  = t2m_ann 
            t2m_sum_2D(:,j)  = t2m_sum 
            pr_ann_2D(:,j)   = pr_ann 
            pr_sum_2D(:,j)   = pr_sum 
            d18O_ann_2D(:,j) = d18O_ann 
            
        end do 

        ! Now update tracers using 3D call 
        call tracer_update(trc,time,x,y,z,z_srf_2D,H_2D,ux_3D,uy_3D,uz_3D, &
                            lon_2D,lat_2D,t2m_ann_2D,t2m_sum_2D,pr_ann_2D,pr_sum_2D,d18O_ann_2D, &
                            dep_now,stats_now,order="ijk")

        return 

    end subroutine tracer2D_update

    subroutine tracer2D_end(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        
        ! Call normal tracer_end subroutine 
        call tracer_end(trc) 

        return 

    end subroutine tracer2D_end


end module tracer2D

