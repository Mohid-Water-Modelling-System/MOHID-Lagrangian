
module tracer_interp 

    use tracer_precision 
    use bspline_module 

    implicit none 

    type lin_interp_par_type 
        integer    :: i1, i2 
        real(prec) :: alpha
    end type 

    type lin3_interp_par_type 
        integer    :: i1, i2, j1, j2, k1, k2  
        real(prec) :: alpha_x, alpha_y, alpha_z
    end type 

    private

    public :: lin_interp_par_type
    public :: lin3_interp_par_type
    public :: interp_bilinear_weights, interp_bilinear 
    public :: interp_trilinear_weights, interp_trilinear  

    public :: interp_bspline3D_weights, interp_bspline3D
    
contains 

    subroutine calc_interp_linear_weights(idx1,idx2,alpha,x,xout)

        implicit none 

        integer,    intent(OUT) :: idx1, idx2  
        real(prec), intent(OUT) :: alpha 
        real(prec), intent(IN)  :: x(:)
        real(prec), intent(IN)  :: xout

        if (x(1) .lt. x(size(x))) then 

            call calc_interp_linear_weights_ascending(idx1,idx2,alpha,x,xout)

        else 
            
            call calc_interp_linear_weights_descending(idx1,idx2,alpha,x,xout)

        end if 

        return 

    end subroutine calc_interp_linear_weights 

    subroutine calc_interp_linear_weights_ascending(idx1,idx2,alpha,x,xout)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        integer,    intent(OUT) :: idx1, idx2  
        real(prec), intent(OUT) :: alpha 
        real(prec), intent(IN)  :: x(:)
        real(prec), intent(IN)  :: xout

        integer :: i, nx

        nx = size(x,1)

        ! By default, set right index out of range and weight to zero
        idx2   = -1 
        alpha = 0.0 

        ! Get x-index corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest (inside range)
        if (xout .gt. x(1) .and. xout .le. x(nx)) then 

            do i = 1, nx 
                if (x(i) .ge. xout) exit 
            end do 

            idx2   = i 
            alpha = (xout - x(i-1)) / (x(i)-x(i-1))

        else if (xout .eq. x(1)) then 

            idx2   = 2 
            alpha = 0.0 

        end if

        ! Define left index too 
        idx1 = idx2 - 1 

        return 

    end subroutine calc_interp_linear_weights_ascending

    subroutine calc_interp_linear_weights_descending(idx1,idx2,alpha,x,xout)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        ! Note: for axis that is decreasing in value (descending)

        implicit none 

        integer,    intent(OUT) :: idx1, idx2
        real(prec), intent(OUT) :: alpha 
        real(prec), intent(IN)  :: x(:)
        real(prec), intent(IN)  :: xout

        integer :: i, nx

        nx = size(x,1)

        ! By default, set right index out of range and weight to zero
        idx2   = -1 
        alpha = 0.0 

        ! Get x-index corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest (inside range)
        if (xout .gt. x(nx) .and. xout .le. x(1)) then 

            do i = nx, 1, -1 
                if (x(i) .ge. xout) exit 
            end do 

            idx2   = i 
            alpha = (xout - x(i+1)) / (x(i)-x(i+1))

        else if (xout .eq. x(nx)) then 

            idx2   = nx-1 
            alpha = 0.0 

        end if

        ! Set left index too 
        idx1 = idx2 + 1 
        return 

    end subroutine calc_interp_linear_weights_descending

    function interp_bilinear_weights(x,y,xout,yout) result(par)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        real(prec), intent(IN) :: x(:), y(:) 
        real(prec), intent(IN) :: xout, yout
        type(lin3_interp_par_type)   :: par 

        call calc_interp_linear_weights(par%i1,par%i2,par%alpha_x,x,xout)
        call calc_interp_linear_weights(par%j1,par%j2,par%alpha_y,y,yout)
        
        return 

    end function interp_bilinear_weights

    function interp_trilinear_weights(x,y,z,xout,yout,zout) result(par)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        real(prec), intent(IN) :: x(:), y(:), z(:) 
        real(prec), intent(IN) :: xout, yout, zout
        type(lin3_interp_par_type)   :: par 

        call calc_interp_linear_weights(par%i1,par%i2,par%alpha_x,x,xout)
        call calc_interp_linear_weights(par%j1,par%j2,par%alpha_y,y,yout)
        call calc_interp_linear_weights(par%k1,par%k2,par%alpha_z,z,zout)
        
        return 

    end function interp_trilinear_weights

    function interp_bilinear(par,var) result (varout)

        implicit none 

        type(lin3_interp_par_type), intent(IN) :: par
        real(prec), intent(IN) :: var(:,:) 
        real(prec) :: varout 
        real(prec) :: p1, p2
        integer    :: i1, i2, j1, j2

        varout = MV 

        if (par%i1 .gt. 0 .and. par%j1 .gt. 0 .and. &
            par%i2 .gt. 0 .and. par%j2 .gt. 0) then

            i1 = par%i1
            i2 = par%i2  
            j1 = par%j1 
            j2 = par%j2 

            ! Lower z-plane 
            p1 = var(i1,j1) + par%alpha_x*(var(i2,j1)-var(i1,j1))
            p2 = var(i1,j2) + par%alpha_x*(var(i2,j2)-var(i1,j2))
            varout = p1 + par%alpha_y*(p2-p1)

        end if 

        return 

    end function interp_bilinear

    function interp_trilinear(par,var) result (varout)

        implicit none 

        type(lin3_interp_par_type), intent(IN) :: par
        real(prec), intent(IN) :: var(:,:,:) 
        real(prec) :: varout 
        real(prec) :: p1, p2, v1, v2 
        integer    :: i1, i2, j1, j2, k1, k2  

        varout = MV 

        if (par%i1 .gt. 0 .and. par%j1 .gt. 0 .and. par%k1 .gt. 0 .and. &
            par%i2 .gt. 0 .and. par%j2 .gt. 0 .and. par%k2 .gt. 0) then

            i1 = par%i1
            i2 = par%i2
            j1 = par%j1
            j2 = par%j2
            k1 = par%k1
            k2 = par%k2

            ! Lower z-plane 
            p1 = var(i1,j1,k1) + par%alpha_x*(var(i2,j1,k1)-var(i1,j1,k1))
            p2 = var(i1,j2,k1) + par%alpha_x*(var(i2,j2,k1)-var(i1,j2,k1))
            v1 = p1 + par%alpha_y*(p2-p1)

            ! Upper z-plane 
            p1 = var(i1,j1,k2) + par%alpha_x*(var(i2,j1,k2)-var(i1,j1,k2))
            p2 = var(i1,j2,k2)   + par%alpha_x*(var(i2,j2,k2)-var(i1,j2,k2))
            v2 = p1 + par%alpha_y*(p2-p1)

            ! Linear z-interpolation 
            varout = v1 + par%alpha_z*(v2-v1)

        end if 

        return 

    end function interp_trilinear


    ! === bspline interface ===


    subroutine interp_bspline3D_weights(bspl3D,x,y,z,var)

        implicit none 

        type(bspline_3d), intent(INOUT) :: bspl3D 
        real(prec),       intent(IN)    :: x(:), y(:), z(:)
        real(prec),       intent(IN)    :: var(:,:,:)

        ! Local variables
        integer :: bspline_flag 

        call bspl3D%initialize(dble(x),dble(y),dble(z),dble(var),kx=4,ky=4,kz=4,iflag=bspline_flag)
        
        if (bspline_flag .ne. 0) then 
            write(*,*) "interp_bspline3D_weights:: error initializing."
            stop 
        end if     

        return 

    end subroutine interp_bspline3D_weights 

    function interp_bspline3D(bspl3D,x,y,z) result(var)

        implicit none 

        type(bspline_3d), intent(INOUT) :: bspl3D 
        real(prec),       intent(IN)    :: x, y, z
        real(prec)                      :: var

        ! Local variables
        integer :: idx, idy, idz 
        integer :: bspline_flag 
        double precision :: var_dbl 

        idx = 0 
        idy = 0 
        idz = 0 

        call bspl3D%evaluate(dble(x),dble(y),dble(z),idx,idy,idz,var_dbl,bspline_flag)

        if (bspline_flag .ne. 0) then 
            write(*,*) "interp_bspline3D:: error interpolating."
            stop 
        end if     

        var = var_dbl 
        
        return 

    end function interp_bspline3D


end module tracer_interp 
