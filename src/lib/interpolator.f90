    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : interpolator
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : September 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines an Interpolator class.
    !------------------------------------------------------------------------------

    module interpolator_mod

    use common_modules
    use AoT_mod
    use background_mod

    implicit none
    private

    type :: interpolator_class        !< Interpolator class
        integer :: interpolatorType = 1     !< Interpolation Algorithm
        type(string) :: name                !< Name of the Interpolation algorithm
    contains
    procedure :: initialize => initInterpolator
    procedure :: print => printInterpolator
    end type interpolator_class

    !Public access vars
    public :: interpolator_class

    contains
    
    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> method to interpolate a particle position in a given data box based
    !> on array coordinates. 4d interpolation is a weighted average of 16
    !> neighbors. Consider the 4D domain between the 16 neighbors. The hypercube is
    !> divided into 16 sub-hypercubes by the point in question. The weight of each
    !> neighbor is given by the volume of the opposite sub-hypercube, as a fraction
    !> of the whole hypercube. Explanation adapted from:
    !>
    !> @param[in] x,y,z  1-d. Array of particles, contains the component positions.
    !> @param[in] t      Time
    !> @param[in] f_in   Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    !> @param[in] n_e    Number of particles
    !> @param[out] f_out  Field evaluated at x,y,z,t
    !---------------------------------------------------------------------------
    subroutine interp4D(self, x, y, z, t, f_in, f_out, n_fv, n_cv, n_pv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    integer, intent(in) :: n_fv, n_cv, n_pv, n_tv, n_e
    real, dimension(n_e),intent(in):: x, y, z
    real, intent(in) :: t
    real, dimension(n_fv, n_cv, n_pv, n_tv), intent(in) :: f_in
    real, dimension(n_e), intent(out) :: f_out
    integer, dimension(n_e) :: x0, y0, z0, x1, y1, z1
    integer :: t0, t1
    real, dimension(n_e) :: xd, yd, zd, c000, c100, c010, c110, c001, c101, c011, c111, c00, c10, c01, c11, c0, c1
    real :: td
    integer :: i, j, k, l    

    ! From x,y,z,t in array coordinates, find the the box inside the field where the partcle is
    x0 = floor(x)
    y0 = floor(y)
    z0 = floor(z)
    t0 = floor(t)
    x1 = ceiling(x)
    y1 = ceiling(y)
    z1 = ceiling(z)
    t1 = ceiling(t)

    ! Compute the "normalized coordinates" of the particle inside the data field box
    xd = (x-x0)/(x1-x0)
    yd = (y-y0)/(y1-y0)
    zd = (z-z0)/(z1-z0)
    td = (t-t0)/(t1-t0)
    
    ! In case that particle is on a point box, we set it to 0 to avoid inf errors
    where (x1 == x0) xd = 0.
    where (y1 == y0) yd = 0.
    where (z1 == z0) zd = 0.
    if (t1 == t0)    td = 0.
    
    ! Interpolation on the first dimension and collapse it to a three dimension problem
    forall(i=1:n_e)
        c000(i) = f_in(x0(i),y0(i),z0(i),t0)*(1.-xd(i)) + f_in(x1(i),y0(i),z0(i),t0)*xd(i) !y0x0z0t0!  y0x1z0t0
        c100(i) = f_in(x0(i),y1(i),z0(i),t0)*(1.-xd(i)) + f_in(x1(i),y1(i),z0(i),t0)*xd(i)
        c010(i) = f_in(x0(i),y0(i),z1(i),t0)*(1.-xd(i)) + f_in(x1(i),y0(i),z1(i),t0)*xd(i)
        c110(i) = f_in(x0(i),y1(i),z1(i),t0)*(1.-xd(i)) + f_in(x1(i),y1(i),z1(i),t0)*xd(i)

        c001(i) = f_in(x0(i),y0(i),z0(i),t1)*(1.-xd(i)) + f_in(x1(i),y0(i),z0(i),t1)*xd(i) !y0x0z0t0!  y0x1z0t0
        c101(i) = f_in(x0(i),y1(i),z0(i),t1)*(1.-xd(i)) + f_in(x1(i),y1(i),z0(i),t1)*xd(i)
        c011(i) = f_in(x0(i),y0(i),z1(i),t1)*(1.-xd(i)) + f_in(x1(i),y0(i),z1(i),t1)*xd(i)
        c111(i) = f_in(x0(i),y1(i),z1(i),t1)*(1.-xd(i)) + f_in(x1(i),y1(i),z1(i),t1)*xd(i)
    end forall

    ! Interpolation on the second dimension and collapse it to a two dimension problem
    c00 = c000*(1.-yd)+c100*yd
    c10 = c010*(1.-yd)+c110*yd
    c01 = c001*(1.-yd)+c101*yd
    c11 = c011*(1.-yd)+c111*yd

    ! Interpolation on the third dimension and collapse it to a one dimension problem
    c0 = c00*(1.-zd)+c10*zd
    c1 = c01*(1.-zd)+c11*zd

    ! Interpolation on the time dimension and get the final result.
    f_out = c0*(1.-td)+c1*td

    end subroutine interp4D


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializer method for the Interpolator class. Sets the type of interpolator
    !> and name of the algorithm this Interpolator will call
    !> @param[in] self, flag, name
    !---------------------------------------------------------------------------
    subroutine initInterpolator(self, flag, name)
    class(interpolator_class), intent(inout) :: self
    integer, intent(in) :: flag
    type(string), intent(in) :: name
    self%interpolatorType = flag
    self%name = name
    end subroutine initInterpolator

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the Interpolator information
    !---------------------------------------------------------------------------
    subroutine printInterpolator(self)
    class(interpolator_class), intent(inout) :: self
    type(string) :: outext, t
    outext = 'Interpolation algorithm is '//self%name
    call Log%put(outext,.false.)
    end subroutine printInterpolator

    end module interpolator_mod