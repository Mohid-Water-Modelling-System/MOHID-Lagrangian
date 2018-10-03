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
    use field_types_mod

    implicit none
    private

    type :: interpolator_class        !< Interpolator class
        integer :: interpType = 1     !< Interpolation Algorithm 1: linear
        type(string) :: name                !< Name of the Interpolation algorithm
    contains
    procedure :: run
    procedure :: initialize => initInterpolator
    procedure :: print => printInterpolator
    procedure :: interp4D
    end type interpolator_class

    !Public access vars
    public :: interpolator_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that runs the chosen interpolator method on the given data.
    !> @param[in] self, aot, bdata, time, dt, var_dt, var_name
    !---------------------------------------------------------------------------
    subroutine run(self, aot, bdata, time, var_dt, var_name)
    class(interpolator_class), intent(in) :: self
    type(aot_class), intent(in) :: aot
    type(background_class), intent(in) :: bdata
    real(prec_time), intent(in) :: time
    real(prec), dimension(:,:), intent(out) :: var_dt
    type(string), dimension(:), intent(out) :: var_name
    real(prec_time) :: newtime
    class(*), pointer :: aField
    integer :: i = 1
    type(string) :: outext

    !Check field extents and what particles will be interpolated
    !interpolate each field to the correspoing slice in var_dt
    call bdata%fields%reset()                   ! reset list iterator
    do while(bdata%fields%moreValues())         ! loop while there are values
        aField => bdata%fields%currentValue()   ! get current value
        select type(aField)
        class is(scalar4d_field_class)          !4D interpolation is possible
            if (self%interpType == 1) then !linear interpolation in space and time
                var_name(i) = aField%name
                call self%interp4D(aot%x, aot%y, aot%z, time, aField%field, var_dt(:,i), size(aField%field,1), size(aField%field,2), size(aField%field,3), size(aField%field,4), size(aot%x))
            end if !add more interpolation types here
        class is(scalar3d_field_class)          !3D interpolation is possible
            if (self%interpType == 1) then !linear interpolation in space and time
                !call self%interp3D(...)
            end if !add more interpolation types here
            !add more field types here
            class default
            outext = '[Interpolator::Run] Unexepected type of field, not correct or supported at the time'
            call Log%put(outext)
            stop
        end select
        call bdata%fields%next()                ! increment the list iterator
        i = i+1 !to select the correct slice of var_dt for the corresponding field
    end do
    call bdata%fields%reset()                   ! reset list iterator

    end subroutine run

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> method to interpolate a particle position in a given data box based
    !> on array coordinates. 4d interpolation is a weighted average of 16
    !> neighbors. Consider the 4D domain between the 16 neighbors. The hypercube is
    !> divided into 16 sub-hypercubes by the point in question. The weight of each
    !> neighbor is given by the volume of the opposite sub-hypercube, as a fraction
    !> of the whole hypercube.
    !> @param[in] self, x, y, z, t, field, f_out, n_fv, n_cv, n_pv, n_tv, n_e
    !---------------------------------------------------------------------------
    subroutine interp4D(self, x, y, z, t, field, f_out, n_fv, n_cv, n_pv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y, z                       !< 1-d. Array of particle component positions
    real(prec_time), intent(in) :: t                                           !< time to interpolate to
    real(prec), dimension(n_fv, n_cv, n_pv, n_tv), intent(in) :: field    !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    real(prec), dimension(n_e), intent(out) :: f_out                      !< Field evaluated at x,y,z,t
    integer, intent(in) :: n_fv, n_cv, n_pv, n_tv                   !< field dimensions
    integer, intent(in) :: n_e                                      !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, z0, x1, y1, z1
    real(prec), dimension(n_e) :: xd, yd, zd, c000, c100, c010, c110, c001
    real(prec), dimension(n_e) :: c101, c011, c111, c00, c10, c01, c11, c0, c1
    real(prec) :: td
    integer :: i, j, k, l, t0, t1

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
        c000(i) = field(x0(i),y0(i),z0(i),t0)*(1.-xd(i)) + field(x1(i),y0(i),z0(i),t0)*xd(i) !y0x0z0t0!  y0x1z0t0
        c100(i) = field(x0(i),y1(i),z0(i),t0)*(1.-xd(i)) + field(x1(i),y1(i),z0(i),t0)*xd(i)
        c010(i) = field(x0(i),y0(i),z1(i),t0)*(1.-xd(i)) + field(x1(i),y0(i),z1(i),t0)*xd(i)
        c110(i) = field(x0(i),y1(i),z1(i),t0)*(1.-xd(i)) + field(x1(i),y1(i),z1(i),t0)*xd(i)

        c001(i) = field(x0(i),y0(i),z0(i),t1)*(1.-xd(i)) + field(x1(i),y0(i),z0(i),t1)*xd(i) !y0x0z0t0!  y0x1z0t0
        c101(i) = field(x0(i),y1(i),z0(i),t1)*(1.-xd(i)) + field(x1(i),y1(i),z0(i),t1)*xd(i)
        c011(i) = field(x0(i),y0(i),z1(i),t1)*(1.-xd(i)) + field(x1(i),y0(i),z1(i),t1)*xd(i)
        c111(i) = field(x0(i),y1(i),z1(i),t1)*(1.-xd(i)) + field(x1(i),y1(i),z1(i),t1)*xd(i)
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
    self%interpType = flag
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