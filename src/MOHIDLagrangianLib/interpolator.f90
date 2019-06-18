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
    use fieldTypes_mod

    implicit none
    private

    type :: interpolator_class        !< Interpolator class
        integer :: interpType = 1     !< Interpolation Algorithm 1: linear
        type(string) :: name                !< Name of the Interpolation algorithm
    contains
    procedure :: run
    procedure :: getArrayCoordRegular
    procedure :: getPointCoordRegular
    procedure :: initialize => initInterpolator
    procedure :: print => printInterpolator
    procedure :: test4D
    procedure, private :: interp4D
    end type interpolator_class

    !Public access vars
    public :: interpolator_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that runs the chosen interpolator method on the given data.
    !> @param[in] self, aot, bdata, time, var_dt, var_name
    !---------------------------------------------------------------------------
    subroutine run(self, aot, bdata, time, var_dt, var_name)
    class(interpolator_class), intent(in) :: self
    type(aot_class), intent(in) :: aot
    type(background_class), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(:,:), intent(out) :: var_dt
    type(string), dimension(:), intent(out) :: var_name
    real(prec) :: newtime
    class(*), pointer :: aField
    integer :: i
    type(string) :: outext

    real(prec), dimension(size(aot%x)) :: xx, yy, zz
    real(prec) :: tt

    !Check field extents and what particles will be interpolated
    !interpolate each field to the correspoing slice in var_dt
    i = 1
    call bdata%fields%reset()                   ! reset list iterator
    do while(bdata%fields%moreValues())         ! loop while there are values
        aField => bdata%fields%currentValue()   ! get current value
        select type(aField)        
        class is(scalar4d_field_class)          !4D interpolation is possible
            if (self%interpType == 1) then !linear interpolation in space and time
                var_name(i) = aField%name
                xx = self%getArrayCoordRegular(aot%x, bdata, Globals%Var%lon)
                yy = self%getArrayCoordRegular(aot%y, bdata, Globals%Var%lat)
                zz = self%getArrayCoordRegular(aot%z, bdata, Globals%Var%level)
                tt = self%getPointCoordRegular(time, bdata, Globals%Var%time, -Globals%SimDefs%dt)
                var_dt(:,i) = self%interp4D(xx, yy, zz, tt, aField%field, size(aField%field,1), size(aField%field,2), size(aField%field,3), size(aField%field,4), size(aot%x))
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
    !> @param[in] self, x, y, z, t, field, n_fv, n_cv, n_pv, n_tv, n_e
    !---------------------------------------------------------------------------
    function interp4D(self, x, y, z, t, field, n_fv, n_cv, n_pv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y, z                       !< 1-d. Array of particle component positions in array coordinates
    real(prec), intent(in) :: t                                      !< time to interpolate to in array coordinates
    real(prec), dimension(n_fv, n_cv, n_pv, n_tv), intent(in) :: field    !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    integer, intent(in) :: n_fv, n_cv, n_pv, n_tv                         !< field dimensions
    integer, intent(in) :: n_e                                            !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, z0, x1, y1, z1
    real(prec), dimension(n_e) :: xd, yd, zd, c000, c100, c010, c110, c001
    real(prec), dimension(n_e) :: c101, c011, c111, c00, c10, c01, c11, c0, c1
    real(prec) :: td
    integer :: i, j, k, l, t0, t1
    real(prec), dimension(n_e) :: interp4D                      !< Field evaluated at x,y,z,t
    
    ! From x,y,z,t in array coordinates, find the the box inside the field where the particle is
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
    interp4D = c0*(1.-td)+c1*td

    end function interp4D
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinates of a set of points, given a coordinate 
    !> array. Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dimName
    !---------------------------------------------------------------------------
    function getArrayCoordRegular(self, xdata, bdata, dimName)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                 !< Background to use
    type(string), intent(in) :: dimName
    integer :: dim                                              !< corresponding background dimension
    real(prec), dimension(size(xdata)) :: getArrayCoordRegular  !< coordinates in array index
    real(prec) :: minBound, maxBound, res
    type(string) :: outext
    dim = bdata%getDimIndex(dimName)
    res = size(bdata%dim(dim)%field)-1
    minBound = bdata%dim(dim)%getFieldMinBound()
    maxBound = bdata%dim(dim)%getFieldMaxBound()
    res = abs(maxBound - minBound)/res
    getArrayCoordRegular = (xdata - minBound)/res + 1
    !if (.not.Utils%isBounded(xdata, minBound, maxBound, res/3.0)) then
    !    outext = '[Interpolator::getArrayCoordRegular] Points not contained in "'//dimName//'" dimension, stoping'
    !    call Log%put(outext)
    !    stop
    !end if
    end function getArrayCoordRegular


    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension. 
    !> @param[in] self, xdata, bdata, dimName
    ! !---------------------------------------------------------------------------
    ! function getPointCoordNonRegular(self, xdata, bdata, dimName, eta)
    ! class(interpolator_class), intent(in) :: self
    ! real(prec), intent(in):: xdata              !< Tracer coordinate component
    ! type(background_class), intent(in) :: bdata !< Background to use
    ! type(string), intent(in) :: dimName
    ! real(prec), intent(in), optional :: eta
    ! integer :: i,idx_1,idx_2                              !< corresponding background dimension
    ! real(prec) :: getPointCoordNonRegular          !< coordinates in array index
    ! real(prec) :: minBound, maxBound, res, ieta
    ! type(string) :: outext

    ! do i=1,size(xdata)
    !     idx_1 = minloc(xdata(i)-bdata%dim(dim)%field)
    !     idx_2 = idx_1+1
    !     getPointCoordNonRegular(i) = idx_1 + xdata(i)*(idx_2-idx_1)/(bdata%dim(dim)%values(idx_2)-bdata%dim(dim)%field(idx_1))
    ! end do
    ! ieta = -res/10.0
    ! if (present(eta)) ieta = eta
    ! if (.not.Utils%isBounded(xdata, minBound, maxBound, ieta)) then
    !     outext = '[Interpolator::getPointCoordNonRegular] Point not contained in "'//dimName//'" dimension, stoping'
    !     print*, xdata
    !     print*, minBound, maxBound+ieta
    !     call Log%put(outext)
    !     stop
    ! end if
    ! end function getPointCoordNonRegular

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension. 
    !> Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dimName
    !---------------------------------------------------------------------------
    function getPointCoordRegular(self, xdata, bdata, dimName, eta)
        class(interpolator_class), intent(in) :: self
        real(prec), intent(in):: xdata              !< Tracer coordinate component
        type(background_class), intent(in) :: bdata !< Background to use
        type(string), intent(in) :: dimName
        real(prec), intent(in), optional :: eta
        integer :: dim                              !< corresponding background dimension
        real(prec) :: getPointCoordRegular          !< coordinates in array index
        real(prec) :: minBound, maxBound, res, ieta
        type(string) :: outext
        dim = bdata%getDimIndex(dimName) 
        res = size(bdata%dim(dim)%field)-1
        minBound = bdata%dim(dim)%getFieldMinBound()
        maxBound = bdata%dim(dim)%getFieldMaxBound()
        res = abs(maxBound - minBound)/res
        getPointCoordRegular = (xdata - minBound)/res+1
        ieta = -res/10.0
        if (present(eta)) ieta = eta
        if (.not.Utils%isBounded(xdata, minBound, maxBound, ieta)) then
            outext = '[Interpolator::getPointCoordRegular] Point not contained in "'//dimName//'" dimension, stoping'
            print*, xdata
            print*, minBound, maxBound+ieta
            call Log%put(outext)
            stop
        end if
        end function getPointCoordRegular

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
    !> Test for interp 4D
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine test4D(self)
    class(interpolator_class), intent(inout) :: self
    real(prec), dimension(:,:,:,:), allocatable :: field
    real(prec), dimension(:), allocatable :: xx, yy, zz 
    real(prec) :: time
    integer :: npts, fieldDims
    real(prec) :: fieldVal
    fieldDims = 10
    fieldVal = 1.9
    allocate(field(fieldDims,fieldDims,fieldDims,fieldDims))
    field = fieldVal
    npts = 1 
    allocate(xx(npts), yy(npts), zz(npts))
    xx = 13.45
    yy = xx
    zz = xx
    time = 1
    print*, 'testing 4D interpolation, expected result is ', fieldVal
    xx = self%interp4D(xx, yy, zz, time, field, fieldDims, fieldDims, fieldDims, fieldDims, npts)
    print*, 'result = ', xx
    if (xx(1) == fieldVal) then 
        print*, 'Test: SUCCESS'
    else
        print*, 'Test: FAILED'
    end if
    read(*,*)
    end subroutine test4D

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