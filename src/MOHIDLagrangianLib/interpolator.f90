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
    use stateVector_mod
    use background_mod
    use fieldTypes_mod
    USE ieee_arithmetic

    implicit none
    private

    type :: interpolator_class        !< Interpolator class
        integer :: interpType = 1     !< Interpolation Algorithm 1: linear
        type(string) :: name                !< Name of the Interpolation algorithm
    contains
    procedure :: run
    procedure, private :: getArrayCoord
    procedure, private :: getArrayCoordRegular
    procedure, private :: getPointCoordRegular
    procedure, private :: getPointCoordNonRegular
    procedure, private :: getArrayCoordNonRegular
    procedure :: initialize => initInterpolator
    procedure :: print => printInterpolator
    procedure, private :: test4D
    procedure, private :: interp4D
    procedure, private :: interp3D
    end type interpolator_class

    !Public access vars
    public :: interpolator_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that runs the chosen interpolator method on the given data.
    !> @param[in] self, state, bdata, time, var_dt, var_name, toInterp
    !---------------------------------------------------------------------------
    subroutine run(self, state, bdata, time, var_dt, var_name, toInterp)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:,:), intent(in) :: state
    type(background_class), intent(in) :: bdata
    real(prec), intent(in) :: time
    real(prec), dimension(:,:), intent(out) :: var_dt
    type(string), dimension(:), intent(out) :: var_name
    type(string), dimension(:), intent(in), optional :: toInterp
    logical :: interp
    real(prec) :: newtime
    class(*), pointer :: aField
    integer :: i
    type(string) :: outext

    real(prec), dimension(size(state,1)) :: xx, yy, zz
    logical, dimension(size(state,1)) :: outOfBounds
    real(prec) :: tt

    !Check field extents and what particles will be interpolated
    !interpolate each field to the correspoing slice in var_dt
    i = 1
    call bdata%fields%reset()                   ! reset list iterator
    do while(bdata%fields%moreValues())         ! loop while there are values
        interp = .true.
        aField => bdata%fields%currentValue()   ! get current value
        select type(aField)
        class is(scalar4d_field_class)          !4D interpolation is possible
            if (self%interpType == 1) then !linear interpolation in space and time
                if (present(toInterp)) then
                    if (.not.(any(toInterp == aField%name))) then
                        interp = .false.
                    end if
                end if
                if (interp) then
                    var_name(i) = aField%name
                    outOfBounds = .false.
                    xx = self%getArrayCoord(state(:,1), bdata, Globals%Var%lon, outOfBounds)
                    yy = self%getArrayCoord(state(:,2), bdata, Globals%Var%lat, outOfBounds)
                    zz = self%getArrayCoord(state(:,3), bdata, Globals%Var%level, outOfBounds)
                    tt = self%getPointCoordNonRegular(time, bdata, Globals%Var%time)
                    var_dt(:,i) = self%interp4D(xx, yy, zz, tt, outOfBounds, aField%field, size(aField%field,1), size(aField%field,2), size(aField%field,3), size(aField%field,4), size(state,1))
                end if
            end if !add more interpolation types here
        class is(scalar3d_field_class)          !3D interpolation is possible
            if (self%interpType == 1) then !linear interpolation in space and time
                if (present(toInterp)) then
                    if (.not.(any(toInterp == aField%name))) then
                        interp = .false.
                    end if
                end if
                if (interp) then
                    var_name(i) = aField%name
                    outOfBounds = .false.
                    xx = self%getArrayCoord(state(:,1), bdata, Globals%Var%lon, outOfBounds)
                    yy = self%getArrayCoord(state(:,2), bdata, Globals%Var%lat, outOfBounds)
                    tt = self%getPointCoordNonRegular(time, bdata, Globals%Var%time)
                    var_dt(:,i) = self%interp3D(xx, yy, tt, outOfBounds, aField%field, size(aField%field,1), size(aField%field,2), size(aField%field,3), size(state,1))
                end if
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
    !> @param[in] self, x, y, z, t, out, field, n_fv, n_cv, n_pv, n_tv, n_e
    !---------------------------------------------------------------------------
    function interp4D(self, x, y, z, t, out, field, n_fv, n_cv, n_pv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y, z                       !< 1-d. Array of particle component positions in array coordinates
    real(prec), intent(in) :: t                                           !< time to interpolate to in array coordinates
    logical, dimension(:), intent(in) :: out
    real(prec), dimension(n_fv, n_cv, n_pv, n_tv), intent(in) :: field    !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    integer, intent(in) :: n_fv, n_cv, n_pv, n_tv                         !< field dimensions
    integer, intent(in) :: n_e                                            !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, z0, x1, y1, z1
    real(prec), dimension(n_e) :: xd, yd, zd, c000, c100, c010, c110, c001
    real(prec), dimension(n_e) :: c101, c011, c111, c00, c10, c01, c11, c0, c1
    real(prec) :: td
    integer :: i, t0, t1
    real(prec), dimension(n_e) :: interp4D                                !< Field evaluated at x,y,z,t

    ! From x,y,z,t in array coordinates, find the the box inside the field where the particle is
    do concurrent(i=1:n_e, .not. out(i))
        x0(i) = floor(x(i))
        x1(i) = ceiling(x(i))
        y0(i) = floor(y(i))
        y1(i) = ceiling(y(i))
        z0(i) = floor(z(i))
        z1(i) = ceiling(z(i))
    end do

    t0 = floor(t)
    t1 = ceiling(t)

    ! If depth layer has one layer
    if (n_pv == 1) then
        z0 = 1
        z1 = 1
    end if

    xd = 0.
    yd = 0.
    zd = 0.
    td = 0.
    
    ! Compute the "normalized coordinates" of the particle inside the data field box
    where (x1 /= x0) xd = (x-x0)/(x1-x0)
    where (y1 /= y0) yd = (y-y0)/(y1-y0)
    where (z1 /= z0) zd = (z-z0)/(z1-z0)
    if (t1 /= t0) td = (t-t0)/(t1-t0)

    ! Interpolation on the first dimension and collapse it to a three dimension problem
    interp4D = 0.0
    
    do concurrent(i=1:n_e, .not. out(i))
        c000(i) = field(x0(i),y0(i),z0(i),t0)*(1.-xd(i)) + field(x1(i),y0(i),z0(i),t0)*xd(i) !y0x0z0t0!  y0x1z0t0
        c100(i) = field(x0(i),y1(i),z0(i),t0)*(1.-xd(i)) + field(x1(i),y1(i),z0(i),t0)*xd(i)
        c010(i) = field(x0(i),y0(i),z1(i),t0)*(1.-xd(i)) + field(x1(i),y0(i),z1(i),t0)*xd(i)
        c110(i) = field(x0(i),y1(i),z1(i),t0)*(1.-xd(i)) + field(x1(i),y1(i),z1(i),t0)*xd(i)

        c001(i) = field(x0(i),y0(i),z0(i),t1)*(1.-xd(i)) + field(x1(i),y0(i),z0(i),t1)*xd(i) !y0x0z0t0!  y0x1z0t0
        c101(i) = field(x0(i),y1(i),z0(i),t1)*(1.-xd(i)) + field(x1(i),y1(i),z0(i),t1)*xd(i)
        c011(i) = field(x0(i),y0(i),z1(i),t1)*(1.-xd(i)) + field(x1(i),y0(i),z1(i),t1)*xd(i)
        c111(i) = field(x0(i),y1(i),z1(i),t1)*(1.-xd(i)) + field(x1(i),y1(i),z1(i),t1)*xd(i)
        
    ! Interpolation on the second dimension and collapse it to a two dimension problem
        c00(i) = c000(i)*(1.-yd(i))+c100(i)*yd(i)
        c10(i) = c010(i)*(1.-yd(i))+c110(i)*yd(i)
        c01(i) = c001(i)*(1.-yd(i))+c101(i)*yd(i)
        c11(i) = c011(i)*(1.-yd(i))+c111(i)*yd(i)
    ! Interpolation on the third dimension and collapse it to a one dimension problem
        c0(i) = c00(i)*(1.-zd(i))+c10(i)*zd(i)
        c1(i) = c01(i)*(1.-zd(i))+c11(i)*zd(i)
    ! Interpolation on the time dimension and get the final result.
        interp4D(i) = c0(i)*(1.-td)+c1(i)*td
    end do
    
    end function interp4D

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> method to interpolate a particle position in a given data box based
    !> on array coordinates. 3d interpolation is a weighted average of 8
    !> neighbors. Consider the 4D domain between the 8 neighbors. The hypercube is
    !> divided into 4 sub-hypercubes by the point in question. The weight of each
    !> neighbor is given by the volume of the opposite sub-hypercube, as a fraction
    !> of the whole hypercube.
    !> @param[in] self, x, y, t, out, field, n_fv, n_cv, n_tv, n_e
    !---------------------------------------------------------------------------
    function interp3D(self, x, y, t, out, field, n_fv, n_cv, n_tv, n_e)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(n_e),intent(in):: x, y                        !< 1-d. Array of particle component positions in array coordinates
    real(prec), intent(in) :: t                                         !< time to interpolate to in array coordinates
    logical, dimension(:), intent(in) :: out
    real(prec), dimension(n_fv, n_cv, n_tv), intent(in) :: field        !< Field data with dimensions [n_fv,n_cv,n_pv,n_tv]
    integer, intent(in) :: n_fv, n_cv,n_tv                              !< field dimensions
    integer, intent(in) :: n_e                                          !< Number of particles to interpolate to
    integer, dimension(n_e) :: x0, y0, x1, y1
    real(prec), dimension(n_e) :: xd, yd, c00, c10, c01, c11
    real(prec), dimension(n_e) :: c0, c1
    real(prec) :: td
    integer :: i, t0, t1
    real(prec), dimension(n_e) :: interp3D                              !< Field evaluated at x,y,z,t

    ! From x,y,z,t in array coordinates, find the the box inside the field where the particle is
    do concurrent(i=1:n_e, .not. out(i))
        x0(i) = floor(x(i))
        x1(i) = ceiling(x(i))
        y0(i) = floor(y(i))
        y1(i) = ceiling(y(i))
    end do

    t0 = floor(t)
    t1 = ceiling(t)

    ! Compute the "normalized coordinates" of the particle inside the data field box
    xd = 0.
    yd = 0.
    td = 0.
    
    ! Compute the "normalized coordinates" of the particle inside the data field box
    where (x1 /= x0) xd = (x-x0)/(x1-x0)
    where (y1 /= y0) yd = (y-y0)/(y1-y0)
    if (t1 /= t0) td = (t-t0)/(t1-t0)


    interp3D = 0.0

    ! Interpolation on the first dimension and collapse it to a three dimension problem
    do concurrent(i=1:n_e, .not. out(i))
        c00(i) = field(x0(i),y0(i),t0)*(1.-xd(i)) + field(x1(i),y0(i),t0)*xd(i) !y0x0z0t0!  y0x1z0t0
        c10(i) = field(x0(i),y1(i),t0)*(1.-xd(i)) + field(x1(i),y1(i),t0)*xd(i)
        c01(i) = field(x0(i),y0(i),t1)*(1.-xd(i)) + field(x1(i),y0(i),t1)*xd(i)
        c11(i) = field(x0(i),y1(i),t1)*(1.-xd(i)) + field(x1(i),y1(i),t1)*xd(i)
    
    ! Interpolation on the second dimension and collapse it to a two dimension problem
        c0(i) = c00(i)*(1.-yd(i))+c10(i)*yd(i)
        c1(i) = c01(i)*(1.-yd(i))+c11(i)*yd(i)
    ! Interpolation on the time dimension and get the final result.
        interp3D(i) = c0(i)*(1.-td)+c1(i)*td
    end do

    end function interp3D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinates of a set of points, given a coordinate
    !> array.
    !> @param[in] self, xdata, bdata, dimName, out
    !---------------------------------------------------------------------------
    function getArrayCoord(self, xdata, bdata, dimName, out)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                 !< Background to use
    type(string), intent(in) :: dimName
    logical, dimension(:), intent(inout) :: out
    integer :: dim                                              !< corresponding background dimension
    real(prec), dimension(size(xdata)) :: getArrayCoord         !< coordinates in array index

    dim = bdata%getDimIndex(dimName)
    if (bdata%regularDim(dim)) getArrayCoord = self%getArrayCoordRegular(xdata, bdata, dim, out)
    if (.not.bdata%regularDim(dim)) getArrayCoord = self%getArrayCoordNonRegular(xdata, bdata, dim, out)

    end function getArrayCoord

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinates of a set of points, given a coordinate
    !> array. Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dim, out
    !---------------------------------------------------------------------------
    function getArrayCoordRegular(self, xdata, bdata, dim, out)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                 !< Background to use
    integer, intent(in) :: dim
    logical, dimension(:), intent(inout) :: out
    real(prec), dimension(size(xdata)) :: getArrayCoordRegular  !< coordinates in array index
    real(prec) :: minBound, maxBound, res
    if(size(bdata%dim(dim)%field) == 1) then
        getArrayCoordRegular = 1
        return
    end if
    minBound = bdata%dim(dim)%getFieldMinBound()
    maxBound = bdata%dim(dim)%getFieldMaxBound()
    res = abs(maxBound - minBound)/(size(bdata%dim(dim)%field)-1.0)
    getArrayCoordRegular = (xdata - minBound)/res + 1.0
    where (xdata < minBound) out = .true.
    where (xdata > maxBound) out = .true.
    end function getArrayCoordRegular


    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension.
    !> @param[in] self, xdata, bdata, dim, out
    ! !---------------------------------------------------------------------------
    function getArrayCoordNonRegular(self, xdata, bdata, dim, out)
    class(interpolator_class), intent(in) :: self
    real(prec), dimension(:), intent(in):: xdata                    !< Tracer coordinate component
    type(background_class), intent(in) :: bdata                     !< Background to use
    integer, intent(in) :: dim
    logical, dimension(:), intent(inout) :: out
    integer :: i                                                
    integer :: id, idx_1, idx_2                                 
    real(prec), dimension(size(xdata)) :: getArrayCoordNonRegular   !< coordinates in array index
    real(prec) :: minBound, maxBound

    if(size(bdata%dim(dim)%field) == 1) then
        getArrayCoordNonRegular = 1
        return
    end if
    
    getArrayCoordNonRegular = 1
    minBound = bdata%dim(dim)%getFieldMinBound()
    maxBound = bdata%dim(dim)%getFieldMaxBound()
    where (xdata < minBound) out = .true.
    where (xdata > maxBound) out = .true.
    do concurrent(id = 1:size(xdata), .not. out(id))
            do i = 2, size(bdata%dim(dim)%field)
                if (bdata%dim(dim)%field(i) >= xdata(id)) then
                    idx_1 = i-1
                    idx_2 = i
                    exit
                end if
            end do
        getArrayCoordNonRegular(id) = idx_1 + abs((xdata(id)-bdata%dim(dim)%field(idx_1))/(bdata%dim(dim)%field(idx_2)-bdata%dim(dim)%field(idx_1)))
    end do
    end function getArrayCoordNonRegular

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension.
    !> Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dimName
    !---------------------------------------------------------------------------
    function getPointCoordNonRegular(self, xdata, bdata, dimName)
    class(interpolator_class), intent(in) :: self
    real(prec), intent(in):: xdata              !< Tracer coordinate component
    type(background_class), intent(in) :: bdata !< Background to use
    type(string), intent(in) :: dimName
    integer :: dim                              !< corresponding background dimension
    real(prec) :: getPointCoordNonRegular       !< coordinates in array index
    type(string) :: outext
    integer :: i
    integer :: idx_1, idx_2, n_idx
    logical :: found

    found = .false.
    dim = bdata%getDimIndex(dimName)
    if(size(bdata%dim(dim)%field) == 1) then
        getPointCoordNonRegular = 1
        return
    end if
    n_idx = size(bdata%dim(dim)%field)
    do i = 2, n_idx
        if (bdata%dim(dim)%field(i) >= xdata) then
            idx_1 = i-1
            idx_2 = i
            found = .true.
            exit
        end if
    end do
    if (.not.found) then
        outext = '[Interpolator::getPointCoordNonRegular] Point not contained in "'//dimName//'" dimension, stoping'
        call Log%put(outext)
        stop
    end if
    getPointCoordNonRegular = idx_1 + abs((xdata-bdata%dim(dim)%field(idx_1))/(bdata%dim(dim)%field(idx_2)-bdata%dim(dim)%field(idx_1)))
    
    end function getPointCoordNonRegular

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns the array coordinate of a point, along a given dimension.
    !> Works only for regularly spaced data.
    !> @param[in] self, xdata, bdata, dimName, eta
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
    integer :: i                                                !< corresponding background dimension
    integer :: id,idx_1,idx_2,n_idx                                 !< corresponding background dimension
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
    logical, dimension(:), allocatable :: out
    real(prec) :: time
    integer :: npts, fieldDims
    real(prec) :: fieldVal
    fieldDims = 10
    fieldVal = 1.9
    allocate(field(fieldDims,fieldDims,fieldDims,fieldDims))
    field = fieldVal
    npts = 1
    allocate(xx(npts), yy(npts), zz(npts))
    allocate(out(npts))
    xx = 13.45
    yy = xx
    zz = xx
    time = 1
    print*, 'testing 4D interpolation, expected result is ', fieldVal
    out = .false.
    xx = self%interp4D(xx, yy, zz, time, out, field, fieldDims, fieldDims, fieldDims, fieldDims, npts)
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