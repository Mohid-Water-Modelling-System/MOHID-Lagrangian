    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : background
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines a background class that describes a solution from which to
    !> interpolate. A background object contains an arbitrary number of scalar or
    !> vectorial fields, in 2, 3 or 4D, indexed to labeled 1D fields of dimensions.
    !> The fields are stored in a linked list, enabling trivial iteration.
    !------------------------------------------------------------------------------

    module background_mod

    use common_modules
    use abstract_LinkedList_mod
    use fieldTypes_mod

    implicit none
    private

    type, extends(linkedlist) :: fieldsList_class !< List of fields class
    contains
    procedure :: print => print_fieldList
    procedure :: printCurrent => print_fieldListCurrent
    end type fieldsList_class

    type :: background_class        !< a background solution class
        logical :: initialized = .false.
        integer :: id = 0                                                       !< ID of the Background
        type(string) :: name                                                    !< Name of the Background
        type(box) :: extents                                                    !< shape::box that defines the extents of the Background solution
        type(scalar1d_field_class), allocatable, dimension(:) :: dim            !< Dimensions of the Background fields (time,lon,lat,level for example)
        type(fieldsList_class) :: fields                                        !< Linked list to store the fields in the Background
    contains
    procedure :: add => addField
    procedure :: getDimIndex
    procedure :: getDimExtents
    procedure :: append => appendFieldByTime
    procedure :: appendBackgroundByTime
    procedure :: getHyperSlab
    procedure, private :: getSlabDim
    procedure, private :: getPointDimIndexes
    procedure :: finalize => cleanBackground
    procedure, private :: setDims
    procedure, private :: setExtents
    procedure, private :: setID
    !clean by dimension range

    procedure :: test
    procedure :: print => printBackground
    end type background_class

    interface Background !< Constructor
    procedure constructor
    end interface

    !Public access vars
    public :: background_class, Background

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that adds a field to the Background object's field list
    !> @param[in] self, gfield
    !---------------------------------------------------------------------------
    subroutine addField(self, gfield)
    class(background_class), intent(inout) :: self
    type(generic_field_class), intent(in) :: gfield
    if (allocated(gfield%scalar1d%field)) call self%fields%add(gfield%scalar1d)
    if (allocated(gfield%scalar2d%field)) call self%fields%add(gfield%scalar2d)
    if (allocated(gfield%scalar3d%field)) call self%fields%add(gfield%scalar3d)
    if (allocated(gfield%scalar4d%field)) call self%fields%add(gfield%scalar4d)
    if (allocated(gfield%vectorial2d%field)) call self%fields%add(gfield%vectorial2d)
    if (allocated(gfield%vectorial3d%field)) call self%fields%add(gfield%vectorial3d)
    if (allocated(gfield%vectorial4d%field)) call self%fields%add(gfield%vectorial4d)
    end subroutine addField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Constructor for Background object
    !> @param[in] id, name, extents, dims
    !---------------------------------------------------------------------------
    function constructor(id, name, extents, dims)
    type(background_class) :: constructor
    integer, intent(in) :: id
    type(string), intent(in) :: name
    type(box), intent(in) :: extents
    type(scalar1d_field_class), dimension(:), intent(in) :: dims
    constructor%initialized = .true.
    call constructor%setID(id, name)
    call constructor%setExtents(extents)
    call constructor%setDims(dims)
    end function constructor

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the index of a given dimension by it's name
    !> @param[in] self, name, mandatory
    !---------------------------------------------------------------------------
    integer function getDimIndex(self, name, mandatory)
    class(background_class), intent(in) :: self
    type(string), intent(in) :: name
    logical, optional, intent(in) :: mandatory
    integer :: i
    type(string) :: outext
    logical found
    found = .false.
    do i=1, size(self%dim)
        if (self%dim(i)%name == name) then
            found = .true.
            getDimIndex = i
            return
        end if
    end do
    if (present(mandatory)) then
        if (mandatory) then            
            if (.not. found) then
                outext = '[background_class::getDimIndex]: Field dimensions dont contain a field called '// name //', stoping'
                call Log%put(outext)
                stop
            end if
        else 
            getDimIndex = MV_INT
        end if
    end if
    end function getDimIndex

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns two reals, min and max of a given dimension
    !> @param[in] self, name, mandatory
    !---------------------------------------------------------------------------
    function getDimExtents(self, name, mandatory) result(dimExtent)
    class(background_class), intent(in) :: self
    type(string), intent(in) :: name
    logical, optional, intent(in) :: mandatory
    logical :: mand
    integer :: i
    real(prec) :: dimExtent(2)
    mand = .true.
    if (present(mandatory)) mand = mandatory
    i = self%getDimIndex(name, mand)
    if ( i/= MV_INT) then
        dimExtent(1) = self%dim(i)%getFieldMinBound()
        dimExtent(2) = self%dim(i)%getFieldMaxBound()
    else
        dimExtent(1) = MV
        dimExtent(2) = MV
    end if
    end function getDimExtents

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> appends a given field to a matching field in the Background object's
    !> field list, if possible
    !> @param[in] self, gfield, dims, done
    !---------------------------------------------------------------------------
    subroutine appendFieldByTime(self, gfield, dims, done)
    class(background_class), intent(inout) :: self
    type(generic_field_class), intent(in) :: gfield
    type(scalar1d_field_class), dimension(:), intent(in) :: dims
    logical, intent(out) :: done
    real(prec), allocatable, dimension(:) :: newTime, oldTime, toAppendTime
    class(*), pointer :: aField
    type(string) :: outext
    integer :: i, j, posiTime
    
    done = .false.
    !check that dimensions are compatible
    !spacial dims must be the same, temporal must be consecutive
    if (size(self%dim) == size(dims)) then !ammount of dimensions is the same
        do i = 1, size(dims)
            j = self%getDimIndex(dims(i)%name) !getting the same dimension for the fields
            if (dims(i)%name /= Globals%Var%time) then  !dimension is not 'time'
                if (size(dims(i)%field) == size(self%dim(j)%field)) then !size of the arrays is the same
                    done = all(dims(i)%field == self%dim(j)%field)  !dimensions array is the same
                end if
            else
                posiTime = i
                done = all(dims(i)%field >= maxval(self%dim(j)%field)) !time arrays are consecutive or the same
            end if
        end do
    end if
    if (.not.done) return
    
    done = .false.
    !check that fields are compatible
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values to print
        aField => self%fields%currentValue()
        select type(aField)
        class is (generic_field_class)
            if (aField%compare(gfield)) then                
                !concatenate the fields on the background
                call aField%concatenate(gfield)
                !concatenate the 'time' dimension of the background
                i = self%getDimIndex(Globals%Var%time)
                allocate(oldTime, source = self%dim(i)%field)
                !allocate(newTime(size(oldTime) + size(dims(posiTime)%field)))
                newTime = [oldTime, dims(posiTime)%field]
                deallocate(self%dim(i)%field)
                allocate(self%dim(i)%field, source = newTime)
                call self%fields%reset()
                done = .true.
                return
            end if
            class default
            outext = '[Background::appendFieldByTime] Unexepected type of content, not a Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
    end do
    call self%fields%reset()               ! reset list iterator
    return
    
    end subroutine appendFieldByTime
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> appends a given field to a matching field in the Background object's
    !> field list, if possible
    !> @param[in] self, gfield, dims, done
    !---------------------------------------------------------------------------
    subroutine appendBackgroundByTime(self, bkg, done)
    class(background_class), intent(inout) :: self
    type(background_class), intent(in) :: bkg
    type(generic_field_class) :: gfield
    type(scalar1d_field_class), allocatable, dimension(:) :: dims
    logical, intent(out) :: done
    real(prec), allocatable, dimension(:) :: newTime, oldTime, toAppendTime
    class(*), pointer :: aField
    type(string) :: outext
    integer :: i, j, posiTime
    
    print*, 'appending backgrounds'
    done = .false.
    !check that dimensions are compatible
    !spacial dims must be the same, temporal must be consecutive
    if (size(self%dim) == size(bkg%dim)) then !ammount of dimensions is the same
        do i = 1, size(bkg%dim)
            j = self%getDimIndex(bkg%dim(i)%name) !getting the same dimension for the fields
            if (dims(i)%name /= Globals%Var%time) then  !dimension is not 'time'
                if (size(bkg%dim(i)%field) == size(self%dim(j)%field)) then !size of the arrays is the same
                    done = all(bkg%dim(i)%field == self%dim(j)%field)  !dimensions array is the same
                end if
            else
                posiTime = i
                done = all(bkg%dim(i)%field >= maxval(self%dim(j)%field)) !time arrays are consecutive or the same
            end if
        end do
    end if
    if (done) print*, 'field dims are compatible'
    if (.not.done) print*, 'field dims are not compatible'
    
    !if (.not.done) return
    
    !done = .false.
    !!check that fields are compatible
    !call self%fields%reset()               ! reset list iterator
    !do while(self%fields%moreValues())     ! loop while there are values to print
    !    aField => self%fields%currentValue()
    !    select type(aField)
    !    class is (generic_field_class)
    !        if (aField%compare(gfield)) then                
    !            !concatenate the fields on the background
    !            call aField%concatenate(gfield)
    !            !concatenate the 'time' dimension of the background
    !            i = self%getDimIndex(Globals%Var%time)
    !            allocate(oldTime, source = self%dim(i)%field)
    !            !allocate(newTime(size(oldTime) + size(dims(posiTime)%field)))
    !            newTime = [oldTime, dims(posiTime)%field]
    !            deallocate(self%dim(i)%field)
    !            allocate(self%dim(i)%field, source = newTime)
    !            call self%fields%reset()
    !            done = .true.
    !            return
    !        end if
    !        class default
    !        outext = '[Background::appendFieldByTime] Unexepected type of content, not a Field'
    !        call Log%put(outext)
    !        stop
    !    end select
    !    call self%fields%next()            ! increment the list iterator
    !end do
    !call self%fields%reset()               ! reset list iterator
    !return
    
    end subroutine appendBackgroundByTime
        
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns a background as a subset of another
    !> @param[in] self, domain, time
    !---------------------------------------------------------------------------
    type(background_class) function getHyperSlab(self, domain, time)
    class(background_class), intent(in) :: self
    type(box), intent(in) :: domain
    real(prec), intent(in), optional :: time(2)
    real(prec) :: ltime(2)
    type(scalar1d_field_class), allocatable, dimension(:) :: backgrounDims
    type(generic_field_class), allocatable, dimension(:) :: gfield
    class(*), pointer :: curr
    type(box) :: extents
    type(vector) :: pt
    real(prec), dimension(3,2) :: dimExtents
    integer, allocatable, dimension(:) :: llbound
    integer, allocatable, dimension(:) :: uubound
    integer :: temp_int
    type(string) :: outext
    integer :: i
    
    ltime = self%getDimExtents(Globals%Var%time)
    if (present(time)) ltime = time
    !finding index bounds of the slicing geometry
    allocate(llbound(size(self%dim)))
    allocate(uubound(size(self%dim)))
    llbound = self%getPointDimIndexes(domain%pt, ltime(1))
    uubound = self%getPointDimIndexes(domain%pt+domain%size, ltime(2))
    !slicing dimensions
    allocate(backgrounDims(size(self%dim)))    
    do i=1, size(self%dim)
        if (llbound(i) > uubound(i)) then !because We're not inverting the dimension and fields - Needs to be corrected 
            temp_int = llbound(i)
            llbound(i) = uubound(i)
            uubound(i) = temp_int
        end if
        llbound(i) = max(1, llbound(i)-1) !adding safety net to index bounds
        uubound(i) = min(uubound(i)+1, size(self%dim(i)%field))
        call backgrounDims(i)%initialize(self%dim(i)%name, self%dim(i)%units, 1, self%getSlabDim(i, llbound(i), uubound(i)))
    end do
    !slicing variables
    allocate(gfield(self%fields%getSize()))
    i=1
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values
        curr => self%fields%currentValue() ! get current value
        select type(curr)
        class is (field_class)
            gfield(i) = curr%getFieldSlice(llbound, uubound)        
        class default
        outext = '[background_class::getHyperSlab] Unexepected type of content, not a scalar Field'
        call Log%put(outext)
        stop
        end select
        call self%fields%next()            ! increment the list iterator
        i = i+1
    end do
    call self%fields%reset()               ! reset list iterator
    !creating bounding box
    dimExtents = 0.0
    do i = 1, size(backgrounDims)
        if (backgrounDims(i)%name == Globals%Var%lon) then
            dimExtents(1,1) = backgrounDims(i)%getFieldMinBound()
            dimExtents(1,2) = backgrounDims(i)%getFieldMaxBound()
        else if (backgrounDims(i)%name == Globals%Var%lat) then
            dimExtents(2,1) = backgrounDims(i)%getFieldMinBound()
            dimExtents(2,2) = backgrounDims(i)%getFieldMaxBound()
        else if (backgrounDims(i)%name == Globals%Var%level) then
            dimExtents(3,1) = backgrounDims(i)%getFieldMinBound()
            dimExtents(3,2) = backgrounDims(i)%getFieldMaxBound()
        end if
    end do
    extents%pt = dimExtents(1,1)*ex + dimExtents(2,1)*ey + dimExtents(3,1)*ez
    pt = dimExtents(1,2)*ex + dimExtents(2,2)*ey + dimExtents(3,2)*ez
    extents%size = pt - extents%pt
    !creating the sliced background
    getHyperSlab = constructor(1, self%name, extents, backgrounDims)
    do i=1, size(gfield)
        call getHyperSlab%add(gfield(i))
    end do
    
    end function getHyperSlab
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns the indexes of the dims for a given point
    !> @param[in] self, pt, time
    !---------------------------------------------------------------------------
    function getPointDimIndexes(self, pt, time)
    class(background_class), intent(in) :: self
    type(vector), intent(in) :: pt
    real(prec), intent(in) :: time
    integer, allocatable, dimension(:) :: getPointDimIndexes
    integer :: i
    allocate(getPointDimIndexes(size(self%dim)))
    do i= 1, size(self%dim)
        if (self%dim(i)%name == Globals%Var%lon) getPointDimIndexes(i) = self%dim(i)%getFieldNearestIndex(pt%x)
        if (self%dim(i)%name == Globals%Var%lat) getPointDimIndexes(i) = self%dim(i)%getFieldNearestIndex(pt%y)
        if (self%dim(i)%name == Globals%Var%level) getPointDimIndexes(i) = self%dim(i)%getFieldNearestIndex(pt%z)
        if (self%dim(i)%name == Globals%Var%time)  getPointDimIndexes(i) = self%dim(i)%getFieldNearestIndex(time)
    end do
    end function getPointDimIndexes
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns an array witn a sliced dimension
    !> @param[in] self, domain, time
    !---------------------------------------------------------------------------
    function getSlabDim(self, numDim, llbound, uubound)
    class(background_class), intent(in) :: self
    integer, intent(in) :: numDim, llbound, uubound
    real(prec), allocatable, dimension(:) :: getSlabDim
    allocate(getSlabDim, source = self%dim(numDim)%field(llbound:uubound))    
    end function getSlabDim    
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that cleans all data in the Background object
    !---------------------------------------------------------------------------
    subroutine cleanBackground(self)
    class(background_class), intent(inout) :: self
    class(*), pointer :: curr
    type(string) :: outext
    !self%initialized = .false.
    self%id = MV_INT
    self%name = ''
    deallocate(self%dim)
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values
        curr => self%fields%currentValue() ! get current value
        select type(curr)
        class is (scalar1d_field_class)
            call curr%finalize()
        class is (scalar2d_field_class)
            call curr%finalize()
        class is (scalar3d_field_class)
            call curr%finalize()
        class is (scalar4d_field_class)
            call curr%finalize()
        class default
        outext = '[background_class::cleanBackground] Unexepected type of content, not a scalar Field'
        call Log%put(outext)
        stop
        end select
        call self%fields%next()            ! increment the list iterator
    end do
    call self%fields%reset()               ! reset list iterator
    call self%fields%finalize()    
    end subroutine cleanBackground

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and sets the dimensions of the Background object
    !> @param[in] self, dims
    !---------------------------------------------------------------------------
    subroutine setDims(self, dims)
    class(background_class), intent(inout) :: self
    type(scalar1d_field_class), dimension(:), intent(in) :: dims
    allocate(self%dim, source = dims)
    end subroutine setDims

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the extents (bounding box) of the Background object
    !> @param[in] self, bbox
    !---------------------------------------------------------------------------
    subroutine setExtents(self, bbox)
    class(background_class), intent(inout) :: self
    type(box), intent(in) :: bbox
    self%extents = bbox
    end subroutine setExtents

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the ID and name of the Background object
    !> @param[in] self, id, name
    !---------------------------------------------------------------------------
    subroutine setID(self, id, name)
    class(background_class), intent(inout) :: self
    integer, intent(in) :: id
    type(string), intent(in) :: name
    self%id = id
    self%name = name
    end subroutine setID

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> A class 'unit' test for the background_class
    !---------------------------------------------------------------------------
    subroutine test(self)
    class(background_class), intent(inout) :: self
    type(background_class) :: background1
    type(generic_field_class) :: gfield1, gfield2, gfield3
    real(prec), allocatable, dimension(:) :: field1
    real(prec), allocatable, dimension(:,:) :: field2
    type(vector), allocatable, dimension(:,:,:) :: field3
    real(prec), allocatable, dimension(:,:,:) :: any_field_3d
    type(string) :: name1, name2, name3, bname
    type(string) :: units1, units2, units3
    type(box) :: backgroundbbox
    type(scalar1d_field_class), allocatable, dimension(:) :: backgroundims
    !generating fields
    !inquire nc dimensions
    !allocate approptiate real matrix (1d, 2d, 3d ..)
    allocate(field1(50))
    allocate(field2(20,60))
    allocate(field3(2,3,4))
    !inquire nc field name and units
    name1 = 'testfield1d'
    name2 = 'testfield2d'
    name3 = 'testfield3d'
    units1 = 'm/s'
    units2 = 'km'
    units3 = 'ms-1'
    !nc_get_var - stores data in allocated matrix
    !put data in generic field
    call gfield1%initialize(name1, units1, field1)
    call gfield2%initialize(name2, units2, field2)
    call gfield3%initialize(name3, units3, field3)
    !assembling our Background
    !nc inquire dimensions names
    bname = 'TestBackground'
    name1 = 'lon'
    name2 = 'lat'
    !make background bounding box (this should come from the block + a margin)
    backgroundbbox%pt = 1*ex + 2*ey + 3*ez
    backgroundbbox%size = 4*ex + 5*ey + 6*ez
    !allocate space for the dimensions vectors of the background
    allocate(backgroundims(2))
    !inquire dimensions units, and data
    call backgroundims(1)%initialize(name1,units2,1, field1)
    call backgroundims(2)%initialize(name2,units2,1, field1)
    !construct background
    background1 = Background(5, bname, backgroundbbox, backgroundims)
    call background1%add(gfield1)
    call background1%add(gfield2)
    call background1%add(gfield3)
    call background1%print()
    end subroutine test

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the Background object
    !---------------------------------------------------------------------------
    subroutine printBackground(self)
    class(background_class), intent(inout) :: self
    type(string) :: outext, t
    integer :: i
    t = self%id
    outext = 'Background['//t//', '//self%name//'] is a'
    call Log%put(outext,.false.)
    call Geometry%print(self%extents)
    outext = 'The dimensions fields are:'
    call Log%put(outext,.false.)
    do i=1, size(self%dim)
        call self%dim(i)%print()
    end do
    outext = 'The data fields are:'
    call Log%put(outext,.false.)
    call self%fields%print()
    end subroutine printBackground

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the links of the list
    !---------------------------------------------------------------------------
    subroutine print_fieldList(this)
    class(fieldsList_class), intent(in) :: this
    call this%reset()               ! reset list iterator
    do while(this%moreValues())     ! loop while there are values to print
        call this%printCurrent()
        call this%next()            ! increment the list iterator
    end do
    call this%reset()               ! reset list iterator
    end subroutine print_fieldList

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the current link of the list
    !---------------------------------------------------------------------------
    subroutine print_fieldListCurrent(this)
    class(fieldsList_class), intent(in) :: this
    class(*), pointer :: curr
    type(string) :: outext
    curr => this%currentValue() ! get current value
    select type(curr)
    class is (field_class)
        call curr%print()
        class default
        outext = '[fieldsList_class::print] Unexepected type of content, not a Field'
        call Log%put(outext)
        stop
    end select
    end subroutine print_fieldListCurrent


    end module background_mod