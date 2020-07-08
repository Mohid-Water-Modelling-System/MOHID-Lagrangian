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
        logical, allocatable, dimension(:) :: regularDim                        !< Flag that indicates if the respective dimension is regular or irregular
        type(fieldsList_class) :: fields                                        !< Linked list to store the fields in the Background
        type(stringList_class) :: variables
    contains
    procedure :: add => addField
    procedure :: getDimIndex
    procedure :: getDimExtents
    procedure :: append => appendBackgroundByTime
    procedure :: getHyperSlab
    procedure :: ShedMemory
    procedure :: makeLandMaskField
    procedure :: makeResolutionField
    procedure :: makeBathymetryField
    procedure :: copy
    procedure :: hasVars
    procedure, private :: getSlabDim
    procedure, private :: getPointDimIndexes
    procedure :: finalize => cleanBackground
    procedure, private :: cleanFields
    procedure, private :: setDims
    procedure, private :: setExtents
    procedure, private :: setID
    procedure :: check
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
    logical :: added
    added = .false.
    if (allocated(gfield%scalar1d%field)) then
        call self%fields%add(gfield%scalar1d)
        added = .true.
    end if
    if (allocated(gfield%scalar2d%field)) then
        call self%fields%add(gfield%scalar2d)
        added = .true.
    end if
    if (allocated(gfield%scalar3d%field)) then
        call self%fields%add(gfield%scalar3d)
        added = .true.
    end if
    if (allocated(gfield%scalar4d%field)) then
        call self%fields%add(gfield%scalar4d)
        added = .true.
    end if
    !if (allocated(gfield%vectorial2d%field)) call self%fields%add(gfield%vectorial2d)
    !if (allocated(gfield%vectorial3d%field)) call self%fields%add(gfield%vectorial3d)
    !if (allocated(gfield%vectorial4d%field)) call self%fields%add(gfield%vectorial4d)
    if (added) then
        if(self%variables%notRepeated(gfield%name)) call self%variables%add(gfield%name)
    end if
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
    !> @param[in] self, bkg, done
    !---------------------------------------------------------------------------
    subroutine appendBackgroundByTime(self, bkg, done)
    class(background_class), intent(inout) :: self
    type(background_class), intent(in) :: bkg
    type(generic_field_class) :: tempGField, tempGField2
    type(generic_field_class), allocatable, dimension(:) :: gField
    logical, intent(out) :: done
    real(prec), allocatable, dimension(:) :: newTime
    class(*), pointer :: aField, bField
    type(string) :: outext, name, units
    integer :: i, j, k
    logical, allocatable, dimension(:) :: usedTime

    done = .false.
    !check that dimensions are compatible
    !spacial dims must be the same, temporal must be consecutive
    if (size(self%dim) == size(bkg%dim)) then !ammount of dimensions is the same
        do i = 1, size(bkg%dim)
            j = self%getDimIndex(bkg%dim(i)%name) !getting the same dimension for the fields
            if (bkg%dim(i)%name /= Globals%Var%time) then  !dimension is not 'time'
                if (size(bkg%dim(i)%field) == size(self%dim(j)%field)) then !size of the arrays is the same
                    done = all(bkg%dim(i)%field == self%dim(j)%field)  !dimensions array is the same
                end if
            else
                !done = all(bkg%dim(i)%field >= maxval(self%dim(j)%field)) !time arrays are consecutive or the same
                allocate(newTime, source = self%dim(j)%field)
                call Utils%appendArraysUniqueReal(newTime, bkg%dim(i)%field, usedTime)
                !check if new time dimension is consistent (monotonic and not repeating)
                done = all(newTime(2:)-newTime(1:size(newTime)-1) > 0)
                name = self%dim(j)%name
                units = self%dim(j)%units
                call self%dim(j)%finalize()
                call self%dim(j)%initialize(name, units, 1, newTime)
            end if
        end do
    end if
    if (.not.done) return

    allocate(gField(self%fields%getSize()))
    i=1
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values to process
        aField => self%fields%currentValue()
        select type(aField)
        class is (field_class)
            call gField(i)%getGField(aField)
            call bkg%fields%reset()
            do while(bkg%fields%moreValues())
                bField => bkg%fields%currentValue()
                select type(bField)
                class is (field_class)
                    if (bField%name == aField%name) then
                        call tempGField%getGField(bField)
                        !append the new time instances of the field
                        call gField(i)%concatenate(tempGField, usedTime)
                        call tempGField%finalize()
                        exit
                    end if
                end select
                call bkg%fields%next()
                nullify(bField)
            end do
            call bkg%fields%reset()
            class default
            outext = '[Background::appendBackgroundByTime] Unexepected type of content, not a Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
        i = i+1
        nullify(aField)
    end do
    call self%fields%reset()               ! reset list iterator

    call self%cleanFields()
    call self%fields%finalize()

    do i=1, size(gField)
        call self%add(gField(i))
    end do

    if(.not.self%check()) then
        outext = '[Background::appendBackgroundByTime]: non-conformant Background, stoping '
        call Log%put(outext)
        stop
    end if

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
        llbound(i) = max(1, llbound(i)-2) !adding safety net to index bounds
        uubound(i) = min(uubound(i)+2, size(self%dim(i)%field))
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
        nullify(curr)
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

    if(.not.self%check()) then
        outext = '[Background::getHyperSlab]: non-conformant Background, stoping '
        call Log%put(outext)
        stop
    end if

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
    !> @param[in] self, numDim, llbound, uubound
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
    !> Method that cleans data in the Background object that is already
    !> out of scope
    !---------------------------------------------------------------------------
    subroutine ShedMemory(self)
    class(background_class), intent(inout) :: self
    integer, allocatable, dimension(:) :: llbound
    integer, allocatable, dimension(:) :: uubound
    real(prec), allocatable, dimension(:) :: newTime
    type(string) :: outext, name, units
    type(generic_field_class), allocatable, dimension(:) :: gField
    class(*), pointer :: aField
    logical :: done
    integer :: i, j

    if (self%initialized) then
        done = .false.
        allocate(llbound(size(self%dim)))
        allocate(uubound(size(self%dim)))
        !select valid time coordinate array elements
        do i= 1, size(self%dim)
            llbound(i) = 1
            uubound(i) = size(self%dim(i)%field)

            if (self%dim(i)%name == Globals%Var%time) then
                if (self%dim(i)%field(1) < Globals%SimTime%CurrTime - Globals%Parameters%BufferSize) then
                    llbound(i) = self%dim(i)%getFieldNearestIndex(Globals%SimTime%CurrTime - Globals%Parameters%BufferSize/2.0)
                    if (llbound(i) == self%dim(i)%getFieldNearestIndex(Globals%SimTime%CurrTime)) llbound(i) = llbound(i) - 1
                    if (llbound(i) > 2) then
                        uubound(i) = size(self%dim(i)%field)
                        j=size(self%dim(i)%field(llbound(i):uubound(i)))
                        allocate(newTime(j))
                        newTime = self%dim(i)%field(llbound(i):uubound(i))
                        !allocate(newTime, source = self%getSlabDim(i, llbound(i), uubound(i))) !gfortran doesn't like this...
                        name = self%dim(i)%name
                        units = self%dim(i)%units
                        call self%dim(i)%finalize()
                        call self%dim(i)%initialize(name, units, 1, newTime)
                        done  = .true.
                    end if
                    exit
                end if
            end if
        end do
        if (done) then
            !slice variables accordingly
            allocate(gField(self%fields%getSize()))
            i=1
            call self%fields%reset()               ! reset list iterator
            do while(self%fields%moreValues())     ! loop while there are values to process
                aField => self%fields%currentValue()
                select type(aField)
                class is (field_class)
                    gfield(i) = aField%getFieldSlice(llbound, uubound)
                    class default
                    outext = '[Background::ShedMemory] Unexepected type of content, not a Field'
                    call Log%put(outext)
                    stop
                end select
                call self%fields%next()            ! increment the list iterator
                i = i+1
                nullify(aField)
            end do
            call self%fields%reset()               ! reset list iterator

            call self%cleanFields()
            call self%fields%finalize()

            do i=1, size(gField)
                call self%add(gField(i))
            end do

            if(.not.self%check()) then
                outext = '[Background::ShedMemory]: non-conformant Background, stoping '
                call Log%put(outext)
                stop
            end if
        end if
    end if

    end subroutine ShedMemory

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to use a stored binary field to make a land mask used for land
    !> interaction: bed settling, beaching, ...
    !---------------------------------------------------------------------------
    subroutine makeLandMaskField(self)
    class(background_class), intent(inout) :: self
    class(*), pointer :: curr
    logical, allocatable, dimension(:,:,:) :: shiftleftlon3d, shiftuplat3d, shiftrigthlon3d, shiftdownlat3d, beach3d
    logical, allocatable, dimension(:,:,:,:) :: shiftleftlon4d, shiftuplat4d, shiftrigthlon4d, shiftdownlat4d, shiftUpLevel, shiftDownLevel, beach4d, bed4d
    type(string) :: outext
    integer :: dimIndx, k
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values
        curr => self%fields%currentValue() ! get current value
        select type(curr)
        class is (scalar3d_field_class)
            if (curr%name == Globals%Var%landIntMask) then
                allocate(shiftleftlon3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                allocate(shiftuplat3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                allocate(shiftrigthlon3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                allocate(shiftdownlat3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                allocate(beach3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                shiftleftlon3d = .false.
                shiftleftlon3d(:size(curr%field,1)-1,:,:) = abs(curr%field(:size(curr%field,1)-1,:,:) - curr%field(2:,:,:)) /= 0.0
                shiftrigthlon3d = .false.
                shiftrigthlon3d(2:,:,:) = abs(curr%field(2:,:,:) - curr%field(:size(curr%field,1)-1,:,:)) /= 0.0
                shiftuplat3d = .false.
                shiftuplat3d(:,:size(curr%field,2)-1,:) = abs(curr%field(:,:size(curr%field,2)-1,:) - curr%field(:,2:,:)) /= 0.0
                shiftdownlat3d = .false.
                shiftdownlat3d(:, 2:,:) = abs(curr%field(:,2:,:) - curr%field(:,:size(curr%field,2)-1,:)) /= 0.0
                beach3d = .false.
                beach3d = shiftleftlon3d .or. shiftrigthlon3d .or. shiftuplat3d .or. shiftdownlat3d !colapsing all the shifts
                beach3d = beach3d .and. (curr%field == Globals%Mask%waterVal) !just points that were already wet
                where(beach3d) curr%field = Globals%Mask%beachVal
                !searching for areas that are beach and land and mark them as beach only 
                do k=1, size(curr%field,3)
                    beach3d(:,:,1) = beach3d(:,:,1) .or. beach3d(:,:,k)
                    beach3d(:,:,k) = beach3d(:,:,1)
                end do
                where ((curr%field == Globals%Mask%landVal) .and. beach3d) curr%field = Globals%Mask%beachVal
            end if
        class is (scalar4d_field_class)
            if (curr%name == Globals%Var%landIntMask) then
                allocate(shiftleftlon4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(shiftuplat4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(shiftrigthlon4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(shiftdownlat4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(shiftUpLevel(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(shiftDownLevel(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(beach4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(bed4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                shiftleftlon4d = .false.
                shiftleftlon4d(:size(curr%field,1)-1,:,:,:) = abs(curr%field(:size(curr%field,1)-1,:,:,:) - curr%field(2:,:,:,:)) /= 0.0
                shiftrigthlon4d = .false.
                shiftrigthlon4d(2:,:,:,:) = abs(curr%field(2:,:,:,:) - curr%field(:size(curr%field,1)-1,:,:,:)) /= 0.0
                shiftdownlat4d = .false.
                shiftdownlat4d(:,:size(curr%field,2)-1,:,:) = abs(curr%field(:,:size(curr%field,2)-1,:,:) - curr%field(:,2:,:,:)) /= 0.0
                shiftuplat4d = .false.
                shiftuplat4d(:,2:,:,:) = abs(curr%field(:,2:,:,:) - curr%field(:,:size(curr%field,2)-1,:,:)) /= 0.0
                shiftUpLevel = .false.
                shiftUpLevel(:,:,2:,:) = abs(curr%field(:,:,2:,:) - curr%field(:,:,:size(curr%field,3)-1,:)) /= 0.0
                shiftDownLevel = .false.
                shiftDownLevel(:,:,:size(curr%field,3)-1,:) = abs(curr%field(:,:,:size(curr%field,3)-1,:) - curr%field(:,:,2:,:)) /= 0.0
                beach4d = .false.
                beach4d = shiftleftlon4d .or. shiftrigthlon4d .or. shiftuplat4d .or. shiftdownlat4d .or. shiftDownLevel .or. shiftUpLevel !colapsing all the shifts
                beach4d = beach4d .and. (curr%field == Globals%Mask%waterVal) !just points that were already wet
                bed4d = beach4d
                dimIndx = self%getDimIndex(Globals%Var%level)
                dimIndx = minloc(abs(self%dim(dimIndx)%field - Globals%Constants%BeachingLevel),1)
                beach4d(:,:,:dimIndx,:) = .false. !this must be above a certain level only
                bed4d(:,:,dimIndx:,:) = .false.   !bellow a certain level
                where(beach4d) curr%field = Globals%Mask%beachVal
                where(bed4d) curr%field = Globals%Mask%bedVal
                !searching for areas that are beach and land and mark them as beach only 
                do k=1, size(curr%field,4)
                    beach4d(:,:,:,1) = beach4d(:,:,:,1) .or. beach4d(:,:,:,k)
                    beach4d(:,:,:,k) = beach4d(:,:,:,1)
                end do
                where ((curr%field == Globals%Mask%landVal) .and. beach4d) curr%field = Globals%Mask%beachVal
            end if
            class default
            outext = '[background_class::makeLandMaskField] Unexepected type of content, not a 3D or 4D scalar Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
        nullify(curr)
    end do
    call self%fields%reset()               ! reset list iterator

    end subroutine makeLandMaskField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to use a stored binary field to make a resolution proxy for the
    !> underlying mesh
    !---------------------------------------------------------------------------
    subroutine makeResolutionField(self)
    class(background_class), intent(inout) :: self
    class(*), pointer :: curr
    real(prec), allocatable, dimension(:,:,:) :: xx3d, yy3d
    real(prec), allocatable, dimension(:,:,:,:) :: xx4d, yy4d, zz4d
    type(string) :: outext
    integer :: xIndx, yIndx, zIndx, i, j, k, t
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values
        curr => self%fields%currentValue() ! get current value
        select type(curr)
        class is (scalar3d_field_class)
            if (curr%name == Globals%Var%resolution) then
                allocate(xx3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                allocate(yy3d(size(curr%field,1), size(curr%field,2), size(curr%field,3)))
                xIndx = self%getDimIndex(Globals%Var%lon)
                yIndx = self%getDimIndex(Globals%Var%lat)
                do j=1, size(xx3d,2)
                    xx3d(:,j,1) = Utils%geo2m(abs(self%dim(xIndx)%field(:size(curr%field,1)-1) - self%dim(xIndx)%field(2:)), self%dim(yIndx)%field(j), .false.)
                end do
                yy3d(1,:,1) = Utils%geo2m(abs(self%dim(yIndx)%field(:size(curr%field,2)-1) - self%dim(yIndx)%field(2:)), self%dim(yIndx)%field(1), .true.)
                do i=2, size(yy3d,1)
                    yy3d(i,:,1) = yy3d(1,:,1)
                end do
                do k=2, size(curr%field,3)
                    yy3d(:,:,k) = yy3d(:,:,1)
                    xx3d(:,:,k) = xx3d(:,:,1)
                end do
                curr%field = (xx3d + yy3d)/2.0
            end if
        class is (scalar4d_field_class)
            if (curr%name == Globals%Var%resolution) then
                allocate(xx4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(yy4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(zz4d(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                xx4d = -99999.0
                xIndx = self%getDimIndex(Globals%Var%lon)
                yIndx = self%getDimIndex(Globals%Var%lat)
                zIndx = self%getDimIndex(Globals%Var%level)
                do j=1, size(xx4d,2)
                    xx4d(:,j,1,1) = Utils%geo2m(abs(self%dim(xIndx)%field(:size(curr%field,1)-1) - self%dim(xIndx)%field(2:)), self%dim(yIndx)%field(j), .false.)
                end do
                yy4d(1,:,1,1) = Utils%geo2m(abs(self%dim(yIndx)%field(:size(curr%field,2)-1) - self%dim(yIndx)%field(2:)), self%dim(yIndx)%field(1), .true.)
                do i=2, size(yy4d,1)
                    yy4d(i,:,1,1) = yy4d(1,:,1,1)
                end do
                do k=2, size(curr%field,3)
                    xx4d(:,:,k,1) = xx4d(:,:,1,1)
                    yy4d(:,:,k,1) = yy4d(:,:,1,1)
                end do
                zz4d(1,1,2:,1) = abs(self%dim(zIndx)%field(:size(curr%field,3)-1) - self%dim(zIndx)%field(2:))
                do i=2, size(zz4d,1)
                    zz4d(i,1,:,1) = zz4d(1,1,:,1)
                end do
                do j=2, size(zz4d,2)
                    zz4d(:,j,:,1) = zz4d(:,1,:,1)
                end do
                do t=2, size(curr%field,4)
                    xx4d(:,:,:,t) = xx4d(:,:,:,1)
                    yy4d(:,:,:,t) = yy4d(:,:,:,1)
                    zz4d(:,:,:,t) = zz4d(:,:,:,1)
                end do
                curr%field = xx4d!(xx4d + yy4d + zz4d)/3.0
            end if
            class default
            outext = '[background_class::makeResolutionField] Unexepected type of content, not a 3D or 4D scalar Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
        nullify(curr)
    end do
    call self%fields%reset()               ! reset list iterator

    end subroutine makeResolutionField
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to use a stored binary field to make a bathymetry field - depends on fill values
    !---------------------------------------------------------------------------
    subroutine makeBathymetryField(self)
    class(background_class), intent(inout) :: self
    class(*), pointer :: curr
    logical, allocatable, dimension(:,:,:,:) :: shiftUpLevel
    real(prec), allocatable, dimension(:,:,:,:) :: bathymetry
    type(string) :: outext
    integer :: dimIndx, i, j, t, k
    
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values
        curr => self%fields%currentValue() ! get current value
        select type(curr)
        class is (scalar3d_field_class)
            if (curr%name == Globals%Var%bathymetry) then
                curr%field = MV
            end if
        class is (scalar4d_field_class)               
            if (curr%name == Globals%Var%bathymetry) then         
                allocate(shiftUpLevel(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                allocate(bathymetry(size(curr%field,1), size(curr%field,2), size(curr%field,3), size(curr%field,4)))
                dimIndx = self%getDimIndex(Globals%Var%level)
                bathymetry = self%dim(dimIndx)%field(size(curr%field,3))
                if (size(curr%field,3) == 1) then   ! if netcdf data has one depth layer.
                    curr%field = bathymetry
                else
                    dimIndx = self%getDimIndex(Globals%Var%level)
                    bathymetry = self%dim(dimIndx)%field(size(curr%field,3))
                    shiftUpLevel = .false.
                    shiftUpLevel(:,:,2:,:) = abs(curr%field(:,:,2:,:) - curr%field(:,:,:size(curr%field,3)-1,:)) == 0.0
                    shiftUpLevel(:,:,1,:) = abs(curr%field(:,:,2,:) - curr%field(:,:,1,:)) == 0.0
                    shiftUpLevel(:,:,size(curr%field,3),:) = abs(curr%field(:,:,size(curr%field,3),:) - curr%field(:,:,size(curr%field,3)-1,:)) == 0.0
                    do t=1, size(curr%field,4)
                        do j=1, size(curr%field,2)
                            do i=1, size(curr%field,1)
                                do k=1, size(shiftUpLevel(i,j,:,t))
                                    if (all(shiftUpLevel(i,j,k:size(curr%field,3),t))) then 
                                        bathymetry(i,j,:,t) = self%dim(dimIndx)%field(k)
                                        exit
                                    end if
                                end do
                            end do
                        end do
                    end do
                    curr%field = bathymetry
                end if
            end if
            class default
            outext = '[background_class::makeBathymetryField] Unexepected type of content, not a 3D or 4D scalar Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
        nullify(curr)
    end do
    call self%fields%reset()               ! reset list iterator

    end subroutine makeBathymetryField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that copies all data in another Background object
    !---------------------------------------------------------------------------
    subroutine copy(self, bkg)
    class(background_class), intent(inout) :: self
    type(background_class), intent(in) :: bkg
    integer :: i
    class(*), pointer :: aField
    type(generic_field_class) :: gField
    type(string) :: outext

    if(self%initialized) call self%finalize()
    self%initialized = bkg%initialized
    self%id = bkg%id
    self%name = bkg%name
    self%extents%pt = bkg%extents%pt
    self%extents%size = bkg%extents%size
    if (allocated(bkg%dim)) then
        allocate(self%dim(size(bkg%dim)))
        do i=1, size(bkg%dim)
            call self%dim(i)%initialize(bkg%dim(i)%name, bkg%dim(i)%units, 1, bkg%dim(i)%field)
        end do
    end if
    if (allocated(bkg%regularDim)) then
        allocate(self%regularDim(size(bkg%regularDim)))
        self%regularDim = bkg%regularDim
    end if
    call bkg%fields%reset()               ! reset list iterator
    do while(bkg%fields%moreValues())     ! loop while there are values to process
        aField => bkg%fields%currentValue()
        select type(aField)
        class is (field_class)
            call gField%getGField(aField)
            class default
            outext = '[Background::copy] Unexepected type of content, not a Field'
            call Log%put(outext)
            stop
        end select
        call self%add(gField)
        call gField%finalize()
        nullify(aField)
        call bkg%fields%next()            ! increment the list iterator
    end do
    call self%fields%reset()               ! reset list iterator

    end subroutine copy
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns true if a given set of variable names in a string array is 
    !> contained in the variable list of the background
    !> @param[in] self, reqVars
    !---------------------------------------------------------------------------
    logical function hasVars(self, reqVars)
    class(background_class), intent(in) :: self
    type(string), dimension(:), intent(in) :: reqVars
    integer :: i
    hasVars = .true.
    do i=1, size(reqVars)
        if (self%variables%notRepeated(reqVars(i))) then
            hasVars = .false.
            return
        end if
    end do
    end function hasVars

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that cleans all data in the Background object
    !---------------------------------------------------------------------------
    subroutine cleanBackground(self)
    class(background_class), intent(inout) :: self
    integer :: i
    self%initialized = .false.
    self%id = MV_INT
    self%name = ''
    self%extents%size = 0.0
    if (allocated(self%dim)) then
        do i=1, size(self%dim)
            call self%dim(i)%finalize()
        end do
        deallocate(self%dim)
    end if
    if (allocated(self%regularDim)) deallocate(self%regularDim)
    call self%cleanFields()
    call self%fields%finalize()
    call self%variables%finalize()
    end subroutine cleanBackground

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that cleans all data in the Background object fields
    !---------------------------------------------------------------------------
    subroutine cleanFields(self)
    class(background_class), intent(inout) :: self
    class(*), pointer :: curr
    type(string) :: outext
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
            outext = '[background_class::cleanFields] Unexepected type of content, not a scalar Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
        nullify(curr)
    end do
    call self%fields%reset()               ! reset list iterator
    end subroutine cleanFields

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and sets the dimensions of the Background object
    !> @param[in] self, dims
    !---------------------------------------------------------------------------
    subroutine setDims(self, dims)
    class(background_class), intent(inout) :: self
    type(scalar1d_field_class), dimension(:), intent(in) :: dims
    real(prec), allocatable, dimension(:) :: rest
    integer :: i
    real(prec) ::fmin, fmax, eta,dreg
    integer :: f_1,f_N
    allocate(self%dim, source = dims)
    allocate(self%regularDim(size(dims)))
    self%regularDim = .false.
    do i=1, size(dims)
        fmin = self%dim(i)%getFieldMinBound() 
        fmax = self%dim(i)%getFieldMaxBound()
        eta = (fmax-fmin)/(10.0*size(self%dim(i)%field))
        dreg = (fmax-fmin)/(size(self%dim(i)%field))
        allocate(rest(size(self%dim(i)%field)-1))
        rest = dims(i)%field(2:)-dims(i)%field(:size(self%dim(i)%field)-1)
        self%regularDim(i) = all(abs(rest - dreg) < abs(eta))
        deallocate(rest)
    end do
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

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that checks the internal consistency of the fields within the
    !> background
    !---------------------------------------------------------------------------
    logical function check(self)
    class(background_class), intent(in) :: self
    integer, allocatable, dimension(:) :: dimSize, fieldShape
    class(*), pointer :: aField
    logical :: equal
    integer :: i
    type(string) :: outext
    check = .true.
    equal = .true.
    if (.not.self%initialized) check = .false.
    allocate(dimSize(size(self%dim)))
    allocate(fieldShape(size(self%dim)))
    do i=1, size(self%dim)
        dimSize(i) = size(self%dim(i)%field)
    end do
    call self%fields%reset()               ! reset list iterator
    do while(self%fields%moreValues())     ! loop while there are values to process
        aField => self%fields%currentValue()
        select type(aField)
        class is (field_class)
            equal = ALL(dimSize.eq.aField%getFieldShape())
            if (.not.equal) check = .false.
            class default
            outext = '[background_class::check] Unexepected type of content, not a scalar Field'
            call Log%put(outext)
            stop
        end select
        call self%fields%next()            ! increment the list iterator
    end do
    call self%fields%reset()               ! reset list iterator

    end function check


    end module background_mod