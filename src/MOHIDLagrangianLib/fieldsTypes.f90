    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : field_types
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines classes for 'fields': 1, 2, 3 and 4D labeled data.
    !> Valid for both scalar and vectorial (real) data. Defines a
    !> generic wrapper for these classes, that abstracts the user from having to
    !> choose their data dimensionality or type to create a field.
    !------------------------------------------------------------------------------

    module fieldTypes_mod

    use common_modules

    implicit none
    private

    !base field class

    type :: field_class      !< a base field class
        type(string) :: name            !< name of the field
        type(string) :: units           !< units of the field, preferably SI please
        integer :: dim                  !< dimensions of the field (1, 2, 3 or 4D)
    contains
    procedure :: setFieldMetadata
    procedure :: print => printField    !< Method that prints the field information
    procedure :: getFieldType           !< Method that returns the field type (scalar or vectorial), in a string
    procedure :: getFieldSlice
    procedure :: getFieldMinBound
    procedure :: getFieldMaxBound
    procedure :: getFieldShape
    end type field_class

    !scalar fields
    !1D (x)
    !2D (t,x)
    !3D (t,x,y)
    !4D (t,x,y,z)

    type, extends(field_class) :: scalar_field_class      !< a scalar field class
    contains
    procedure :: getFieldNearestIndex
    end type scalar_field_class
    
    !Because fortran for some reason does not like me to allocate generic fields in some places
    type, extends(scalar_field_class) :: scalar_dimfield_class      !< a scalar field class
        real(prec), allocatable, dimension(:) :: field1D !< the data on the scalar data 1Dfield
        real(prec), allocatable, dimension(:,:) :: field2D !< the data on the scalar data 2Dfield
        real(prec), allocatable, dimension(:,:,:) :: field3D !< the data on the scalar data 3Dfield
        real(prec), allocatable, dimension(:,:,:,:) :: field4D !< the data on the scalar data 4Dfield
    contains
    procedure :: initScalardim1dField
    procedure :: initScalardim2dField
    procedure :: initScalardim3dField
    procedure :: initScalardim4dField
    generic   :: initialize => initScalardim1dField, initScalardim2dField, initScalardim3dField, initScalardim4dField
    procedure :: finalize => cleanScalardimField
    end type scalar_dimfield_class

    type, extends(scalar_field_class) :: scalar1d_field_class      !< a 1D scalar field class
        real(prec), allocatable, dimension(:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar1dField
    !procedure :: getFieldNearestIndex
    procedure :: finalize => cleanScalar1dField
    end type scalar1d_field_class

    type, extends(scalar_field_class) :: scalar2d_field_class      !< a 2D scalar field class
        real(prec), allocatable, dimension(:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar2dField
    !procedure :: getFieldNearestIndex2D
    procedure :: finalize => cleanScalar2dField
    end type scalar2d_field_class

    type, extends(scalar_field_class) :: scalar3d_field_class      !< a 3D scalar field class
        real(prec), allocatable, dimension(:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar3dField
    procedure :: finalize => cleanScalar3dField
    end type scalar3d_field_class
    
    type, extends(scalar_field_class) :: scalar3d_int_field_class      !< a 3D scalar field class
        integer, allocatable, dimension(:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initIntScalar3dField
    procedure :: finalize => cleanIntScalar3dField
    end type scalar3d_int_field_class

    type, extends(scalar_field_class) :: scalar4d_field_class      !< a 4D scalar field class
        real(prec), allocatable, dimension(:,:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar4dField
    procedure :: finalize => cleanScalar4dField
    end type scalar4d_field_class
    
    type, extends(scalar_field_class) :: scalar4d_int_field_class      !< a 4D scalar field class
        integer, allocatable, dimension(:,:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initIntScalar4dField
    procedure :: finalize => cleanIntScalar4dField
    end type scalar4d_int_field_class

    !Vectorial fields
    !2D (t,vx)
    !3D (t,vx,vy)
    !4D (t,vx,vy,vz)

    type, extends(field_class) :: vectorial_field_class      !< a vectorial field class
    contains
    end type vectorial_field_class

    type, extends(vectorial_field_class) :: vectorial2d_field_class      !< a 2D vectorial field class
        type(vector), allocatable, dimension(:,:) :: field !< the data on the 2D vectorial data field
    contains
    procedure :: initialize => initVectorial2dField
    end type vectorial2d_field_class

    type, extends(vectorial_field_class) :: vectorial3d_field_class      !< a 3D vectorial field class
        type(vector), allocatable, dimension(:,:,:) :: field !< the data on the 3D vectorial data field
    contains
    procedure :: initialize => initVectorial3dField
    end type vectorial3d_field_class

    type, extends(vectorial_field_class) :: vectorial4d_field_class      !< a 4D vectorial field class
        type(vector), allocatable, dimension(:,:,:,:) :: field !< the data on the 4D vectorial data field
    contains
    procedure :: initialize => initVectorial4dField
    end type vectorial4d_field_class

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type, extends(field_class) :: generic_field_class !< generic field class. This works as a wrapper for a generic initialization routine.
        type(scalar1d_field_class) :: scalar1d       !< 1D scalar field
        type(scalar2d_field_class) :: scalar2d       !< 2D scalar field
        type(scalar3d_field_class) :: scalar3d       !< 3D scalar field
        type(scalar4d_field_class) :: scalar4d       !< 4D scalar field
        type(vectorial2d_field_class) :: vectorial2d !< 2D vectorial field
        type(vectorial3d_field_class) :: vectorial3d !< 3D vectorial field
        type(vectorial4d_field_class) :: vectorial4d !< 4D vectorial field
        type(scalar3d_int_field_class) :: intScalar3D !< 4D vectorial field
        type(scalar4d_int_field_class) :: intScalar4D !< 4D vectorial field
    contains
    procedure, private :: test
    procedure :: initS1D, initS2D, initS3D, initS4D
    procedure :: initV2D, initV3D, initV4D
    procedure :: initIS3D, initIS4D
    generic   :: initialize => initS1D, initS2D, initS3D, initS4D, initV2D, initV3D, initV4D, initIS3D, initIS4D
    procedure :: replaceMetaData
    procedure :: compare
    procedure :: concatenate
    procedure :: getGFieldType
    procedure :: getGField
    procedure :: finalize => cleanField
    final :: delGfield
    procedure :: print => printGenericField
    end type generic_field_class

    !Public access vars
    public :: field_class, generic_field_class
    public :: scalar_field_class, scalar_dimfield_class, scalar1d_field_class, scalar2d_field_class, scalar3d_field_class, scalar4d_field_class
    public :: vectorial_field_class, vectorial2d_field_class, vectorial3d_field_class, vectorial4d_field_class
    public :: initIntScalar3dField, initIntScalar4dField, cleanIntScalar3dField, cleanIntScalar4dField

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method concatenates two fields by their last dimension
    !> @param[in] self, gfield, usedPosi
    !---------------------------------------------------------------------------
    subroutine concatenate(self, gfield, usedPosi)
    class(generic_field_class), intent(inout) :: self
    type(generic_field_class), intent(in) :: gfield
    logical, dimension(:), optional, intent(in) :: usedPosi 
    logical, allocatable, dimension(:) :: usedPos
    integer :: fDim, nUsedPos, i, j
    type(string) :: fType
    real(prec), allocatable, dimension(:) :: field1d
    real(prec), allocatable, dimension(:,:) :: field2d
    real(prec), allocatable, dimension(:,:,:) :: field3d
    real(prec), allocatable, dimension(:,:,:,:) :: field4d
    
    !write(*,*)"Entrada concatenate" 

    fDim = self%dim
    fType = self%getGFieldType()
    !write(*,*)"Got fieldtype" 
    if (gfield%dim /= fDim) return
    if (gfield%getGFieldType() /= fType) return
    
    if (present(usedPosi)) then
        allocate(usedPos, source = usedPosi)
        nUsedPos = count(usedPos)
    end if 
    !write(*,*)"fType : ", fType 
    if (fType == 'Scalar') then
        if (fDim == 1) then
            !write(*,*)"fDim = 1" 
            if (.not.allocated(usedPos)) then
                nUsedPos = size(gfield%scalar1d%field)
                allocate(usedPos(nUsedPos))
                usedPos = .true.
            end if
            allocate(field1d(size(self%scalar1d%field) + nUsedPos))
            field1d(1:size(self%scalar1d%field)) = self%scalar1d%field
            i= size(self%scalar1d%field)+1
            do j= 1, size(usedPos)
                if (usedPos(j)) then
                    field1d(i) = gfield%scalar1d%field(j)
                    i= i+1
                end if
            end do
            deallocate(self%scalar1d%field)
            allocate(self%scalar1d%field(size(field1d)))
            self%scalar1d%field = field1d
        end if
        if (fDim == 2) then
            !write(*,*)"fDim = 2"
            if (.not.allocated(usedPos)) then
                nUsedPos = size(gfield%scalar2d%field,2)
                allocate(usedPos(nUsedPos))
                usedPos = .true.
            end if
            allocate(field2d(size(self%scalar2d%field,1),size(self%scalar2d%field,2) + nUsedPos))
            field2d(:,1:size(self%scalar2d%field,2)) = self%scalar2d%field
            i= size(self%scalar2d%field, 2) +1
            do j= 1, size(usedPos)
                if (usedPos(j)) then
                    field2d(:, i) = gfield%scalar2d%field(:, j)
                    i= i+1
                end if
            end do            
            deallocate(self%scalar2d%field)            
            allocate(self%scalar2d%field(size(field2d,1), size(field2d,2)))
            self%scalar2d%field = field2d
        end if
        if (fDim == 3) then
            !write(*,*)"fDim = 3"
            if (.not.allocated(usedPos)) then
                nUsedPos = size(gfield%scalar3d%field,3)
                allocate(usedPos(nUsedPos))
                usedPos = .true.
            end if
            allocate(field3d(size(self%scalar3d%field,1),size(self%scalar3d%field,2),size(self%scalar3d%field,3) + nUsedPos))
            field3d(:,:,1:size(self%scalar3d%field,3)) = self%scalar3d%field
            i= size(self%scalar3d%field, 3) +1
            do j= 1, size(usedPos)
                if (usedPos(j)) then
                    field3d(:,:, i) = gfield%scalar3d%field(:,:, j)
                    i= i+1
                end if
            end do
            deallocate(self%scalar3d%field)
            allocate(self%scalar3d%field(size(field3d,1), size(field3d,2), size(field3d,3)))
            self%scalar3d%field = field3d
        end if
        if (fDim == 4) then
            !write(*,*)"fDim = 4"
            if (.not.allocated(usedPos)) then
                nUsedPos = size(gfield%scalar4d%field,4)
                !write(*,*)"allocating used Pos, nUsedPos = ", nUsedPos
                allocate(usedPos(nUsedPos))
                usedPos = .true.
            end if
            !write(*,*)"allocating field 4D"
            allocate(field4d(size(self%scalar4d%field,1),size(self%scalar4d%field,2),size(self%scalar4d%field,3),size(self%scalar4d%field,4) + nUsedPos))
            !write(*,*)"allocated field 4D"
            field4d(:,:,:,1:size(self%scalar4d%field,4)) = self%scalar4d%field
            !write(*,*)"made equal to self%scalar4d%field. Size field4d = ", size(field4d, 4)
            i = size(self%scalar4d%field, 4) + 1
            !write(*,*)"i = ", i
            !write(*,*)"size(usedPos) = ", size(usedPos)
            !write(*,*)"size field (1) ", size(gfield%scalar4d%field, 1)
            !write(*,*)"size field (2) ", size(gfield%scalar4d%field, 2)
            !write(*,*)"size field (3) ", size(gfield%scalar4d%field, 3)
            !write(*,*)"size field (4) ", size(gfield%scalar4d%field, 4)
            !write(*,*)"size field4d(1) ", size(field4d, 1)
            !write(*,*)"size field4d(2) ", size(field4d, 2)
            !write(*,*)"size field4d(3) ", size(field4d, 3)
            !write(*,*)"size field4d(4) ", size(field4d, 4)
            do j= 1, size(usedPos)
                !write(*,*)"j = ", j
                if (usedPos(j)) then
                    !write(*,*)"used position . i = , j = ", i, j
                    field4d(:,:,:, i) = gfield%scalar4d%field(:,:,:, j)
                    i= i+1
                end if
            end do
            !write(*,*)"cheguei aqui"
            deallocate(self%scalar4d%field)
            allocate(self%scalar4d%field(size(field4d,1), size(field4d,2), size(field4d,3), size(field4d,4)))
            self%scalar4d%field = field4d
        end if
    else if (fType == 'Vectorial') then
        return
    end if

    end subroutine concatenate
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method concatenates two fields by their last dimension
    !> @param[in] self, gfield, usedPosi
    !---------------------------------------------------------------------------
    subroutine concatenate_Int(self, gfield, usedPosi)
    class(generic_field_class), intent(inout) :: self
    type(generic_field_class), intent(in) :: gfield
    logical, dimension(:), optional, intent(in) :: usedPosi 
    logical, allocatable, dimension(:) :: usedPos
    integer :: fDim, nUsedPos, i, j
    !integer, allocatable, dimension(:) :: field1d
    !integer, allocatable, dimension(:,:) :: field2d
    integer, allocatable, dimension(:,:,:) :: field3d
    integer, allocatable, dimension(:,:,:,:) :: field4d
    
    !write(*,*)"Entrada concatenate" 

    fDim = self%dim
    !write(*,*)"Got fieldtype" 
    if (gfield%dim /= fDim) return
    
    if (present(usedPosi)) then
        allocate(usedPos, source = usedPosi)
        nUsedPos = count(usedPos)
    end if 
    !write(*,*)"fType : ", fType 
    if (fDim == 3) then
        !write(*,*)"fDim = 3"
        if (.not.allocated(usedPos)) then
            nUsedPos = size(gfield%intScalar3D%field,3)
            allocate(usedPos(nUsedPos))
            usedPos = .true.
        end if
        allocate(field3d(size(self%intScalar3D%field,1),size(self%intScalar3D%field,2),size(self%intScalar3D%field,3) + nUsedPos))
        field3d(:,:,1:size(self%intScalar3D%field,3)) = self%intScalar3D%field
        i= size(self%intScalar3D%field, 3) +1
        do j= 1, size(usedPos)
            if (usedPos(j)) then
                field3d(:,:, i) = gfield%intScalar3D%field(:,:, j)
                i= i+1
            end if
        end do
        deallocate(self%intScalar3D%field)
        allocate(self%intScalar3D%field(size(field3d,1), size(field3d,2), size(field3d,3)))
        self%intScalar3D%field = field3d
    end if
    if (fDim == 4) then
        !write(*,*)"fDim = 4"
        if (.not.allocated(usedPos)) then
            nUsedPos = size(gfield%intScalar4D%field,4)
            !write(*,*)"allocating used Pos, nUsedPos = ", nUsedPos
            allocate(usedPos(nUsedPos))
            usedPos = .true.
        end if
        !write(*,*)"allocating field 4D"
        allocate(field4d(size(self%intScalar4D%field,1),size(self%intScalar4D%field,2),size(self%intScalar4D%field,3),size(self%intScalar4D%field,4) + nUsedPos))
        !write(*,*)"allocated field 4D"
        field4d(:,:,:,1:size(self%intScalar4D%field,4)) = self%intScalar4D%field
        !write(*,*)"made equal to self%intScalar3D%field. Size field4d = ", size(field4d, 4)
        i = size(self%intScalar4D%field, 4) + 1
        !write(*,*)"i = ", i
        !write(*,*)"size(usedPos) = ", size(usedPos)
        !write(*,*)"size field (1) ", size(gfield%intScalar3D%field, 1)
        !write(*,*)"size field (2) ", size(gfield%intScalar3D%field, 2)
        !write(*,*)"size field (3) ", size(gfield%intScalar3D%field, 3)
        !write(*,*)"size field (4) ", size(gfield%intScalar3D%field, 4)
        !write(*,*)"size field4d(1) ", size(field4d, 1)
        !write(*,*)"size field4d(2) ", size(field4d, 2)
        !write(*,*)"size field4d(3) ", size(field4d, 3)
        !write(*,*)"size field4d(4) ", size(field4d, 4)
        do j= 1, size(usedPos)
            !write(*,*)"j = ", j
            if (usedPos(j)) then
                !write(*,*)"used position . i = , j = ", i, j
                field4d(:,:,:, i) = gfield%intScalar4D%field(:,:,:, j)
                i= i+1
            end if
        end do
        !write(*,*)"cheguei aqui"
        deallocate(self%intScalar4D%field)
        allocate(self%intScalar4D%field(size(field4d,1), size(field4d,2), size(field4d,3), size(field4d,4)))
        self%intScalar4D%field = field4d
    end if

    end subroutine concatenate_Int

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that compares the metadata of a scalar 1D field to the object field
    !> @param[in] self, gfield
    !---------------------------------------------------------------------------
    logical function compare(self, gfield) result(comp)
    class(generic_field_class), intent(inout) :: self
    class(generic_field_class), intent(in) :: gfield
    comp = .true.
    if (self%name /= gfield%name) comp = .false.
    if (self%units /= gfield%units) comp = .false. !might be asking too much...
    if (self%dim /= gfield%dim) comp = .false.
    if (allocated(self%scalar1d%field) .and. .not. allocated(gfield%scalar1d%field)) comp = .false.
    if (allocated(self%scalar2d%field) .and. .not. allocated(gfield%scalar2d%field)) comp = .false.
    if (allocated(self%scalar3d%field) .and. .not. allocated(gfield%scalar3d%field)) comp = .false.
    if (allocated(self%scalar4d%field) .and. .not. allocated(gfield%scalar4d%field)) comp = .false.
    if (allocated(self%vectorial2d%field) .and. .not. allocated(gfield%vectorial2d%field)) comp = .false.
    if (allocated(self%vectorial3d%field) .and. .not. allocated(gfield%vectorial3d%field)) comp = .false.
    if (allocated(self%vectorial4d%field) .and. .not. allocated(gfield%vectorial4d%field)) comp = .false.
    end function compare
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> replaces metadata on a generic field
    !> @param[in] self, name, units
    !---------------------------------------------------------------------------
    subroutine replaceMetaData(self, name, units)
    class(generic_field_class), intent(inout) :: self
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    self%name = name
    self%units = units
    end subroutine replaceMetaData

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> cleans a generic field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanField(self)
    class(generic_field_class), intent(inout) :: self
    self%name = ''
    self%units = ''
    self%dim = MV
    if (allocated(self%scalar1d%field)) deallocate(self%scalar1d%field)
    if (allocated(self%scalar2d%field)) deallocate(self%scalar2d%field)
    if (allocated(self%scalar3d%field)) deallocate(self%scalar3d%field)
    if (allocated(self%scalar4d%field)) deallocate(self%scalar4d%field)
    if (allocated(self%vectorial2d%field)) deallocate(self%vectorial2d%field)
    if (allocated(self%vectorial3d%field)) deallocate(self%vectorial3d%field)
    if (allocated(self%vectorial4d%field)) deallocate(self%vectorial4d%field)
    end subroutine cleanField

    subroutine delGfield(self)
    type(generic_field_class), intent(inout) :: self
    self%name = ''
    self%units = ''
    self%dim = MV
    if (allocated(self%scalar1d%field)) deallocate(self%scalar1d%field)
    if (allocated(self%scalar2d%field)) deallocate(self%scalar2d%field)
    if (allocated(self%scalar3d%field)) deallocate(self%scalar3d%field)
    if (allocated(self%scalar4d%field)) deallocate(self%scalar4d%field)
    if (allocated(self%vectorial2d%field)) deallocate(self%vectorial2d%field)
    if (allocated(self%vectorial3d%field)) deallocate(self%vectorial3d%field)
    if (allocated(self%vectorial4d%field)) deallocate(self%vectorial4d%field)
    end subroutine delGfield

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 1D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initS1D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%scalar1d%field)) then
        outext = '[generic_field_class::initialize]: scalar 1D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%scalar1d%initialize(name, units, 1, field)
        call self%setFieldMetadata(name, units, 1)
    end if
    end subroutine initS1D

    !---------------------------------------------------------------------------
    !!> @author Ricardo Birjukovs Canelas - MARETEC
    !!> @brief
    !!> Method that returns the field's nearest value index (scalar)
    !!> @param[in] self, value
    !!---------------------------------------------------------------------------
    !integer function getFieldNearestIndex(self, value)
    !class(scalar1d_field_class), intent(in) :: self
    !real(prec), intent(in) :: value
    !real(prec), allocatable, dimension(:) :: comp
    !allocate(comp(size(self%field)))
    !comp = value
    !getFieldNearestIndex = minloc(abs(comp - self%field), DIM=1)
    !end function getFieldNearestIndex
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the field's nearest value index (scalar)
    !> @param[in] self, value
    !---------------------------------------------------------------------------
    integer function getFieldNearestIndex(self, value, dimID)
    class(scalar_field_class), intent(in) :: self
    real(prec), intent(in) :: value
    integer, intent(in), optional :: dimID !lat = 1 or lon = 2
    real(prec), allocatable, dimension(:) :: comp1d
    real(prec), allocatable, dimension(:,:) :: comp2d
    real(prec), allocatable, dimension(:,:,:,:) :: comp4d
    integer, allocatable, dimension(:) :: aux
    integer, allocatable, dimension(:,:) :: aux_2D
    type(string) :: outext
    !Begin-----------------------------------------------
    select type(self)
    class is (scalar1d_field_class)
        allocate(comp1d(size(self%field)))
        comp1d = value
        getFieldNearestIndex = minloc(abs(comp1d - self%field), DIM=1)
    class is (scalar_dimfield_class)
        if (allocated(self%field1D)) then 
            allocate(comp1d(size(self%field1D)))
            comp1d = value
            getFieldNearestIndex = minloc(abs(comp1d - self%field1D), DIM=1)
        elseif (allocated(self%field2D)) then
            if (.not. present(dimID)) then
                outext = '[field_class::getFieldNearestIndex]: for a 2D field, dimID variable must be provided'
                call Log%put(outext)
                stop
            endif
            allocate(comp2d(size(self%field2D,1), size(self%field2D,2)))
            !minloc will give an array the size of the perpendicular direction of dimID
            if (dimID == 1) allocate(aux(size(self%field2D,2)))
            if (dimID == 2) allocate(aux(size(self%field2D,1)))
            comp2d = value
            aux = minloc(abs(comp2d - self%field2D), DIM=dimID) !output is a 1D array
            getFieldNearestIndex = minval(aux) !get the smallest row or column ID
        elseif (allocated(self%field4D)) then
            if (.not. present(dimID)) then
                outext = '[field_class::getFieldNearestIndex]: for a 4D field, dimID variable must be provided'
                call Log%put(outext)
                stop
            endif
            allocate(comp4d(size(self%field4D,1), size(self%field4D,2), size(self%field4D,3), size(self%field4D,4)))
            !minloc will give an array the size of the perpendicular direction of dimID
            allocate(aux(4))
            comp4d = value
            aux = minloc(abs(comp4d - self%field4D))
            getFieldNearestIndex = max(aux(3) - 1, 1) !get the smallest vertical layer ID
        else
            outext = '[field_class::getFieldNearestIndex]: scalar field must be 1D or 2D or 4D'
            call Log%put(outext)
            stop
        endif
    end select
    end function getFieldNearestIndex
    
    !!> @author Joao Sobrinho - Colab Atlantic
    !!> @brief
    !!> Method that returns the 2D field's nearest value index (scalar)
    !!> @param[in] self, value
    !!---------------------------------------------------------------------------
    !integer function getFieldNearestIndex2D(self, value, dimID)
    !class(scalar2d_field_class), intent(in) :: self
    !real(prec), intent(in) :: value
    !real(prec), allocatable, dimension(:,:) :: comp
    !integer, intent(in) :: dimID !lat = 1 or lon = 2
    !integer, dimension(size(self%field,dimID)) :: aux
    !type(string) :: outext
    !!begin--------------------------------------------------------------
    !allocate(comp(size(self%field,1), size(self%field,2)))
    !comp = value
    !aux = minloc(abs(comp - self%field), DIM=dimID) !output is a 1D array
    !getFieldNearestIndex2D = minval(aux) !get the smallest row or column ID
    !end function getFieldNearestIndex2D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 2D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initS2D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%scalar2d%field)) then
        outext = '[generic_field_class::initialize]: scalar 2D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%scalar2d%initialize(name, units, 2, field)
        call self%setFieldMetadata(name, units, 2)
    end if
    end subroutine initS2D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 3D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initS3D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%scalar3d%field)) then
        outext = '[generic_field_class::initialize]: scalar 3D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%scalar3d%initialize(name, units, 3, field)
        call self%setFieldMetadata(name, units, 3)
    end if
    end subroutine initS3D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 4D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initS4D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:,:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%scalar4d%field)) then
        outext = '[generic_field_class::initialize]: scalar 4D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%scalar4d%initialize(name, units, 4, field)
        call self%setFieldMetadata(name, units, 4)
    end if
    end subroutine initS4D
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that allocates and initializes a scalar integer 3D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initIS3D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    integer, intent(in), dimension(:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%intScalar3D%field)) then
        outext = '[generic_field_class::initialize]: scalar integer 3D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%intScalar3D%initialize(name, units, 3, field)
        call self%setFieldMetadata(name, units, 3)
    end if
    end subroutine initIS3D

        !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that allocates and initializes a scalar integer 3D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initIS4D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    integer, intent(in), dimension(:,:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%intScalar4D%field)) then
        outext = '[generic_field_class::initialize]: scalar integer 4D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%intScalar4D%initialize(name, units, 4, field)
        call self%setFieldMetadata(name, units, 4)
    end if
    end subroutine initIS4D

    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a vectorial 2D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initV2D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    type(vector), intent(in), dimension(:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%vectorial2d%field)) then
        outext = '[generic_field_class::initialize]: vectorial 2D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%vectorial2d%initialize(name, units, 2, field)
        call self%setFieldMetadata(name, units, 2)
    end if
    end subroutine initV2D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a vectorial 3D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initV3D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    type(vector), intent(in), dimension(:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%vectorial3d%field)) then
        outext = '[generic_field_class::initialize]: vectorial 3D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%vectorial3d%initialize(name, units, 3, field)
        call self%setFieldMetadata(name, units, 3)
    end if
    end subroutine initV3D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a vectorial 4D field in a generic field
    !> @param[in] self, name, units, field
    !---------------------------------------------------------------------------
    subroutine initV4D(self, name, units, field)
    class(generic_field_class), intent(inout) :: self
    type(vector), intent(in), dimension(:,:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    type(string) :: outext
    if (allocated(self%vectorial4d%field)) then
        outext = '[generic_field_class::initialize]: vectorial 4D field already allocated'
        call Log%put(outext)
        stop
    else
        call self%vectorial4d%initialize(name, units, 4, field)
        call self%setFieldMetadata(name, units, 4)
    end if
    end subroutine initV4D
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - colab atlantic
    !> @brief
    !> Method that initializes a scalar dim 1D field
    !> @param[in] self, name, units, dim, field1D
    !---------------------------------------------------------------------------
    subroutine initScalardim1dField(self, name, units, dim, field1D)
    class(scalar_dimfield_class), intent(inout) :: self
    real(prec), dimension(:), intent(in) :: field1D
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field1D, source = field1D)
    end subroutine initScalardim1dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - colab atlantic
    !> @brief
    !> Method that initializes a scalar dim 2D field
    !> @param[in] self, name, units, dim, field2D
    !---------------------------------------------------------------------------
    subroutine initScalardim2dField(self, name, units, dim, field2D)
    class(scalar_dimfield_class), intent(inout) :: self
    real(prec), dimension(:,:), intent(in) :: field2D
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field2D, source = field2D)
    end subroutine initScalardim2dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 1D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initScalar1dField(self, name, units, dim, field)
    class(scalar1d_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initScalar1dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that initializes a scalar dim 3D field
    !> @param[in] self, name, units, dim, field3D
    !---------------------------------------------------------------------------
    subroutine initScalardim3dField(self, name, units, dim, field3D)
    class(scalar_dimfield_class), intent(inout) :: self
    real(prec), dimension(:,:,:), intent(in) :: field3D
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field3D, source = field3D)
    end subroutine initScalardim3dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that initializes a scalar dim 4D field
    !> @param[in] self, name, units, dim, field3D
    !---------------------------------------------------------------------------
    subroutine initScalardim4dField(self, name, units, dim, field4D)
    class(scalar_dimfield_class), intent(inout) :: self
    real(prec), dimension(:,:,:,:), intent(in) :: field4D
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field4D, source = field4D)
    end subroutine initScalardim4dField

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Method that deallocates a scalar dim field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanScalardimField(self)
    class(scalar_dimfield_class), intent(out) :: self
    if (allocated(self%field1D)) deallocate(self%field1D)
    if (allocated(self%field2D)) deallocate(self%field2D)
    end subroutine cleanScalardimField
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 1D field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanScalar1dField(self)
    class(scalar1d_field_class), intent(out) :: self
    if (allocated(self%field)) deallocate(self%field)
    end subroutine cleanScalar1dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 2D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initScalar2dField(self, name, units, dim, field)
    class(scalar2d_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    write(*,*) "Sizes = ", SIZE(self%field, 1), SIZE(self%field, 2)
    end subroutine initScalar2dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that deallocates a scalar 2D field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanScalar2dField(self)
    class(scalar2d_field_class), intent(out) :: self
    if (allocated(self%field)) deallocate(self%field)
    end subroutine cleanScalar2dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 3D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initScalar3dField(self, name, units, dim, field)
    class(scalar3d_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initScalar3dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that initializes a scalar integer 3D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initIntScalar3dField(self, name, units, dim, field)
    class(scalar3d_int_field_class), intent(inout) :: self
    integer, intent(in), dimension(:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initIntScalar3dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that deallocates a scalar 3D field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanScalar3dField(self)
    class(scalar3d_field_class), intent(out) :: self
    if (allocated(self%field)) deallocate(self%field)
    end subroutine cleanScalar3dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that deallocates a scalar int 3D field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanIntScalar3dField(self)
    class(scalar3d_int_field_class), intent(out) :: self
    if (allocated(self%field)) deallocate(self%field)
    end subroutine cleanIntScalar3dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 4D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initScalar4dField(self, name, units, dim, field)
    class(scalar4d_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:,:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initScalar4dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that initializes a scalar integer 4D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initIntScalar4dField(self, name, units, dim, field)
    class(scalar4d_int_field_class), intent(inout) :: self
    integer, intent(in), dimension(:,:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initIntScalar4dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that deallocates a scalar 4D field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanScalar4dField(self)
    class(scalar4d_field_class), intent(out) :: self
    if (allocated(self%field)) deallocate(self%field)
    end subroutine cleanScalar4dField
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Method that deallocates a scalar integer 4D field
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine cleanIntScalar4dField(self)
    class(scalar4d_int_field_class), intent(out) :: self
    if (allocated(self%field)) deallocate(self%field)
    end subroutine cleanIntScalar4dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a vectorial 2D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initVectorial2dField(self, name, units, dim, field)
    class(vectorial2d_field_class), intent(inout) :: self
    type(vector), intent(in), dimension(:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initVectorial2dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a vectorial 3D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initVectorial3dField(self, name, units, dim, field)
    class(vectorial3d_field_class), intent(inout) :: self
    type(vector), intent(in), dimension(:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initVectorial3dField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a vectorial 4D field
    !> @param[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initVectorial4dField(self, name, units, dim, field)
    class(vectorial4d_field_class), intent(inout) :: self
    type(vector), intent(in), dimension(:,:,:,:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    call self%setFieldMetadata(name, units, dim)
    allocate(self%field, source = field)
    end subroutine initVectorial4dField


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a base field object by filling metadata
    !> @param[in] self, name, units, dim
    !---------------------------------------------------------------------------
    subroutine setFieldMetadata(self, name, units, dim)
    class(field_class), intent(inout) :: self
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    self%name = name
    self%units = units
    self%dim = dim
    end subroutine setFieldMetadata

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the generic field information
    !---------------------------------------------------------------------------
    subroutine printGenericField(self)
    class(generic_field_class), intent(in) :: self
    if (allocated(self%scalar1d%field)) call self%scalar1d%print()
    if (allocated(self%scalar2d%field)) call self%scalar2d%print()
    if (allocated(self%scalar3d%field)) call self%scalar3d%print()
    if (allocated(self%scalar4d%field)) call self%scalar4d%print()
    if (allocated(self%vectorial2d%field)) call self%vectorial2d%print()
    if (allocated(self%vectorial3d%field)) call self%vectorial3d%print()
    if (allocated(self%vectorial4d%field)) call self%vectorial4d%print()
    end subroutine printGenericField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the field type (scalar or vectorial), in a string,
    !> of a generic field
    !---------------------------------------------------------------------------
    type(string) function getGFieldType(self)
    class(generic_field_class), intent(in) :: self
    if (allocated(self%scalar1d%field)) getGFieldType = self%scalar1d%getFieldType()
    if (allocated(self%scalar2d%field)) getGFieldType = self%scalar2d%getFieldType()
    if (allocated(self%scalar3d%field)) getGFieldType = self%scalar3d%getFieldType()
    if (allocated(self%scalar4d%field)) getGFieldType = self%scalar4d%getFieldType()
    if (allocated(self%vectorial2d%field)) getGFieldType = self%vectorial2d%getFieldType()
    if (allocated(self%vectorial3d%field)) getGFieldType = self%vectorial3d%getFieldType()
    if (allocated(self%vectorial4d%field)) getGFieldType = self%vectorial4d%getFieldType()
    end function getGFieldType

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> A class 'unit' test for the generic_field_class
    !---------------------------------------------------------------------------
    subroutine test(self)
    class(generic_field_class), intent(inout) :: self
    type(generic_field_class) :: gfield1, gfield2, gfield3
    real(prec), allocatable, dimension(:) :: field1
    real(prec), allocatable, dimension(:,:) :: field2
    type(vector), allocatable, dimension(:,:,:) :: field3
    type(string) :: name1, name2, name3
    type(string) :: units1, units2, units3
    allocate(field1(50))
    allocate(field2(20,60))
    allocate(field3(2,3,4))
    name1 = 'testfield1d'
    name2 = 'testfield2d'
    name3 = 'testfield3d'
    units1 = 'm/s'
    units2 = 'km'
    units3 = 'ms-1'
    call gfield1%initialize(name1, units1, field1)
    call gfield2%initialize(name2, units2, field2)
    call gfield3%initialize(name3, units3, field3)
    call gfield1%print()
    call gfield2%print()
    call gfield3%print()
    end subroutine test

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the field information
    !---------------------------------------------------------------------------
    subroutine printField(self)
    class(field_class), intent(in) :: self
    type(string) :: outext, t(5)
    t(1) = self%dim
    t(2) = self%getFieldType()
    t(3) = self%getFieldMinBound()
    t(4) = self%getFieldMaxBound()
    !t(5) = self%getFieldShape()
    outext = t(2)//' field['//self%name//'] has dimensionality '//t(1)//' and is in '//self%units//' - ['//t(3)//';'//t(4)//']'!, ('//t(5)//')'
    call Log%put(outext,.false.)
    end subroutine printField

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the field type (scalar or vectorial), in a string
    !---------------------------------------------------------------------------
    function getFieldType(self)
    class(field_class), intent(in) :: self
    type(string) :: getFieldType
    type(string) :: outext
    select type(self)
    class is (scalar_field_class)
        getFieldType = 'Scalar'
    class is (vectorial_field_class)
        getFieldType = 'Vectorial'
        class default
        outext = '[field_class::getFieldType]: Unexepected type of content, not a scalar or vectorial Field'
        call Log%put(outext)
        stop
    end select
    end function getFieldType

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a slice of a field stored on a generic field
    !---------------------------------------------------------------------------
    function getFieldSlice(self, llbound, uubound)
    class(field_class), intent(in) :: self
    integer, dimension(:), intent(in) :: llbound, uubound
    type(generic_field_class) :: getFieldSlice
    logical :: sliced
    type(string) :: outext
    sliced = .false.
    select type(self)
    class is (scalar_field_class)
    class is (scalar1d_field_class)
        if (size(llbound) == 1) then
            call getFieldSlice%initialize(self%name, self%units, self%field(llbound(1):uubound(1)))
            sliced = .true.
        end if
    class is (scalar2d_field_class)
        if (size(llbound) == 2) then
            call getFieldSlice%initialize(self%name, self%units, self%field(llbound(1):uubound(1),llbound(2):uubound(2)))
            sliced = .true.
        end if
    class is (scalar3d_field_class)
        if (size(llbound) == 3) then
            !write(*,*)"entrei no slice 3D. "
            !write(*,*)"llbound(1):uubound(1) = ", llbound(1), uubound(1)
            !write(*,*)"llbound(2):uubound(2) = ", llbound(2), uubound(2)
            !write(*,*)"llbound(3):uubound(3) = ", llbound(3), uubound(3)
            !write(*,*)"size field 1 = ", size(self%field,1)
            !write(*,*)"size field 2 = ", size(self%field,2)
            !write(*,*)"size field 3 = ", size(self%field,3)
            call getFieldSlice%initialize(self%name, self%units, self%field(llbound(1):uubound(1), llbound(2):uubound(2), llbound(3):uubound(3)))
            sliced = .true.
        end if
    class is (scalar4d_field_class)
        if (size(llbound) == 4) then
            call getFieldSlice%initialize(self%name, self%units, self%field(llbound(1):uubound(1), llbound(2):uubound(2), llbound(3):uubound(3), llbound(4):uubound(4)))
            sliced = .true.
        end if
        class default
        outext = '[field_class::getFieldSlice]: Unexepected type of content, not a scalar Field'
        call Log%put(outext)
        stop
    end select
    if (.not.sliced) then
        outext = '[field_class::getFieldSlice]: type of Field and dimensions do not match'
        call Log%put(outext)
        stop
    end if
    end function getFieldSlice

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a slice of a field stored on a generic field
    !---------------------------------------------------------------------------
    function getFieldShape(self)
    class(field_class), intent(in) :: self
    integer, allocatable, dimension(:) :: getFieldShape
    type(string) :: outext

    select type(self)
    class is (scalar_field_class)
    class is (scalar1d_field_class)
        allocate(getFieldShape(1))
        getFieldShape(1) = size(self%field)
    class is (scalar2d_field_class)
        allocate(getFieldShape(2))
        getFieldShape(1) = size(self%field,1)
        getFieldShape(2) = size(self%field,2)
    class is (scalar3d_field_class)
        allocate(getFieldShape(3))
        getFieldShape(1) = size(self%field,1)
        getFieldShape(2) = size(self%field,2)
        getFieldShape(3) = size(self%field,3)
    class is (scalar4d_field_class)
        allocate(getFieldShape(4))
        getFieldShape(1) = size(self%field,1)
        getFieldShape(2) = size(self%field,2)
        getFieldShape(3) = size(self%field,3)
        getFieldShape(4) = size(self%field,4)
        class default
        outext = '[field_class::getFieldShape]: Unexepected type of content, not a scalar Field'
        call Log%put(outext)
        stop
    end select
    end function getFieldShape

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Updated in Dec 2022 by Joao Sbrinho - Colab Atlantic
    !> Method that returns the field's maximum value (scalar)
    !---------------------------------------------------------------------------
    real(prec) function getFieldMaxBound(self, arrayDim)
    class(field_class), intent(in) :: self
    type(string) :: outext
    integer, intent(in), optional :: arrayDim !Dimension (columns, rows, depth or time)
    !Begin-----------------------------------------------
    select type(self)
    class is (scalar_field_class)
    class is (scalar_dimfield_class)
        if (allocated(self%field1D)) then
            getFieldMaxBound = maxval(self%field1D)
        elseif (allocated(self%field2D)) then
            getFieldMaxBound = maxval(self%field2D)
        end if
    class is (scalar1d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (scalar2d_field_class)
        if (present(arrayDim)) then
            if (arrayDim == 1) getFieldMaxBound = minval(self%field(:,1)) !Lon
            if (arrayDim == 2) getFieldMaxBound = minval(self%field(1,:)) !Lat
        else
            getFieldMaxBound = maxval(self%field)
        end if
        getFieldMaxBound = maxval(self%field)
    class is (scalar3d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (scalar4d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (vectorial_field_class)
        getFieldMaxBound = MV
    class is (generic_field_class)
        if (allocated(self%scalar1d%field)) getFieldMaxBound = maxval(self%scalar1d%field)
        if (allocated(self%scalar2d%field)) getFieldMaxBound = maxval(self%scalar2d%field)
        if (allocated(self%scalar3d%field)) getFieldMaxBound = maxval(self%scalar3d%field)
        if (allocated(self%scalar4d%field)) getFieldMaxBound = maxval(self%scalar4d%field)
        if (allocated(self%vectorial2d%field)) getFieldMaxBound = MV
        if (allocated(self%vectorial3d%field)) getFieldMaxBound = MV
        if (allocated(self%vectorial4d%field)) getFieldMaxBound = MV
        class default
        outext = '[field_class::getFieldMaxBound]: Unexepected type of content, not a scalar or vectorial Field'
        call Log%put(outext)
        stop
    end select
    end function getFieldMaxBound
    
    !!---------------------------------------------------------------------------
    !!> @author Ricardo Birjukovs Canelas - MARETEC
    !!> @brief
    !!> Method that returns the field's maximum value (scalar)
    !!---------------------------------------------------------------------------
    !real(prec) function getFieldMaxBound(self)
    !class(field_class), intent(in) :: self
    !type(string) :: outext
    !select type(self)
    !class is (scalar_field_class)
    !class is (scalar1d_field_class)
    !    getFieldMaxBound = maxval(self%field)
    !class is (scalar2d_field_class)
    !    getFieldMaxBound = maxval(self%field)
    !class is (scalar3d_field_class)
    !    getFieldMaxBound = maxval(self%field)
    !class is (scalar4d_field_class)
    !    getFieldMaxBound = maxval(self%field)
    !class is (vectorial_field_class)
    !    getFieldMaxBound = MV
    !    class default
    !    outext = '[field_class::getFieldMaxBound]: Unexepected type of content, not a scalar or vectorial Field'
    !    call Log%put(outext)
    !    stop
    !end select
    !end function getFieldMaxBound

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Updated in Dec 2022 by Joao Sbrinho - Colab Atlantic
    !> Method that returns the field's minimum value (scalar)
    !---------------------------------------------------------------------------
    real(prec) function getFieldMinBound(self, arrayDim)
    class(field_class), intent(in) :: self
    type(string) :: outext
    integer, intent(in), optional :: arrayDim !Dimension (columns, rows, depth or time)
    select type(self)
    class is (scalar_field_class)
    class is (scalar_dimfield_class)
        if (allocated(self%field1D)) then
            getFieldMinBound = minval(self%field1D)
        elseif (allocated(self%field2D)) then
            getFieldMinBound = minval(self%field2D)
        end if
    class is (scalar1d_field_class)
        getFieldMinBound = minval(self%field)
    class is (scalar2d_field_class)
        if (present(arrayDim)) then
            if (arrayDim == 1) getFieldMinBound = minval(self%field(:,1)) !Lon
            if (arrayDim == 2) getFieldMinBound = minval(self%field(1,:)) !Lat
        else
            getFieldMinBound = minval(self%field)
        end if
    class is (scalar3d_field_class)
        getFieldMinBound = minval(self%field)
    class is (scalar4d_field_class)
        getFieldMinBound = minval(self%field)
    class is (vectorial_field_class)
        getFieldMinBound = MV
    class is (generic_field_class)
        if (allocated(self%scalar1d%field)) getFieldMinBound = minval(self%scalar1d%field)
        if (allocated(self%scalar2d%field)) getFieldMinBound = minval(self%scalar2d%field)
        if (allocated(self%scalar3d%field)) getFieldMinBound = minval(self%scalar3d%field)
        if (allocated(self%scalar4d%field)) getFieldMinBound = minval(self%scalar4d%field)
        if (allocated(self%vectorial2d%field)) getFieldMinBound = MV
        if (allocated(self%vectorial3d%field)) getFieldMinBound = MV
        if (allocated(self%vectorial4d%field)) getFieldMinBound = MV
        class default
        outext = '[field_class::getFieldMinBound]: Unexepected type of content, not a scalar or vectorial Field'
        call Log%put(outext)
        stop
    end select
    end function getFieldMinBound

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that initializes a generic field object given a field object
    !> @param[in] self, aField
    !---------------------------------------------------------------------------
    subroutine getGField(self, aField)
    class(generic_field_class), intent(inout) :: self
    class(field_class), intent(in) :: aField
    logical :: done
    type(string) :: outext
    done = .false.
    select type(aField)
    class is (scalar1d_field_class)
        call self%initialize(aField%name, aField%units, aField%field)
        done = .true.
        return
    class is (scalar2d_field_class)
        call self%initialize(aField%name, aField%units, aField%field)
        done = .true.
        return
    class is (scalar3d_field_class)
        call self%initialize(aField%name, aField%units, aField%field)
        done = .true.
        return
    class is (scalar4d_field_class)
        call self%initialize(aField%name, aField%units, aField%field)
        done = .true.
        return
    end select
    if (.not.done) then
        outext = '[generic_field_class::getGField] Unexepected type of content, not a scalar Field'
        call Log%put(outext)
        stop
    end if
    end subroutine getGField


    end module fieldTypes_mod