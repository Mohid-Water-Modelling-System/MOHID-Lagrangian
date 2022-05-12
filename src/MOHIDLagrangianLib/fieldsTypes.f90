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
    end type scalar_field_class

    type, extends(scalar_field_class) :: scalar1d_field_class      !< a 1D scalar field class
        real(prec), allocatable, dimension(:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar1dField
    procedure :: getFieldNearestIndex
    procedure :: finalize => cleanScalar1dField
    end type scalar1d_field_class

    type, extends(scalar_field_class) :: scalar2d_field_class      !< a 2D scalar field class
        real(prec), allocatable, dimension(:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar2dField
    procedure :: finalize => cleanScalar2dField
    end type scalar2d_field_class

    type, extends(scalar_field_class) :: scalar3d_field_class      !< a 3D scalar field class
        real(prec), allocatable, dimension(:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar3dField
    procedure :: finalize => cleanScalar3dField
    end type scalar3d_field_class

    type, extends(scalar_field_class) :: scalar4d_field_class      !< a 4D scalar field class
        real(prec), allocatable, dimension(:,:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar4dField
    procedure :: finalize => cleanScalar4dField
    end type scalar4d_field_class

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
    contains
    procedure, private :: test
    procedure :: initS1D, initS2D, initS3D, initS4D
    procedure :: initV2D, initV3D, initV4D
    generic   :: initialize => initS1D, initS2D, initS3D, initS4D, initV2D, initV3D, initV4D
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
    public :: scalar_field_class, scalar1d_field_class, scalar2d_field_class, scalar3d_field_class, scalar4d_field_class
    public :: vectorial_field_class, vectorial2d_field_class, vectorial3d_field_class, vectorial4d_field_class

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

    fDim = self%dim
    fType = self%getGFieldType()
    if (gfield%dim /= fDim) return
    if (gfield%getGFieldType() /= fType) return
    
    if (present(usedPosi)) then
        allocate(usedPos, source = usedPosi)
        nUsedPos = count(usedPos)
    end if
    
    if (fType == 'Scalar') then
        if (fDim == 1) then
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
            if (.not.allocated(usedPos)) then
                nUsedPos = size(gfield%scalar4d%field,4)
                allocate(usedPos(nUsedPos))
                usedPos = .true.
            end if
            allocate(field4d(size(self%scalar4d%field,1),size(self%scalar4d%field,2),size(self%scalar4d%field,3),size(self%scalar4d%field,4) + nUsedPos))
            field4d(:,:,:,1:size(self%scalar4d%field,4)) = self%scalar4d%field
             i= size(self%scalar4d%field, 4) +1
            do j= 1, size(usedPos)
                if (usedPos(j)) then
                    field4d(:,:,:, i) = gfield%scalar4d%field(:,:,:, j)
                    i= i+1
                end if
            end do
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
    !> Method that compares the metadata of a scalar 1D field to the object field
    !> @param[in] self, gfield
    !---------------------------------------------------------------------------
    logical function compare(self, gfield) result(comp)
    class(generic_field_class), intent(inout) :: self
    class(generic_field_class), intent(in) :: gfield
    type(string) :: outext
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
        write(*,*) "Entrei no initS1D: ", name
        call self%scalar1d%initialize(name, units, 1, field)
        call self%setFieldMetadata(name, units, 1)
        write(*,*) "Sai do initS1D: ", name
    end if
    end subroutine initS1D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the field's nearest value index (scalar)
    !> @param[in] self, value
    !---------------------------------------------------------------------------
    integer function getFieldNearestIndex(self, value)
    class(scalar1d_field_class), intent(in) :: self
    real(prec), intent(in) :: value
    real(prec), allocatable, dimension(:) :: comp
    allocate(comp(size(self%field)))
    comp = value
    getFieldNearestIndex = minloc(abs(comp - self%field), DIM=1)
    end function getFieldNearestIndex

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
        write(*,*) "Entrei no initS2D: ", name
        call self%scalar2d%initialize(name, units, 2, field)
        call self%setFieldMetadata(name, units, 2)
        write(*,*) "Sai do initS2D: ", name
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
    write(*,*) "A iniciar variavel: ", name
    call self%setFieldMetadata(name, units, dim)
    write(*,*) "finalizei inicializacao variavel: ", name
    allocate(self%field, source = field)
    end subroutine initScalar1dField

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
    !> Method that returns the field's maximum value (scalar)
    !---------------------------------------------------------------------------
    real(prec) function getFieldMaxBound(self)
    class(field_class), intent(in) :: self
    type(string) :: outext
    select type(self)
    class is (scalar_field_class)
    class is (scalar1d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (scalar2d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (scalar3d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (scalar4d_field_class)
        getFieldMaxBound = maxval(self%field)
    class is (vectorial_field_class)
        getFieldMaxBound = MV
        class default
        outext = '[field_class::getFieldMaxBound]: Unexepected type of content, not a scalar or vectorial Field'
        call Log%put(outext)
        stop
    end select
    end function getFieldMaxBound

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the field's minimum value (scalar)
    !---------------------------------------------------------------------------
    real(prec) function getFieldMinBound(self)
    class(field_class), intent(in) :: self
    type(string) :: outext
    select type(self)
    class is (scalar_field_class)
    class is (scalar1d_field_class)
        getFieldMinBound = minval(self%field)
    class is (scalar2d_field_class)
        getFieldMinBound = minval(self%field)
    class is (scalar3d_field_class)
        getFieldMinBound = minval(self%field)
    class is (scalar4d_field_class)
        getFieldMinBound = minval(self%field)
    class is (vectorial_field_class)
        getFieldMinBound = MV
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
        write(*,*) "Entrei no getGField: ", aField%name
        call self%initialize(aField%name, aField%units, aField%field)
        write(*,*) "Sai do getGField: ", aField%name
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