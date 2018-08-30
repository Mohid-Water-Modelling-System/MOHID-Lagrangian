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
    !>
    !------------------------------------------------------------------------------

    module field_types_mod

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
    end type scalar1d_field_class

    type, extends(scalar_field_class) :: scalar2d_field_class      !< a 2D scalar field class
        real(prec), allocatable, dimension(:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar2dField
    end type scalar2d_field_class

    type, extends(scalar_field_class) :: scalar3d_field_class      !< a 3D scalar field class
        real(prec), allocatable, dimension(:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar3dField
    end type scalar3d_field_class

    type, extends(scalar_field_class) :: scalar4d_field_class      !< a 4D scalar field class
        real(prec), allocatable, dimension(:,:,:,:) :: field !< the data on the scalar data field
    contains
    procedure :: initialize => initScalar4dField
    end type scalar4d_field_class

    !Vectorial fields
    !2D (t,vx)
    !3D (t,vx,vy)
    !4D (t,vx,vy,vz)

    type, extends(field_class) :: vectorial_field_class      !< a vectorial field class
    contains
    !procedure :: initialize => initVectorialField
    end type vectorial_field_class

    type, extends(vectorial_field_class) :: vectorial2d_field_class      !< a 2D vectorial field class
        type(vector), allocatable, dimension(:,:) :: field !< the data on the 2D vectorial data field
    contains
    !procedure :: initialize => initScalarField
    end type vectorial2d_field_class

    type, extends(vectorial_field_class) :: vectorial3d_field_class      !< a 3D vectorial field class
        type(vector), allocatable, dimension(:,:,:) :: field !< the data on the 3D vectorial data field
    contains
    !procedure :: initialize => initScalarField
    end type vectorial3d_field_class

    type, extends(vectorial_field_class) :: vectorial4d_field_class      !< a 4D vectorial field class
        type(vector), allocatable, dimension(:,:,:,:) :: field !< the data on the 4D vectorial data field
    contains
    !procedure :: initialize => initScalarField
    end type vectorial4d_field_class
    
    type, extends(field_class) :: generic_field_class !< generic field class. This works as a wrapper for a generic initialization routine.
        type(scalar1d_field_class) :: scalar1d       !< 1D scalar field
        type(scalar2d_field_class) :: scalar2d       !< 2D scalar field
        type(scalar3d_field_class) :: scalar3d       !< 3D scalar field
        type(scalar4d_field_class) :: scalar4d       !< 4D scalar field
        type(vectorial2d_field_class) :: vectorial2d !< 2D vectorial field
        type(vectorial3d_field_class) :: vectorial3d !< 3D vectorial field
        type(vectorial4d_field_class) :: vectorial4d !< 4D vectorial field
    contains
    procedure :: test
    procedure :: initS1D, initS2D, initS3D, initS4D
    generic   :: initialize => initS1D, initS2D, initS3D, initS4D
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
    !> Method that allocates and initializes a scalar 1D field in a generic field
    !> @parm[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initS1D(self, name, units, dim, field)
    class(generic_field_class), intent(inout) :: self
    real(prec), intent(in), dimension(:) :: field
    type(string), intent(in) :: name
    type(string), intent(in) :: units
    integer, intent(in) :: dim
    if (allocated(self%scalar1d%field)) then
        stop '[generic_field_class::initialize]: scalar 1D field already allocated'
    else
        call self%scalar1d%initialize(name, units, dim, field)
    end if
    end subroutine initS1D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 2D field in a generic field
    !> @parm[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
    subroutine initS2D(self, name, units, dim, field)
        class(generic_field_class), intent(inout) :: self
        real(prec), intent(in), dimension(:,:) :: field
        type(string), intent(in) :: name
        type(string), intent(in) :: units
        integer, intent(in) :: dim
        if (allocated(self%scalar2d%field)) then
            stop '[generic_field_class::initialize]: scalar 2D field already allocated'
        else
            call self%scalar2d%initialize(name, units, dim, field)
        end if
        end subroutine initS2D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 3D field in a generic field
    !> @parm[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
        subroutine initS3D(self, name, units, dim, field)
            class(generic_field_class), intent(inout) :: self
            real(prec), intent(in), dimension(:,:,:) :: field
            type(string), intent(in) :: name
            type(string), intent(in) :: units
            integer, intent(in) :: dim
            if (allocated(self%scalar3d%field)) then
                stop '[generic_field_class::initialize]: scalar 3D field already allocated'
            else
                call self%scalar3d%initialize(name, units, dim, field)
            end if
            end subroutine initS3D

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates and initializes a scalar 4D field in a generic field
    !> @parm[in] self, name, units, dim, field
    !---------------------------------------------------------------------------
            subroutine initS4D(self, name, units, dim, field)
                class(generic_field_class), intent(inout) :: self
                real(prec), intent(in), dimension(:,:,:,:) :: field
                type(string), intent(in) :: name
                type(string), intent(in) :: units
                integer, intent(in) :: dim
                if (allocated(self%scalar4d%field)) then
                    stop '[generic_field_class::initialize]: scalar 4D field already allocated'
                else
                    call self%scalar4d%initialize(name, units, dim, field)
                end if
                end subroutine initS4D
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 1D field
    !> @parm[in] self, name, units, dim, field
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
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that initializes a scalar 2D field
    !> @parm[in] self, name, units, dim, field
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
    !> Method that initializes a scalar 3D field
    !> @parm[in] self, name, units, dim, field
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
    !> Method that initializes a scalar 4D field
    !> @parm[in] self, name, units, dim, field
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
    !> Method that initializes a base field object by filling metadata
    !> @parm[in] self, name, units, dim
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
    !> A class 'unit' test for the generic_field_class
    !---------------------------------------------------------------------------
    subroutine test(self)
        class(generic_field_class), intent(inout) :: self
        type(generic_field_class) :: gfield1, gfield2
        real(prec), allocatable, dimension(:) :: field1
        real(prec), allocatable, dimension(:,:) :: field2
        type(string) :: name1, name2
        type(string) :: units1, units2
        integer :: dim1, dim2
        allocate(field1(50))
        allocate(field2(20,60))
        name1 = 'testfield1d'
        name2 = 'testfield2d'
        units1 = 'm/s'
        units2 = 'km'
        dim1 = 1
        dim2 = 2
        call gfield1%initialize(name1, units1, dim1, field1)
        call gfield2%initialize(name2, units2, dim2, field2)
        call gfield1%print()
        call gfield2%print()
        end subroutine test
        
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the field information
    !---------------------------------------------------------------------------
    subroutine printField(self)
    class(field_class), intent(in) :: self
    type(string) :: outext, t(2)
    t(1) = self%dim
    t(2) = self%getFieldType()
    outext = t(2)//' field['//self%name//'] has dimensionality '//t(1)//' and is in '//self%units
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


    end module field_types_mod
