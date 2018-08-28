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
    type(string) :: name        !< name of the field
    type(string) :: units       !< units of the field, preferably SI please
    integer :: dim              !< dimensions of the field (2,3 or 4)
    contains
    end type field_class

    !scalar fields
    !2D (t,x)
    !3D (t,x,y)
    !4D (t,x,y,z)
    
    type, extends(field_class) :: scalar_field_class      !< a scalar field class
    contains
    !procedure :: initialize => initScalarField
    end type scalar_field_class

    type, extends(scalar_field_class) :: scalar2d_field_class      !< a 2D scalar field class
    real(prec), allocatable, dimension(:,:) :: field !< the data on the scalar data field
    contains
    !procedure :: initialize => initScalarField
    end type scalar2d_field_class

    type, extends(scalar_field_class) :: scalar3d_field_class      !< a 3D scalar field class
    real(prec), allocatable, dimension(:,:,:) :: field !< the data on the scalar data field
    contains
    !procedure :: initialize => initScalarField
    end type scalar3d_field_class
    
    type, extends(scalar_field_class) :: scalar4d_field_class      !< a 4D scalar field class
    real(prec), allocatable, dimension(:,:,:,:) :: field !< the data on the scalar data field
    contains
    !procedure :: initialize => initScalarField
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

    !Public access vars
    public :: field_class
    public :: scalar_field_class, scalar2d_field_class, scalar3d_field_class, scalar4d_field_class
    public :: vectorial_field_class, vectorial2d_field_class, vectorial3d_field_class, vectorial4d_field_class

    contains
    
    
end module field_types_mod
