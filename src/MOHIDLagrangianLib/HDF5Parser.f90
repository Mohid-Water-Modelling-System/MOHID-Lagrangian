    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : hdf5Parser_mod
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : September 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a hdf5 file model class, responsible for abstracting
    !> the parsing of hdf5 MOHID files. Each object of this class is responsible for
    !> one file, effectivelly representing it.
    !> The object should return the necessary data fields with corresponding meta
    !> data, such as names and units.
    !------------------------------------------------------------------------------

    module hdf5Parser_mod

    use common_modules
    use background_mod
    use fieldTypes_mod

    use ModuleHDF5
    use ModuleTime
    use ModuleGlobalData

    implicit none
    private
    
    type :: dim_t
        type(string) :: name
        type(string) :: simName
        integer :: dimid
        integer :: varid
        type (string) :: units
        integer :: length
        logical :: reverse_axis = .false., negate = .false.
        logical :: reverse_data = .false.
    contains
    !procedure :: print => printDimsNC
    end type dim_t

    type :: var_t
        type(string) :: name
        type(string) :: simName
        integer :: varid
        type (string) :: units
        integer :: ndims
        integer, allocatable, dimension(:) :: dimids
        integer :: natts
        real(prec) :: offset, scale, fillvalue
    contains
    !procedure :: print => printVarsNC
    end type var_t

    type :: hdf5file_class !< A class that models a hdf5 mohid file
        type(string) :: filename   !< name of the file to read
        integer :: ncID             !< ID of the file
        integer :: nDims, nVars, nAtt, uDimID   !< number of dimensions, variables, attributes and dim IDs on the file
        type(dim_t), allocatable, dimension(:) :: dimData !< metadata from the dimensions on the file
        type(var_t), allocatable, dimension(:) :: varData   !< metadata from the variables on the file
        integer :: status
    contains
    !procedure :: initialize => getFile
    !procedure :: getVarDimensions
    !procedure :: getVar
    !procedure :: finalize => closeFile
    !procedure, private :: check
    !procedure, private :: getNCid
    !procedure, private :: getNCglobalMetadata
    !procedure, private :: getNCDimMetadata
    !procedure, private :: getNCVarMetadata
    !procedure, private :: getDimByDimID
    !procedure :: print => printNcInfo
    end type hdf5file_class

    !Public access vars
    public :: hdf5file_class

    contains
    
    

    end module hdf5Parser_mod