    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : hdf5parser
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Joao Sobrinho
    !
    ! DESCRIPTION:
    !> Module that defines a hdf5 file model class, responsible for abstracting
    !> the parsing of hdf5 files. Each object of this class is responsible for
    !> one file, effectivelly representing it.
    !> The object should return the necessary data fields with corresponding meta
    !> data, such as names and units.
    !------------------------------------------------------------------------------

    module hdf5Parser_mod

    use hdf5
    
#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif

    use ModuleGlobalData
    use common_modules
    use background_mod
    use fieldTypes_mod
    use xmlParser_mod

    use FoX_dom

    implicit none
    private

    type :: dim_t
        type(string) :: name
        type(string) :: simName
        integer :: dimid
        integer :: varid
        type (string) :: units
        integer :: length
        integer :: ndims = 0
        logical :: reverse_axis = .false., negate = .false.
        logical :: reverse_data = .false.
    contains
    procedure :: print => printDimsHDF5
    end type dim_t

    type :: var_t
        type(string) :: name
        type(string) :: simName
        type(string) :: hdf5GroupName
        integer :: varid
        type (string) :: units
        integer :: ndims = 0
        logical :: isDimensionVar
        integer, allocatable, dimension(:) :: dimids
        integer :: natts
        real(prec) :: offset, scale, fillvalue
    contains
    procedure :: print => printVarsHDF5
    end type var_t

    type :: hdf5file_class !< A class that models a hdf5 file
        type(string) :: filename   !< name of the file to read
        integer :: hdf5ID             !< ID of the file
        integer :: nDims = 0, nVars = 0, nAtt, uDimID   !< number of dimensions, variables, attributes and dim IDs on the file
        type(dim_t), allocatable, dimension(:) :: dimData !< metadata from the dimensions on the file
        type(var_t), allocatable, dimension(:) :: varData   !< metadata from the variables on the file
        type(string) :: mapVarName !< mapping var name
        integer :: mDims
        logical :: grid_2D
        integer :: status
    contains
    procedure :: initialize => getFile
    procedure :: getVarDimensions_variable
    procedure :: getVar
    procedure :: finalize => closeFile
    procedure, private :: check
    procedure, private :: gethdf5ID
    procedure, private :: gethdfglobalMetadata
    procedure, private :: gethdfNumberOfVarsAndMaxDims
    !procedure, private :: getNCDimMetadata
    procedure, private :: getHdfVarMetadata
    procedure, private :: hdfReadAllVariables
    procedure, private :: getDimByDimID
    procedure, private :: getMappingVar
    procedure, private :: check2dGrid
    procedure, private :: readHDFVariable
    procedure, private :: readHDFIntVariable
    procedure, private :: readHDFTime
    !procedure :: print => printNcInfo
    end type hdf5file_class
    
    type :: hdf5Reader_class
    contains
    procedure :: getFullFile
    end type hdf5Reader_class

    public :: hdf5file_class, hdf5Reader_class

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a NC file
    !> @param[in] self, fileName, varList, syntecticVar
    !---------------------------------------------------------------------------
    type(background_class) function getFullFile(self, fileName, varList, syntecticVar)
    class(hdf5Reader_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(string), dimension(:), intent(in) :: varList
    logical, dimension(:), intent(in) :: syntecticVar    
    type(hdf5file_class) :: hdf5File
    type(generic_field_class), allocatable, dimension(:) :: backgrounDims
    type(generic_field_class), allocatable, dimension(:) :: gfield
    type(string) :: name, units
    type(box) :: extents
    type(vector) :: pt
    real(prec), dimension(3,2) :: dimExtents
    integer :: i, realVarIdx, mapVarID
    type(string) :: outext
    
    allocate(gfield(size(syntecticVar)))
    realVarIdx = 0
    units = '-'
    outext = '->Reading '//fileName
    call Log%put(outext,.false.)
    
    call hdf5File%initialize(fileName)
    do i=1, size(syntecticVar)
        if(.not.syntecticVar(i)) then !finding the first real variable to extract dimension arrays
            call hdf5File%getVarDimensions_variable(varList(i), backgrounDims)
            
            if (allocated(backgrounDims)) then
                realVarIdx = i
                exit
            end if 
        end if
    end do

    if (realVarIdx /= 0) then
        write(*,*)"Getting MapVar"
        !Allocates and saves the mapping var using velU as source for dimensions
        mapVarID = size(syntecticVar)
        call hdf5File%getMappingVar(varList(realVarIdx), gfield(mapVarID), varList(mapVarID), units)
        do i=1, size(syntecticVar)-1
            write(*,*)"Getting Var name = ", varList(i)
            call hdf5File%getVar(varList(mapVarID), gfield(i), gfield(mapVarID), .true., varList(i), units)
        end do
    end if
    
    call hdf5File%finalize()
    
    dimExtents = 0.0
    
    do i = 1, size(backgrounDims)
        !write(*,*) "backgroundDims name = ", trim(backgrounDims(i)%name)
        if (backgrounDims(i)%name == Globals%Var%lon) then
            dimExtents(1,1) = backgrounDims(i)%getFieldMinBound()
            !write(*,*) "dimExtents(1,1) = ", dimExtents(1,1)
            dimExtents(1,2) = backgrounDims(i)%getFieldMaxBound()
            !write(*,*) "dimExtents(1,2) = ", dimExtents(1,2)
        else if (backgrounDims(i)%name == Globals%Var%lat) then
            dimExtents(2,1) = backgrounDims(i)%getFieldMinBound()
            !write(*,*) "dimExtents(2,1) = ", dimExtents(2,1)
            dimExtents(2,2) = backgrounDims(i)%getFieldMaxBound()
            !write(*,*) "dimExtents(2,2) = ", dimExtents(2,2)
        else if (backgrounDims(i)%name == Globals%Var%level) then
            dimExtents(3,1) = backgrounDims(i)%getFieldMinBound()
            !write(*,*) "dimExtents(3,1) = ", dimExtents(3,1)
            dimExtents(3,2) = backgrounDims(i)%getFieldMaxBound()
            !write(*,*) "dimExtents(3,2) = ", dimExtents(3,2)
        end if
    end do
    extents%pt = dimExtents(1,1)*ex + dimExtents(2,1)*ey + dimExtents(3,1)*ez
    pt = dimExtents(1,2)*ex + dimExtents(2,2)*ey + dimExtents(3,2)*ez
    extents%size = pt - extents%pt
    name = fileName%basename(strip_last_extension=.true.)
    !For some reason this is not working
    getFullFile = Background(1, name, extents, backgrounDims)
    !write(*,*)"size gfield = ", size(gfield)
    do i = 1, size(gfield)
        !write(*,*) "Entrei gfield add = ", i
        call getFullFile%add(gfield(i))
        !write(*,*) "Sai gfield add = ", i
    end do
    !write(*,*) "Sai do getFullFile"
    end function getFullFile

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Parse the hdf5 file headers and assembles the file model
    !> @param[in] self, filename
    !---------------------------------------------------------------------------
    subroutine getFile(self, filename)
    class(hdf5file_class), intent(inout) :: self
    type(string), intent(in) :: filename
    self%filename = filename
    call self%gethdf5ID()
    call self%gethdfglobalMetadata()
    call self%getHdfVarMetadata()
    call self%check2dGrid()
    end subroutine getFile

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Opens the hdf5 file, assigns the hdf5ID and checks for errors
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine gethdf5ID(self)
    class(hdf5file_class), intent(inout) :: self
    call h5fopen_f (trim(self%filename%chars()), access_flags = H5F_ACC_RDONLY_F,               &
                    file_id = self%hdf5ID, HDFERR = self%status)
    
    call self%check()
    end subroutine gethdf5ID

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Inquires the hdf file for global metadata - number of dimensions,
    !> and variables
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine gethdfglobalMetadata(self)
    class(hdf5file_class), intent(inout) :: self
    integer(HID_T)                     :: gr_idIn
    character(len=1)                   :: GroupNameIn
    integer                            :: STAT
    
    GroupNameIn = "/"
            
    call h5gopen_f        (self%hdf5ID, trim(GroupNameIn), gr_idIn, STAT)
    if (STAT /= SUCCESS_) then
        stop 'GetHDF5AllDataSetsOK - ModuleHDF5 - ERR10'
    endif
    
    call self%gethdfNumberOfVarsAndMaxDims (IDIn = gr_idIn, GroupNameIn = GroupNameIn)
            
    call h5gclose_f       (gr_idIn, STAT)
    if (STAT /= SUCCESS_) then
        stop 'GetHDF5AllDataSetsOK - ModuleHDF5 - ERR20'
    endif 
    
    end subroutine gethdfglobalMetadata

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Inquires the hdf file for variable metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getHdfVarMetadata(self)
    class(hdf5file_class), intent(inout) :: self
    integer(HID_T)                       :: gr_idIn
    character(len=1)                     :: GroupNameIn
    integer                              :: STAT, VarCounter, DimCounter
    !Begin--------------------------------------------------------------------------------
    
    allocate(self%varData(self%nVars))
    allocate(self%dimData(self%nDims+1)) !SOBRINHO - VER SE self%nDims ou self%nDims+1 para incluir tempo
    GroupNameIn = "/"
            
    call h5gopen_f        (self%hdf5ID, trim(GroupNameIn), gr_idIn, STAT)
    if (STAT /= SUCCESS_) then
        stop 'GetHDF5AllDataSetsOK - ModuleHDF5 - ERR10'
    endif
    
    VarCounter = 1
    DimCounter = 1
    call self%hdfReadAllVariables (gr_idIn, GroupNameIn, VarCounter, DimCounter)

    call h5gclose_f       (gr_idIn, STAT)
    if (STAT /= SUCCESS_) then
        stop 'GetHDF5AllDataSetsOK - ModuleHDF5 - ERR20'
    endif 
    
    end subroutine getHdfVarMetadata

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Inquires the nc file for 2D grid coordinates
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine check2dGrid(self)
    class(hdf5file_class), intent(inout) :: self
    integer :: i
    self%grid_2D = .false.
    do i=1, self%nVars
        if ((self%varData(i)%ndims == 2) .and. (Globals%Var%checkDimensionName(self%varData(i)%name))) then
            !Current variable is 2D and is a dimension variable
            self%grid_2D = .true.
            return
        end if
    end do
    end subroutine check2dGrid

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Reads the dimension fields from the nc file for a given variable.
    !> returns an array of scalar 1D fields, each with a name, units and data (includes Lat, Lon and Depth)
    !> @param[in] self, varName, dimsArrays
    !---------------------------------------------------------------------------
    subroutine getVarDimensions_variable(self, varName, dimsArrays)
    class(hdf5file_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(generic_field_class), allocatable, dimension(:), intent(out) :: dimsArrays
    real(prec), allocatable, dimension(:) :: tempRealArray1D
    real(prec), allocatable, dimension(:,:) :: tempRealArray2D
    real(prec), allocatable, dimension(:,:,:) :: tempRealArray3D
    real(prec), allocatable, dimension(:,:,:,:) :: tempRealArray4D
    type(string) :: dimName, dimUnits, testeString
    integer, allocatable, dimension(:) :: varShape
    integer :: i, j, k, ndim, var, t
    type(dim_t) :: tempDim
    logical :: foundDimVar, found2dDimVar, isSpaceDimension
    !Begin-----------------------------------------------------------------------
    write(*,*)"entrei getVarDimensions_variable"
    do i=1, self%nVars !going trough all variables
        write(*,*)"i = ", i
        write(*,*)"var name = ", varName
        write(*,*)"sim name = ", self%varData(i)%simName
        if (self%varData(i)%simName == varName) then   !found the requested var
            write(*,*)"number of dims in varname = ", self%varData(i)%ndims
            allocate(dimsArrays(self%varData(i)%ndims)) !allocating dimension fields
            allocate(varShape(self%varData(i)%ndims))
            write(*,*)"allocate done"
            !Get vector dimension vectors
            do ndim =1, self%varData(i)%ndims   !going trough all of the variable dimensions
                write(*,*)"ndim =  ", ndim
                tempDim = self%getDimByDimID(self%varData(i)%dimids(ndim))
                varShape(ndim) = tempDim%length
            end do
            write(*,*)"varshape done"
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                write(*,*)"dims dimension = ", self%nDims
                foundDimVar = .false.
                do k=1, self%nDims  !going trough the dimensions of the file
                    write(*,*)"enteri novo k "
                    write(*,*)"dim ID da dimensao atual da variavel =  ", self%varData(i)%dimids(j)
                    write(*,*)"dim ID do dimData =  ", self%dimData(k)%dimid
                    write(*,*)"nome da dimensao (k) =  ", trim(self%dimData(k)%name)
                    if (self%varData(i)%dimids(j) == self%dimData(k)%dimid) then    !found a corresponding dimension between the variable and the file
                        found2dDimVar = .false.
                        write(*,*)"enteri if self%varData(i)%dimids(j) == self%dimData(k)%dimid "
                        isSpaceDimension = Globals%Var%checkDimensionName(self%dimData(k)%name) !Lat and lon have the same dimensions so dont need to check if it is lat or lon
                cd1:    if (self%grid_2D .and. isSpaceDimension) then
                    dvar:   do var = 1, self%nVars !going trough all file variables
                                isSpaceDimension = .false.
                                write(*,*)"Var to check = ", TRIM(self%varData(var)%name)
                                write(*,*)"dimData name = ", TRIM(self%dimData(k)%name)
                                testeString = Globals%Var%getVarSimName(self%dimData(k)%name)
                                write(*,*)"DimData SimName = ", TRIM(testeString)
                                write(*,*)"Var SimName = ", TRIM(self%varData(var)%simName)
                                
                                isSpaceDimension = Globals%Var%checkDimensionName(self%varData(var)%name)
                                write(*,*)"isSpaceDimension = ", isSpaceDimension
                                if ((isSpaceDimension) .and. (Globals%Var%getVarSimName(self%dimData(k)%name)==self%varData(var)%simName)) then
                                    !Lat or Lon found, so allocate 2D field
                                    write(*,*)"varShape(1) = ", varShape(1)
                                    write(*,*)"varShape(2) = ", varShape(2)
                                    if (self%dimData(k)%ndims == 2) then
                                        allocate(tempRealArray2D(varShape(1), varShape(2))) !allocating a place to read the field data to
                                        dimName = self%varData(var)%simName
                                        write(*,*)"simulation name of dimension = ", dimName
                                        dimUnits = self%varData(var)%units
                                        write(*,*)"dimUnits of dimension = ", dimUnits
                                        !Here is where Lat, Lon and Depth are read from the hdf5
                                        call self%readHDFVariable(self%varData(var), array2D = tempRealArray2D)
                                    elseif (self%dimData(k)%ndims == 4) then
                                        !Has to be VerticalZ which is 4D (3D + time)
                                        allocate(tempRealArray4D(varShape(1)-1, varShape(2)-1, varShape(3), varShape(4))) !allocating a place to read the field data to
                                        allocate(tempRealArray3D(varShape(1)-1, varShape(2)-1, varShape(3)))
                                        dimName = self%varData(var)%simName
                                        write(*,*)"simulation name of dimension = ", dimName
                                        dimUnits = self%varData(var)%units
                                        write(*,*)"dimUnits of dimension = ", dimUnits
                                        !Here is where Lat, Lon and Depth are read from the hdf5
                                        do t=1, varShape(4)
                                            !Depth is assumed to be VerticalZ which is a 4D var.
                                            call self%readHDFVariable(self%varData(var), array3D = tempRealArray3D, outputNumber = t)
                                            tempRealArray4D(:,:,:,t) = - tempRealArray3D
                                        enddo
                                        
                                        write(*,*) "Min Val = ", MinVal(tempRealArray4D)
                                        where (tempRealArray4D > -Globals%Parameters%FillValueReal / 2.0)
                                            tempRealArray4D = 0.0
                                        endwhere
                                        write(*,*) "Max Val = ", MaxVal(tempRealArray4D)
                                    else
                                        allocate(tempRealArray1D(self%dimData(k)%length)) !allocating a place to read the field data to
                                        dimName = self%dimData(k)%simName
                                        write(*,*)"Dim name = ", dimName
                                        dimUnits = self%dimData(k)%units
                                        write(*,*)"dimUnits = ", dimUnits
                                        write(*,*)"self%dimData(k)%varid = ", self%dimData(k)%varid
                            
                                        call self%readHDFTime(self%varData(var), tempRealArray1D, varShape(4))
                                    endif
                                    
                                    found2dDimVar = .true.
                                    exit cd1
                                end if
                            end do dvar
                        end if cd1
                        
                        write(*,*)"found2dDimVar = ", found2dDimVar
                        !for 1D time or depth
                        !if (.not. found2dDimVar) then
                        !    allocate(tempRealArray1D(self%dimData(k)%length)) !allocating a place to read the field data to
                        !    dimName = self%dimData(k)%simName
                        !    write(*,*)"Dim name = ", dimName
                        !    dimUnits = self%dimData(k)%units
                        !    write(*,*)"dimUnits = ", dimUnits
                        !    write(*,*)"self%dimData(k)%varid = ", self%dimData(k)%varid
                        !    
                        !    call self%readHDFTime(self%varData(var), tempRealArray1D, varShape(4))
                        !    
                        !    !need to check for 'time' variable specific issues
                        !    !if (self%dimData(k)%simName == Globals%Var%time) then
                        !    !    call correctHDF5Time(dimUnits, tempRealArray1D)
                        !    !end if
                        !end if
                        
                        if (allocated(tempRealArray1D)) then
                            call dimsArrays(j)%initialize(dimName, dimUnits, tempRealArray1D)
                            deallocate(tempRealArray1D)
                            foundDimVar = .true.
                        elseif (allocated(tempRealArray2D)) then
                            call dimsArrays(j)%initialize(dimName, dimUnits, tempRealArray2D)
                            deallocate(tempRealArray2D)
                            foundDimVar = .true.
                        elseif (allocated(tempRealArray4D)) then
                            call dimsArrays(j)%initialize(dimName, dimUnits, tempRealArray4D)
                            deallocate(tempRealArray4D)
                            foundDimVar = .true.
                        end if
                        dimsArrays(j)%name = Globals%Var%getVarSimName(self%dimData(k)%name)
                        write(*,*)"nome gravado = ", trim(dimsArrays(j)%name)
                    end if
                end do
            end do
        end if
    end do
    write(*,*)"sai getVarDimensions_variable"
    end subroutine getVarDimensions_variable
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Reads the fields from the nc file for a given variable.
    !> returns a generic field, with a name, units and data
    !> @param[in] self, varName, varField, mapField, binaryVar, altName, altUnits
    !---------------------------------------------------------------------------
    subroutine getVar(self, varName, varField, mapField, binaryVar, altName, altUnits)
    class(hdf5file_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(generic_field_class), intent(out) :: varField
    type(generic_field_class), intent(in) :: mapField
    logical, optional, intent(in) :: binaryVar
    type(string), optional, intent(in) :: altName, altUnits
    logical :: bVar
    real(prec), allocatable, dimension(:,:) :: tempRealField2D
    real(prec), allocatable, dimension(:,:,:) :: tempRealField3D
    real(prec), allocatable, dimension(:,:,:,:) :: tempRealField4D
    type(string) :: dimName, varUnits
    integer :: i, j, k, t, indx, j2
    type(dim_t) :: tempDim, uDim
    integer, allocatable, dimension(:) :: varShape, u_Shape
    type(string) :: outext
    logical variable_u_is4D
    
    bVar= .false.
    if(present(binaryVar)) bVar = binaryVar
    variable_u_is4D = .false.

    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName ) then   !found the requested var
            write(*,*)"Getting Var for simulation name = ", self%varData(i)%simName
            write(*,*)"Getting Var name = ", varName
            allocate(varShape(self%varData(i)%ndims))
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                tempDim = self%getDimByDimID(self%varData(i)%dimids(j))
                varShape(j) = tempDim%length
            end do
            
            if(self%varData(i)%ndims == 3) then !3D variable
                allocate(tempRealField2D(varShape(1)-1, varShape(2)-1))
                allocate(tempRealField3D(varShape(1)-1, varShape(2)-1,varShape(3)))
                do t=1, varShape(3)
                    call self%readHDFVariable(self%varData(i), array2D = tempRealField2D, outputNumber = t)
                    tempRealField3D(:,:,t) = tempRealField2D
                enddo
                if (.not.bVar) then
                    where (mapField%intScalar3D%field == 0)
                        tempRealField3D = 0.0
                    end where
                else
                    if (self%varData(i)%fillvalue == MV) then
                        outext = '[hdf5Parser::getVar]:WARNING - variables without _fillvalue, you might have some problems in a few moments. Masks will not work properly (beaching, land exclusion,...)'
                    call Log%put(outext)
                    end if
                    where (mapField%intScalar3D%field == 1) tempRealField3D = Globals%Mask%waterVal
                    where (mapField%intScalar3D%field == 0) tempRealField3D = Globals%Mask%landVal
                end if
                if (.not.bVar) then
                    call varField%initialize(varName, self%varData(i)%units, tempRealField3D)
                else
                    dimName = varName
                    if(present(altName)) dimName = altName
                    varUnits = self%varData(i)%units
                    if(present(altUnits)) varUnits = altUnits
                    call varField%initialize(dimName, varUnits, tempRealField3D)
                end if
            else if(self%varData(i)%ndims == 4) then !4D variable                
                allocate(tempRealField4D(varShape(1)-1, varShape(2)-1, varShape(3), varShape(4))) !allocating a place to read the field data to
                allocate(tempRealField3D(varShape(1)-1, varShape(2)-1, varShape(3)))
                do t=1, varShape(4)
                    !Depth is assumed to be VerticalZ which is a 4D var.
                    call self%readHDFVariable(self%varData(i), array3D = tempRealField3D, outputNumber = t)
                    tempRealField4D(:,:,:,t) = tempRealField3D
                enddo
                
                if (.not.bVar) then
                    where (mapField%intScalar4D%field == 0)
                        tempRealField4D = 0.0
                    end where
                else
                    if (self%varData(i)%fillvalue == MV) then
                        outext = '[hdf5Parser::getVar]:WARNING - variables without _fillvalue, you might have some problems in a few moments. Masks will not work properly (beaching, land exclusion,...)'
                        call Log%put(outext)
                    end if
                    !Aqui sera onde se incluirao os openpoints
                    where (mapField%intScalar4D%field == 1) tempRealField4D = Globals%Mask%waterVal
                    where (mapField%intScalar4D%field == 0) tempRealField4D = Globals%Mask%landVal
                end if
                
                if (.not.bVar) then
                    call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
                else
                    dimName = varName
                    if(present(altName)) dimName = altName
                    varUnits = self%varData(i)%units
                    if(present(altUnits)) varUnits = altUnits
                    call varField%initialize(dimName, varUnits, tempRealField4D)
                end if
            elseif(self%varData(i)%ndims == 2) then !2D variable, for now only bathymetry
                if (self%varData(i)%simName == Globals%Var%bathymetry) then
                    allocate(tempRealField2D(varShape(1)-1,varShape(2)-1))
do1:                do indx=1, self%nVars
                        !Find velocity u matrix to get its dimensions
                        if (self%varData(indx)%simName == Globals%Var%u) then
                            allocate(u_Shape(self%varData(indx)%ndims))
                            do j2=1, self%varData(indx)%ndims   !going trough all of the variable dimensions
                                uDim = self%getDimByDimID(self%varData(indx)%dimids(j2))
                                u_Shape(j2) = uDim%length
                            end do
                            if (self%varData(indx)%ndims == 4) then
                                variable_u_is4D = .true.
                                allocate(tempRealField4D(varShape(1)-1,varShape(2)-1, u_Shape(3), u_Shape(4)))
                                exit do1
                            else
                                allocate(tempRealField3D(varShape(1)-1,varShape(2)-1, u_Shape(3)))
                                exit do1
                            end if
                            
                        end if
                    end do do1
                    
                    call self%readHDFVariable(self%varData(i), array2D = tempRealField2D)
                    
                    if (.not.bVar) then
                        !Leaving this commented because I am assuming only bathymetry gets in here and it is better that it has a fillvaluereal
                        ! value instead of 0.0
                        !where (tempRealField2D /= self%varData(i)%fillvalue)
                        !    tempRealField2D = tempRealField2D*self%varData(i)%scale + self%varData(i)%offset
                        !elsewhere (tempRealField2D == self%varData(i)%fillvalue)
                        !    tempRealField2D = 0.0
                        !end where
                        if (variable_u_is4D) then
                            !For bathymetry, converts the 2D input field into a 4D field to be consistent with velocity matrixes
                            do t=1,size(tempRealField4D,4)
                                do k=1, size(tempRealField4D,3)
                                    tempRealField4D(:,:,k,t) = - tempRealField2D(:,:)
                                end do
                            end do
                            call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
                        else
                            do k=1, size(tempRealField4D,3)
                                tempRealField3D(:,:,k) = - tempRealField2D(:,:)
                            end do
                            call varField%initialize(varName, self%varData(i)%units, tempRealField3D)
                        end if
                        
                    else
                        outext = '[hdf5parser::getVar]: Variable '//varName//' is synthetic and cannot have 2D dimensionality. Stopping'
                        call Log%put(outext)
                        stop
                    end if
                else
                    outext = '[hdf5parser::getVar]: Variable '//varName//' is 2D and not bathymetry, so it is not supported. Stopping'
                    call Log%put(outext)
                    stop
                end if
            else
                outext = '[hdf5parser::getVar]: Variable '//varName//' has a non-supported dimensionality. Stopping'
                call Log%put(outext)
                stop
            end if
        end if
    end do

    end subroutine getVar

!    !---------------------------------------------------------------------------
!    !> @author Ricardo Birjukovs Canelas - MARETEC
!    !> @brief
!    !> Reads the fields from the nc file for a given variable.
!    !> returns a generic field, with a name, units and data
!    !> @param[in] self, varName, varField, binaryVar, altName, altUnits
!    !---------------------------------------------------------------------------
!    subroutine getVar(self, varName, varField, binaryVar, altName, altUnits)
!    class(hdf5file_class), intent(inout) :: self
!    type(string), intent(in) :: varName
!    type(generic_field_class), intent(out) :: varField
!    logical, optional, intent(in) :: binaryVar
!    type(string), optional, intent(in) :: altName, altUnits
!    logical :: bVar
!    real(prec), allocatable, dimension(:) :: tempRealField1D
!    real(prec), allocatable, dimension(:,:) :: tempRealField2D
!    real(prec), allocatable, dimension(:,:,:) :: tempRealField3D
!    real(prec), allocatable, dimension(:,:,:,:) :: tempRealField4D
!    type(string) :: dimName, varUnits
!    integer :: i, j, k, id_dim, first,last, t, indx, j2
!    type(dim_t) :: tempDim, uDim
!    integer, allocatable, dimension(:) :: varShape, u_Shape
!    type(string) :: outext
!    logical variable_u_is4D
!    
!    bVar= .false.
!    if(present(binaryVar)) bVar = binaryVar
!    variable_u_is4D = .false.
!
!    do i=1, self%nVars !going trough all variables
!        if (self%varData(i)%simName == varName ) then   !found the requested var
!            write(*,*)"Getting Var for simulation name = ", self%varData(i)%simName
!            write(*,*)"Getting Var name = ", varName
!            allocate(varShape(self%varData(i)%ndims))
!            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
!                tempDim = self%getDimByDimID(self%varData(i)%dimids(j))
!                varShape(j) = tempDim%length
!            end do
!            
!            if(self%varData(i)%ndims == 3) then !3D variable
!                allocate(tempRealField2D(varShape(1)-1, varShape(2)-1))
!                allocate(tempRealField3D(varShape(1)-1, varShape(2)-1,varShape(3)))
!                do t=1, varShape(3)
!                    call self%readHDFVariable(self%varData(i), array2D = tempRealField2D, outputNumber = t)
!                    tempRealField3D(:,:,t) = tempRealField2D
!                enddo
!                if (.not.bVar) then
!                    where (tempRealField3D == self%varData(i)%fillvalue / 2.0)
!                        tempRealField3D = 0.0
!                    end where
!                else
!                    if (self%varData(i)%fillvalue == MV) then
!                        outext = '[hdf5Parser::getVar]:WARNING - variables without _fillvalue, you might have some problems in a few moments. Masks will not work properly (beaching, land exclusion,...)'
!                    call Log%put(outext)
!                    end if
!                    where (tempRealField3D > self%varData(i)%fillvalue / 2.0) tempRealField3D = Globals%Mask%waterVal
!                    where (tempRealField3D < self%varData(i)%fillvalue / 2.0) tempRealField3D = Globals%Mask%landVal
!                end if
!                if (.not.bVar) then
!                    call varField%initialize(varName, self%varData(i)%units, tempRealField3D)
!                else
!                    dimName = varName
!                    if(present(altName)) dimName = altName
!                    varUnits = self%varData(i)%units
!                    if(present(altUnits)) varUnits = altUnits
!                    call varField%initialize(dimName, varUnits, tempRealField3D)
!                end if
!            else if(self%varData(i)%ndims == 4) then !4D variable                
!                allocate(tempRealField4D(varShape(1)-1, varShape(2)-1, varShape(3), varShape(4))) !allocating a place to read the field data to
!                allocate(tempRealField3D(varShape(1)-1, varShape(2)-1, varShape(3)))
!                do t=1, varShape(4)
!                    !Depth is assumed to be VerticalZ which is a 4D var.
!                    call self%readHDFVariable(self%varData(i), array3D = tempRealField3D, outputNumber = t)
!                    tempRealField4D(:,:,:,t) = tempRealField3D
!                enddo
!                
!                if (.not.bVar) then
!                    where (tempRealField4D < self%varData(i)%fillvalue / 2.0)
!                        tempRealField4D = 0.0
!                    end where
!                else
!                    if (self%varData(i)%fillvalue == MV) then
!                        outext = '[hdf5Parser::getVar]:WARNING - variables without _fillvalue, you might have some problems in a few moments. Masks will not work properly (beaching, land exclusion,...)'
!                    call Log%put(outext)
!                    end if
!                    !Aqui sera onde se incluirao os openpoints
!                    where (tempRealField4D > self%varData(i)%fillvalue / 2.0) tempRealField4D = Globals%Mask%waterVal
!                    where (tempRealField4D < self%varData(i)%fillvalue / 2.0) tempRealField4D = Globals%Mask%landVal
!                end if
!                
!                if (.not.bVar) then
!                    call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
!                else
!                    dimName = varName
!                    if(present(altName)) dimName = altName
!                    varUnits = self%varData(i)%units
!                    if(present(altUnits)) varUnits = altUnits
!                    call varField%initialize(dimName, varUnits, tempRealField4D)
!                end if
!            elseif(self%varData(i)%ndims == 2) then !2D variable, for now only bathymetry
!                if (self%varData(i)%simName == Globals%Var%bathymetry) then
!                    allocate(tempRealField2D(varShape(1)-1,varShape(2)-1))
!do1:                do indx=1, self%nVars
!                        !Find velocity u matrix to get its dimensions
!                        if (self%varData(indx)%simName == Globals%Var%u) then
!                            allocate(u_Shape(self%varData(indx)%ndims))
!                            do j2=1, self%varData(indx)%ndims   !going trough all of the variable dimensions
!                                uDim = self%getDimByDimID(self%varData(indx)%dimids(j2))
!                                u_Shape(j2) = uDim%length
!                            end do
!                            if (self%varData(indx)%ndims == 4) then
!                                variable_u_is4D = .true.
!                                allocate(tempRealField4D(varShape(1)-1,varShape(2)-1, u_Shape(3), u_Shape(4)))
!                                exit do1
!                            else
!                                allocate(tempRealField3D(varShape(1)-1,varShape(2)-1, u_Shape(3)))
!                                exit do1
!                            end if
!                            
!                        end if
!                    end do do1
!                    
!                    call self%readHDFVariable(self%varData(i), array2D = tempRealField2D)
!                    
!                    if (.not.bVar) then
!                        !Leaving this commented because I am assuming only bathymetry gets in here and it is better that it has a fillvaluereal
!                        ! value instead of 0.0
!                        !where (tempRealField2D /= self%varData(i)%fillvalue)
!                        !    tempRealField2D = tempRealField2D*self%varData(i)%scale + self%varData(i)%offset
!                        !elsewhere (tempRealField2D == self%varData(i)%fillvalue)
!                        !    tempRealField2D = 0.0
!                        !end where
!                        if (variable_u_is4D) then
!                            !For bathymetry, converts the 2D input field into a 4D field to be consistent with velocity matrixes
!                            do t=1,size(tempRealField4D,4)
!                                do k=1, size(tempRealField4D,3)
!                                    tempRealField4D(:,:,k,t) = - tempRealField2D(:,:)
!                                end do
!                            end do
!                            call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
!                        else
!                            do k=1, size(tempRealField4D,3)
!                                tempRealField3D(:,:,k) = - tempRealField2D(:,:)
!                            end do
!                            call varField%initialize(varName, self%varData(i)%units, tempRealField3D)
!                        end if
!                        
!                    else
!                        outext = '[hdf5parser::getVar]: Variable '//varName//' is synthetic and cannot have 2D dimensionality. Stopping'
!                        call Log%put(outext)
!                        stop
!                    end if
!                else
!                    outext = '[hdf5parser::getVar]: Variable '//varName//' is 2D and not bathymetry, so it is not supported. Stopping'
!                    call Log%put(outext)
!                    stop
!                end if
!            else
!                outext = '[hdf5parser::getVar]: Variable '//varName//' has a non-supported dimensionality. Stopping'
!                call Log%put(outext)
!                stop
!            end if
!        end if
!    end do
!
!    end subroutine getVar
!    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Reads the fields from the hdf file for a map variable.
    !> returns a generic field, with a name, units and data
    !> @param[in] self, varName, varField, binaryVar, altName, mapVarUnits
    !---------------------------------------------------------------------------
    subroutine getMappingVar(self, varName, varField, mapVarName, mapVarUnits)
    class(hdf5file_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(generic_field_class), intent(out) :: varField
    type(string), optional, intent(in) :: mapVarName, mapVarUnits
    integer, allocatable, dimension(:,:) :: tempIntegerField2D
    integer, allocatable, dimension(:,:,:) :: tempIntegerField3D
    integer, allocatable, dimension(:,:,:,:) :: tempIntegerField4D
    type(string) :: dimName, varUnits
    integer :: i, j, t, var
    type(dim_t) :: tempDim
    integer, allocatable, dimension(:) :: varShape
    type(string) :: outext
    !Begin----------------------------------------------------------------------
    
    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName ) then   !found the requested var
            write(*,*)"Getting Var for simulation name = ", self%varData(i)%simName
            write(*,*)"Getting Var name = ", varName
            allocate(varShape(self%varData(i)%ndims))
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                tempDim = self%getDimByDimID(self%varData(i)%dimids(j))
                varShape(j) = tempDim%length
            end do
            
            if(self%varData(i)%ndims == 3) then !3D variable
                allocate(tempIntegerField2D(varShape(1)-1, varShape(2)-1))
                allocate(tempIntegerField3D(varShape(1)-1, varShape(2)-1,varShape(3)))
                
                do var = 1, self%nVars
                    if (self%varData(var)%simName == mapVarName) then
                        if (self%varData(var)%ndims == 3) then
                            !map var has the same dimensions as vel U (probably openpoints) and one of the dims is time
                            do t=1, varShape(3)
                                call self%readHDFIntVariable(self%varData(var), array2D = tempIntegerField2D, outputNumber = t)
                                tempIntegerField3D(:,:,t) = tempIntegerField2D
                            enddo
                        elseif (self%varData(var)%ndims == 2) then
                            !map var has one less dimensions as vel U (probably openpoints). must be a 2D waterpoints matrix
                            call self%readHDFIntVariable(self%varData(i), array2D = tempIntegerField2D)
                            tempIntegerField3D(:,:,1:varShape(3)) = tempIntegerField3D
                        else
                            outext = '[hdf5parser::getMappingVar]: Variable '//mapVarName//' has a non-supported dimensionality. Stopping'
                            call Log%put(outext)
                            stop
                        end if
                    endif
                enddo
                
                varUnits = mapVarUnits
                call varField%initialize(mapVarName, varUnits, tempIntegerField3D)
                
            else if(self%varData(i)%ndims == 4) then !4D variable                
                allocate(tempIntegerField4D(varShape(1)-1, varShape(2)-1, varShape(3), varShape(4))) !allocating a place to read the field data to
                allocate(tempIntegerField3D(varShape(1)-1, varShape(2)-1, varShape(3)))
                
                do var = 1, self%nVars
                    if (self%varData(var)%simName == mapVarName) then
                        if (self%varData(var)%ndims == 3) then
                            !map var has one less dimensions as vel U (probably openpoints) meaning it is the static(in time) waterpoints
                            call self%readHDFIntVariable(self%varData(var), array3D = tempIntegerField3D)
                            do t=1, varShape(4)
                                tempIntegerField4D(:,:,:,t) = tempIntegerField3D
                            enddo
                        elseif (self%varData(var)%ndims == 4) then
                            !map var has the same dimensions as vel U (probably openpoints)
                            do t=1, varShape(4)
                                call self%readHDFIntVariable(self%varData(i), array3D = tempIntegerField3D, outputNumber = t)
                                tempIntegerField4D(:,:,:,t) = tempIntegerField3D
                            enddo
                        else
                            
                            outext = '[hdf5parser::getMappingVar]: Variable '//mapVarName//' has a non-supported dimensionality. Stopping'
                            call Log%put(outext)
                            stop
                        end if
                    endif
                enddo

                dimName = varName
                varUnits = mapVarUnits
                call varField%initialize(mapVarName, varUnits, tempIntegerField4D)
            else
                outext = '[hdf5parser::getMappingVar]: Variable '//varName//' has a non-supported dimensionality. Stopping'
                call Log%put(outext)
                stop
            end if
        end if
    end do

    end subroutine getMappingVar
    
!
!    !---------------------------------------------------------------------------
!    !> @author Ricardo Birjukovs Canelas - MARETEC
!    !> @brief
!    !> returns a dimension field metadata structure, given a dimID
!    !> @param[in] self, dimID
!    !---------------------------------------------------------------------------
    type(dim_t) function getDimByDimID(self, dimID)
    class(hdf5file_class), intent(inout) :: self
    integer, intent(in) :: dimID
    integer :: i
    logical :: found
    type(string) :: outext

    found = .false.
    do i=1, self%nDims
        if (dimID == self%dimData(i)%dimid) then
            found = .true.
            getDimByDimID = self%dimData(i)
            return
        end if
    end do
    if (.not.found) then
        outext = dimID
        outext = '[hdf5parser::getDimByDimID]: dimension with ID='//outext//' not found. Stopping'
        call Log%put(outext)
        stop
    end if

    end function getDimByDimID
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Close the hdf5 file
    !> @param[in] self
    !---------------------------------------------------------------------------
    recursive subroutine gethdfNumberOfVarsAndMaxDims (self, IDIn, GroupNameIn)

    !Arguments-------------------------------------------------------------
    class(hdf5file_class), intent(inout) :: self
    integer(HID_T), intent(in)         :: IDIn
    character(len=*), intent(in)       :: GroupNameIn

    !Local-----------------------------------------------------------------
    integer                                     :: nmembersIn
    character(StringLength)                     :: obj_nameIn
    integer                                     :: obj_type, idx
    integer(HID_T)                              :: gr_idIn, dset_id
    integer(HID_T)                              :: space_id 
    integer                                     :: STAT
    character(StringLength)                     :: NewGroupNameIn
    integer(HID_T)                              :: rank       
    type(string)                                :: nameString, testString
    type(string)                                :: outext, temp_str
    logical                                     :: exitAfterAddVar, testResultsGroup, testTimeGroup
    !Begin-----------------------------------------------------------------
    call h5open_f (STAT)
    if (STAT /= SUCCESS_) then
        outext = 'Failed to use hdf5 lib.'
        call Log%put(outext)
        stop
    endif
    !Get the number of members in the Group
    call h5gn_members_f(IDIn, GroupNameIn, nmembersIn, STAT)
    if (STAT /= SUCCESS_) then
        temp_str = IDIn
        outext = 'Failed to get hdf GroupID. IDIn = ' //temp_str// ' '
        call Log%put(outext)
        outext = 'FileName = ' //trim(self%filename%chars())// ' '
        call Log%put(outext)
        outext = 'GroupName = ' //trim(GroupNameIn)// ' '
        call Log%put(outext)
        temp_str = nmembersIn
        outext = 'nmembersIn = ' //temp_str// '. Stopping'
        call Log%put(outext)
        stop
    endif
        
    do idx = 1, nmembersIn

        call h5gget_obj_info_idx_f(IDIn, GroupNameIn, idx-1, obj_nameIn, obj_type, STAT)
        if (STAT /= SUCCESS_) then
            outext = 'Failed to get GroupName. FileName = ' //trim(self%filename%chars())// ' '
            call Log%put(outext)
            outext = 'GroupName = ' //trim(GroupNameIn)// ' '
            call Log%put(outext)
            outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
            call Log%put(outext)
            stop
        endif
        
        if     (obj_type == H5G_DATASET_F) then

            !Opens data set
            call h5dopen_f      (IDIn, trim(adjustl(obj_nameIn)), dset_id, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to open hdf DataSet. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif
                
            !Opens data space
            call h5dget_space_f (dset_id, space_id, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to open hdf data space. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to get hdf dataset rank. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif
            
            !Check if is a known variable or Dim. If not, skip this dataset and write a warning.
            !This way, only dims and variables are saved into memory
            nameString = trim(obj_nameIn)
            !write(*,*)"nameString ini = ", TRIM(nameString)
            
            !Testing group and obj names
            testResultsGroup = .true.
            testTimeGroup = .true.
            
            testString = trim(GroupNameIn)
            testString = testString%replace(old='/Grid/', new='')
            !write(*,*)"testString = ", TRIM(testString)
            !Because VerticalZ in MOHID is a group name... but the dataset name is Vertical_00001...
            if (Globals%Var%checkVarIsDepth(testString)) then
                testString = testString%replace(old='VerticalZ', new='Vertical')
                testResultsGroup = .false.
            endif
            
            exitAfterAddVar = .false.
            if (testString//"_00001" == trim(obj_nameIn)) then
                nameString = testString
                !write(*,*)"nameString = ", TRIM(nameString)
                exitAfterAddVar = .true.
                testResultsGroup = .false.
                testTimeGroup = .false.
            endif
            
            if (testResultsGroup) then
                testString = trim(GroupNameIn)
                testString = testString%replace(old='/Results/', new='')
                !write(*,*)"testString = ", TRIM(testString)
            
                exitAfterAddVar = .false.
                if (testString//"_00001" == trim(obj_nameIn)) then
                    nameString = testString
                    !write(*,*)"nameString = ", TRIM(nameString)
                    exitAfterAddVar = .true.
                    testTimeGroup = .false.
                endif
            endif
            
            if (testTimeGroup) then
                testString = trim(GroupNameIn)
                testString = testString%replace(old='/', new='')
                !write(*,*)"testString time group = ", TRIM(testString)
                
                if (testString == "GridCorners3D") cycle
            
                exitAfterAddVar = .false.
                if (testString//"_00001" == trim(obj_nameIn)) then
                    nameString = testString
                    !write(*,*)"nameString time group = ", TRIM(nameString)
                    exitAfterAddVar = .true.
                endif
            endif
            
            
            if (.NOT. Globals%Var%checkDimensionName(nameString) .AND. (.NOT. Globals%Var%checkVarSimName(nameString)) .AND. .NOT. exitAfterAddVar) then
                !Skip and let users know
                outext = '[hdf5parser::gethdfNumberOfVarsAndMaxDims]: the hdf file '//trim(self%filename%chars())// ' has a variable that is not recognized by the model: ' // nameString // '. skipping it'
                call Log%put(outext)
                !Closes data space
                call h5sclose_f     (space_id, STAT)
                if (STAT /= SUCCESS_) then
                    outext = 'Failed to close hdf data space = ' //trim(self%filename%chars()) // '. Stopping'
                    call Log%put(outext)
                endif
                cycle
            endif
            
            self%nDims = max(self%nDims, rank+1)
            !write(*,*)"nameString adicionada no addVars = ", TRIM(nameString)
            self%nVars = self%nVars + 1
            
            if (exitAfterAddVar) then
                !Variable already accounted for... dont need to go over the rest of the time instants at this point
                !Closes data space
                call h5sclose_f     (space_id, STAT)
                if (STAT /= SUCCESS_) then
                    outext = 'Failed to close hdf data space = ' //trim(self%filename%chars()) // '. Stopping'
                    call Log%put(outext)
                endif
                exit
            endif

            !Closes data space
            call h5sclose_f     (space_id, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to close hdf data space = ' //trim(self%filename%chars()) // '. Stopping'
                call Log%put(outext)
            endif

        elseif (obj_type ==H5G_GROUP_F) then

            !Looks for futher subgroups
            if (GroupNameIn == "/") then
                NewGroupNameIn = GroupNameIn//trim(adjustl(obj_nameIn))
            else
                NewGroupNameIn = GroupNameIn//"/"//trim(adjustl(obj_nameIn))
            endif

            call h5gopen_f        (IDIn   , trim(adjustl(NewGroupNameIn)), gr_idIn, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to open hdf group = ' //trim(NewGroupNameIn) // '. Stopping'
                call Log%put(outext)
            endif
            
                               
            call self%gethdfNumberOfVarsAndMaxDims (gr_idIn, trim(adjustl(NewGroupNameIn)))
            
            call h5gclose_f       (gr_idIn, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to close hdf group = ' //trim(NewGroupNameIn) // '. Stopping'
                call Log%put(outext)
            endif
                               
        endif
            
    enddo
        
    end subroutine gethdfNumberOfVarsAndMaxDims
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Close the hdf5 file
    !> @param[in] self
    !---------------------------------------------------------------------------
    recursive subroutine hdfReadAllVariables (self, IDIn, GroupNameIn, VarCounter, DimCounter)

    !Arguments-------------------------------------------------------------
    class(hdf5file_class), intent(inout) :: self
    integer(HID_T), intent(in)           :: IDIn
    character(len=*), intent(in)         :: GroupNameIn
    integer, intent(inout)               :: VarCounter, DimCounter
    !Local-----------------------------------------------------------------
    integer                                     :: nmembersIn
    character(StringLength)                     :: obj_nameIn
    integer                                     :: obj_type, idx
    integer(HID_T)                              :: gr_idIn, dset_id
    integer(HID_T)                              :: space_id 
    integer                                     :: STAT
    character(StringLength)                     :: NewGroupNameIn
    integer(HSIZE_T), dimension(:), allocatable :: dims, maxdims
    integer(HID_T)                              :: rank      
    type(string)                                :: outext, temp_str, nameString, testString
    logical                                     :: exitAfterAddVar, testResultsGroup, testTimeGroup, isStaticVariable

    !Begin-----------------------------------------------------------------
    call h5open_f (STAT)
    if (STAT /= SUCCESS_) then
        outext = 'Failed to use hdf5 lib.'
        call Log%put(outext)
        stop
    endif
    
    !Get the number of members in the Group
    call h5gn_members_f(IDIn, GroupNameIn, nmembersIn, STAT)
    if (STAT /= SUCCESS_) then
        temp_str = IDIn
        outext = 'Failed to get hdf GroupID. IDIn = ' //temp_str// ' '
        call Log%put(outext)
        outext = 'FileName = ' //trim(self%filename%chars())// ' '
        call Log%put(outext)
        outext = 'GroupName = ' //trim(GroupNameIn)// ' '
        call Log%put(outext)
        temp_str = nmembersIn
        outext = 'nmembersIn = ' //temp_str// '. Stopping'
        call Log%put(outext)
        stop
    endif
        
    do idx = 1, nmembersIn

        call h5gget_obj_info_idx_f(IDIn, GroupNameIn, idx-1, obj_nameIn, obj_type, STAT)
        if (STAT /= SUCCESS_) then
            outext = 'Failed to get GroupName. FileName = ' //trim(self%filename%chars())// ' '
            call Log%put(outext)
            outext = 'GroupName = ' //trim(GroupNameIn)// ' '
            call Log%put(outext)
            outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
            call Log%put(outext)
            stop
        endif
        
        if     (obj_type == H5G_DATASET_F) then

            !Opens data set
            call h5dopen_f      (IDIn, trim(adjustl(obj_nameIn)), dset_id, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to open hdf DataSet. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif
                

            !Opens data space
            call h5dget_space_f (dset_id, space_id, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to open hdf data space. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif

            !Gets the rank
            call h5sget_simple_extent_ndims_f (space_id, rank, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to get hdf dataset rank. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif
                
            allocate(dims(rank), maxdims(rank))
            
            !Gets the size
            call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, STAT)
            if (STAT < SUCCESS_) then
                outext = 'Failed to get hdf extents. FileName = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// ' '
                call Log%put(outext)
                outext = 'DataSet = ' //trim(obj_nameIn ) // '. Stopping'
                call Log%put(outext)
                stop
            endif
  
            !Check if is a known variable or Dim. If not, skip this dataset and write a warning.
            !This way, only dims and variables are saved into memory
            nameString = trim(obj_nameIn)
            !write(*,*)"nameString ini = ", TRIM(nameString)
            
            !Testing group and obj names
            testResultsGroup = .true.
            testTimeGroup    = .true.
            isStaticVariable = .true.
            
            testString = trim(GroupNameIn)
            testString = testString%replace(old='/Grid/', new='')
            !write(*,*)"testString e  group name in = ", TRIM(testString), trim(GroupNameIn)
            !Because VerticalZ in MOHID is a group name... but the dataset name is Vertical_00001...
            if (Globals%Var%checkVarIsDepth(testString)) then
                testString = testString%replace(old='VerticalZ', new='Vertical')
                testResultsGroup = .false.
            endif
            
            exitAfterAddVar = .false.
            if (testString//"_00001" == trim(obj_nameIn)) then
                nameString = testString
                !write(*,*)"nameString = ", TRIM(nameString)
                exitAfterAddVar = .true.
                testResultsGroup = .false.
                testTimeGroup = .false.
                isStaticVariable = .false.
            endif
            
            if (testResultsGroup) then
                testString = trim(GroupNameIn)
                testString = testString%replace(old='/Results/', new='')
                !write(*,*)"testString = ", TRIM(testString)
            
                exitAfterAddVar = .false.
                if (testString//"_00001" == trim(obj_nameIn)) then
                    nameString = testString
                    !write(*,*)"nameString = ", TRIM(nameString)
                    exitAfterAddVar = .true.
                    testTimeGroup = .false.
                    isStaticVariable = .false.
                endif
            endif
            
            if (testTimeGroup) then
                testString = trim(GroupNameIn)
                testString = testString%replace(old='/', new='')
                !write(*,*)"testString time group = ", TRIM(testString)
            
                if (testString == "GridCorners3D") then
                    deallocate(dims, maxdims)
                    cycle
                endif
                
                exitAfterAddVar = .false.
                if (testString//"_00001" == trim(obj_nameIn)) then
                    nameString = testString
                    !write(*,*)"nameString time group = ", TRIM(nameString)
                    exitAfterAddVar = .true.
                endif
            endif

            if (.NOT. Globals%Var%checkDimensionName(nameString) .AND. (.NOT. Globals%Var%checkVarSimName(nameString)) .AND. .NOT. exitAfterAddVar) then
                !Skip and let users know
                !write(*,*)"nameString desconhecida = ", TRIM(nameString)
                outext = '[hdf5parser::gethdfNumberOfVarsAndMaxDims]: the hdf file '//trim(self%filename%chars())// ' has a variable that is not recognized by the model: ' // nameString // '. skipping it'
                call Log%put(outext)
                !Closes data space
                call h5sclose_f     (space_id, STAT)
                if (STAT /= SUCCESS_) then
                    outext = 'Failed to close hdf data space = ' //trim(self%filename%chars()) // '. Stopping'
                    call Log%put(outext)
                    stop
                endif
                deallocate(dims, maxdims)
                cycle
            endif
            
            self%varData(VarCounter)%name = trim(nameString)
            self%varData(VarCounter)%simName = Globals%Var%getVarSimName(trim(nameString))
            self%varData(VarCounter)%hdf5GroupName = trim(GroupNameIn)
            !write(*,*)"name e simName Var = ", trim(nameString), TRIM(self%varData(VarCounter)%simName)
            self%varData(VarCounter)%varid = VarCounter
                
            self%varData(VarCounter)%units = 'none'
            self%varData(VarCounter)%scale = 1.0
            self%varData(VarCounter)%offset = 0.0
            self%varData(VarCounter)%fillvalue = Globals%Parameters%FillValueReal
                
            if (Globals%Var%checkDimensionName(nameString)) then
                    
                !This is a dimension type variable
                self%varData(VarCounter)%ndims = rank
                allocate(self%varData(VarCounter)%dimids(rank))
                    
                self%dimData(DimCounter)%name    = trim(nameString)
                self%dimData(DimCounter)%simName = Globals%Var%getVarSimName(trim(nameString))
                self%dimData(DimCounter)%nDims = rank
                !write(*,*)"name e simName Dim = ", TRIM(nameString), trim(self%dimData(DimCounter)%simName)
                !Distinguir entre lat e lon para obter o dimLength correto. (1 = Lon, 2 = Lat)
                if (Globals%Var%checkVarIsLon(nameString)) then
                    self%dimData(DimCounter)%length  = maxdims(1)
                    self%dimData(DimCounter)%dimid   = 1
                elseif (Globals%Var%checkVarIsLat(nameString)) then
                    self%dimData(DimCounter)%length  = maxdims(2)
                    self%dimData(DimCounter)%dimid   = 2
                elseif (Globals%Var%checkVarIsDepth(nameString)) then
                    self%dimData(DimCounter)%length  = maxdims(3)
                    self%dimData(DimCounter)%dimid   = 3
                    self%dimData(DimCounter)%nDims = rank + 1
                elseif (Globals%Var%checkVarIsTime(nameString)) then
                    self%dimData(DimCounter)%length  = nmembersIn
                    self%dimData(DimCounter)%dimid   = 4
                endif
                
                DimCounter = DimCounter + 1
            elseif(isStaticVariable) then
                self%varData(VarCounter)%ndims = rank
                allocate(self%varData(VarCounter)%dimids(rank))
                if (rank == 2) then !2D
                    self%varData(VarCounter)%dimids(1) = 1 !Lon
                    self%varData(VarCounter)%dimids(2) = 2 !Lat
                elseif (rank == 3) then
                    self%varData(VarCounter)%dimids(1) = 1 !Lon
                    self%varData(VarCounter)%dimids(2) = 2 !Lat
                    self%varData(VarCounter)%dimids(3) = 3 !Depth (VerticalZ)
                endif
            else
                self%varData(VarCounter)%ndims = rank + 1
                allocate(self%varData(VarCounter)%dimids(rank + 1))
                if (rank == 2) then !2D
                    self%varData(VarCounter)%dimids(1) = 1 !Lon
                    self%varData(VarCounter)%dimids(2) = 2 !Lat
                    self%varData(VarCounter)%dimids(3) = 4 !Time
                elseif (rank == 3) then
                    self%varData(VarCounter)%dimids(1) = 1 !Lon
                    self%varData(VarCounter)%dimids(2) = 2 !Lat
                    self%varData(VarCounter)%dimids(3) = 3 !Depth (VerticalZ)
                    self%varData(VarCounter)%dimids(4) = 4 !Time
                endif
            endif
            
            VarCounter = VarCounter + 1
            
            if (exitAfterAddVar) then
                !Variable already accounted for... dont need to go over the rest of the time instants at this point
                deallocate(dims, maxdims)

                !Closes data space
                call h5sclose_f     (space_id, STAT)
                if (STAT /= SUCCESS_) then
                    outext = 'Failed to close hdf data space = ' //trim(self%filename%chars()) // '. Stopping'
                    call Log%put(outext)
                    stop
                endif
                exit
            endif
            
            deallocate(dims, maxdims)

            !Closes data space
            call h5sclose_f     (space_id, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to close hdf data space = ' //trim(self%filename%chars()) // '. Stopping'
                call Log%put(outext)
                stop
            endif           

        elseif (obj_type ==H5G_GROUP_F) then

            !Looks for futher subgroups
            if (GroupNameIn == "/") then
                NewGroupNameIn = GroupNameIn//trim(adjustl(obj_nameIn))
            else
                NewGroupNameIn = GroupNameIn//"/"//trim(adjustl(obj_nameIn))
            endif

            call h5gopen_f        (IDIn   , trim(adjustl(NewGroupNameIn)), gr_idIn, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to open hdf group. Filename = ' //trim(self%filename%chars())// ' '
                call Log%put(outext)
                outext = 'GroupName = ' //trim(GroupNameIn)// '. Stopping'
                stop
            endif
                               
            call self%hdfReadAllVariables (gr_idIn, trim(adjustl(NewGroupNameIn)), VarCounter, DimCounter)
            call h5gclose_f       (gr_idIn, STAT)
            if (STAT /= SUCCESS_) then
                outext = 'Failed to close hdf group. Group Name = ' //trim(self%varData(VarCounter)%hdf5GroupName%chars())// ' '
                call Log%put(outext)
                stop
            endif
                               
        endif
            
    enddo
    
    end subroutine hdfReadAllVariables
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> reads a variable from the current hdf file
    !> @param[in] self, var, array2D, array3D, outputNumber
    !---------------------------------------------------------------------------
    subroutine readHDFVariable(self, var, array2D, array3D, outputNumber)
    
    !Arguments-------------------------------------------------------------
    class(hdf5file_class), intent(inout)                               :: self
    type(var_t), intent(in)                                            :: var
    real(prec), allocatable, dimension(:,:), intent(inout), optional   :: array2D
    real(prec), allocatable, dimension(:,:,:), intent(inout), optional :: array3D
    integer, intent(in), optional                                      :: outputNumber
    !Local-----------------------------------------------------------------
    integer(HID_T)                                                     :: dset_id, gr_id
    integer(HSIZE_T), dimension(7)                                     :: dims
    integer(HID_T)                                                     :: NumType
    type(string)                                                       :: outext
    character(StringLength)                                            :: AuxChar
    integer(HID_T)                                                     :: STAT_CALL
    !Begin------------------------------------------------------------------
    
    NumType = H5T_NATIVE_DOUBLE

    !Opens the Group
    call h5gopen_f (self%hdf5ID, trim(var%hdf5GroupName%chars()), gr_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_)then
        outext = 'Failed to open hdf group = ' //trim(var%hdf5GroupName) // '. Stopping'
        call Log%put(outext)
        stop
    endif

    !Opens the DataSet
    if (present(outputNumber)) then
        call ConstructDSName (var%name%chars(), outputNumber, AuxChar)
    else
        AuxChar = var%name%chars()
    endif

    call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to open hdf DataSet. FileName = ' //trim(self%filename%chars())// ' '
        call Log%put(outext)
        outext = 'var name = ' //trim(var%name%chars())// ' '
        call Log%put(outext)
        stop
    endif

    !Read the data to the file
    if (PRESENT(array2D)) then
        dims(1) = SIZE(array2D,1)
        dims(2) = SIZE(array2D,2)
        !dims(3) = Me%Limits%KUB - Me%Limits%KLB + 1
        call h5dread_f   (dset_id, NumType, array2D, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to read hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif
    elseif (PRESENT(array3D)) then
        dims(1) = SIZE(array3D,1)
        dims(2) = SIZE(array3D,2)
        dims(3) = SIZE(array3D,3)
        call h5dread_f   (dset_id, NumType, array3D, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to read hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif
    endif
                             
    !End access to the dataset
    call h5dclose_f  (dset_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to close hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
        call Log%put(outext)
        stop
    endif

    !Closes group
    call h5gclose_f( gr_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to close hdf group. Group Name = ' //trim(var%hdf5GroupName%chars())// ' '
        call Log%put(outext)
        stop
    endif
    
    end subroutine readHDFVariable
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> reads a variable from the current hdf file
    !> @param[in] self, var, array2D, array3D, outputNumber
    !---------------------------------------------------------------------------
    subroutine readHDFIntVariable(self, var, array2D, array3D, outputNumber)
    
    !Arguments-------------------------------------------------------------
    class(hdf5file_class), intent(inout)                               :: self
    type(var_t), intent(in)                                            :: var
    integer, allocatable, dimension(:,:), intent(inout), optional      :: array2D
    integer, allocatable, dimension(:,:,:), intent(inout), optional    :: array3D
    integer, intent(in), optional                                      :: outputNumber
    !Local-----------------------------------------------------------------
    integer(HID_T)                                                     :: dset_id, gr_id
    integer(HSIZE_T), dimension(7)                                     :: dims
    integer(HID_T)                                                     :: NumType
    type(string)                                                       :: outext
    character(StringLength)                                            :: AuxChar
    integer(HID_T)                                                     :: STAT_CALL
    !Begin------------------------------------------------------------------
    
    NumType = H5T_NATIVE_DOUBLE

    !Opens the Group
    call h5gopen_f (self%hdf5ID, trim(var%hdf5GroupName%chars()), gr_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_)then
        outext = 'Failed to open hdf group = ' //trim(var%hdf5GroupName) // '. Stopping'
        call Log%put(outext)
        stop
    endif

    !Opens the DataSet
    if (present(outputNumber)) then
        call ConstructDSName (var%name%chars(), outputNumber, AuxChar)
    else
        AuxChar = var%name%chars()
    endif

    call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to open hdf DataSet. FileName = ' //trim(self%filename%chars())// ' '
        call Log%put(outext)
        outext = 'var name = ' //trim(var%name%chars())// ' '
        call Log%put(outext)
        stop
    endif

    !Read the data to the file
    if (PRESENT(array2D)) then
        dims(1) = SIZE(array2D,1)
        dims(2) = SIZE(array2D,2)
        call h5dread_f   (dset_id, NumType, array2D, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to read hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif
    elseif (PRESENT(array3D)) then
        dims(1) = SIZE(array3D,1)
        dims(2) = SIZE(array3D,2)
        dims(3) = SIZE(array3D,3)
        call h5dread_f   (dset_id, NumType, array3D, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to read hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif
    endif
                             
    !End access to the dataset
    call h5dclose_f  (dset_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to close hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
        call Log%put(outext)
        stop
    endif

    !Closes group
    call h5gclose_f( gr_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to close hdf group. Group Name = ' //trim(var%hdf5GroupName%chars())// ' '
        call Log%put(outext)
        stop
    endif
    
    end subroutine readHDFIntVariable
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> reads time variable from the current hdf file
    !> @param[in] self, var, array1D, instances
    !---------------------------------------------------------------------------
    subroutine readHDFTime(self, var, array1D, instances)
    
    !Arguments-------------------------------------------------------------
    class(hdf5file_class), intent(inout)                               :: self
    type(var_t), intent(in)                                            :: var
    real(prec), allocatable, dimension(:), intent(inout)               :: array1D
    integer, intent(in), optional                                      :: instances
    !Local-----------------------------------------------------------------
    real(prec), dimension(6)                                           :: timeVector
    integer(HID_T)                                                     :: dset_id, gr_id
    integer(HSIZE_T), dimension(7)                                     :: dims
    integer(HID_T)                                                     :: NumType
    type(string)                                                       :: outext
    character(StringLength)                                            :: AuxChar
    integer(HID_T)                                                     :: STAT_CALL
    real(prec)                                                         :: timeInSeconds
    integer                                                            :: i
    !Begin------------------------------------------------------------------
    
    NumType = H5T_NATIVE_DOUBLE

    dims(1) = 6
    !Opens the Group
    call h5gopen_f (self%hdf5ID, trim(var%hdf5GroupName%chars()), gr_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_)then
        outext = 'Failed to open hdf group = ' //trim(var%hdf5GroupName) // '. Stopping'
        call Log%put(outext)
        stop
    endif

    do i=1, instances
        !Opens the DataSet
        call ConstructDSName (var%name%chars(), i, AuxChar)

        call h5dopen_f (gr_id, trim(adjustl(AuxChar)), dset_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to open hdf DataSet. FileName = ' //trim(self%filename%chars())// ' '
            call Log%put(outext)
            outext = 'var name = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif

        !Read the data to the file
        call h5dread_f   (dset_id, NumType, timeVector, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to read hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif
        timeInSeconds = correctHdf5Time(timeVector)
        array1D(i) = timeInSeconds  
        !End access to the dataset
        call h5dclose_f  (dset_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            outext = 'Failed to close hdf DataSet. VarName = ' //trim(var%name%chars())// ' '
            call Log%put(outext)
            stop
        endif
    enddo
    !Closes group
    call h5gclose_f( gr_id, STAT_CALL)
    if (STAT_CALL /= SUCCESS_) then
        outext = 'Failed to close hdf group. Group Name = ' //trim(var%hdf5GroupName%chars())// ' '
        call Log%put(outext)
        stop
    endif
    
    end subroutine readHDFTime
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> constructs the dataset name according to the time instant. copied from MOHID hdf5 module
    !> @param[in] Name, OutputNumber, AuxChar
    subroutine ConstructDSName (Name, OutputNumber, AuxChar)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Name
        integer                                     :: OutputNumber
        character(len=*)                            :: AuxChar

        !Local-----------------------------------------------------------------
        character(StringLength)                     :: AuxNum

        write(AuxNum, fmt=*)OutputNumber

        if     (OutputNumber < 10     ) then
            AuxChar = trim(adjustl(Name))//"_0000"//trim(adjustl(AuxNum))
        elseif (OutputNumber < 100    ) then
            AuxChar = trim(adjustl(Name))//"_000" //trim(adjustl(AuxNum))
        elseif (OutputNumber < 1000   ) then
            AuxChar = trim(adjustl(Name))//"_00"  //trim(adjustl(AuxNum))
        elseif (OutputNumber < 10000  ) then
            AuxChar = trim(adjustl(Name))//"_0"   //trim(adjustl(AuxNum))
        else
            AuxChar = trim(adjustl(Name))//"_"    //trim(adjustl(AuxNum))
        endif

    end subroutine ConstructDSName
    
    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> Close the hdf5 file
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine closeFile(self)
    class(hdf5file_class),intent(inout) :: self
    integer(HID_T)                      :: STAT_
    type(string)                        :: outext
    !Begin----------------------------------------------------------------------
    call h5fclose_f(self%hdf5ID,   HDFERR = STAT_)
    if (STAT_ /= SUCCESS_) then
        outext = 'Failed to get GroupName. FileName = ' //trim(self%filename%chars())// ' '
        call Log%put(outext)
        stop
    endif
    
    end subroutine closeFile

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Debug the hdf5 error after a hdf5 command
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine check(self)
    class(hdf5file_class), intent(inout) :: self
    type(string) :: outext
    if(self%status /= NF90_noerr) then
        outext = '[hdf5parser::check]: '//trim(NF90_strerror(self%status))//', stoping'
        call Log%put(outext)
        stop
    end if
    end subroutine check

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> print variable metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printVarsHDF5(self)
    class(var_t), intent(inout) :: self
    type(string) :: outext, temp
    integer :: i
    outext = '--->Variable: '// self%name //new_line('a')
    temp = self%varid
    outext = outext//'       varid = '//temp//new_line('a')
    outext = outext//'       units = '//self%units//new_line('a')
    temp = self%dimids(1)
    outext = outext//'       dimids = '//temp
    do i=2, size(self%dimids)
        temp = self%dimids(i)
        outext = outext//', '//temp
    end do
    call Log%put(outext,.false.)
    end subroutine printVarsHDF5

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> print dimensions metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printDimsHDF5(self)
    class(dim_t), intent(inout) :: self
    type(string) :: outext, temp
    outext = '--->Dimension: '// self%name //new_line('a')
    temp = self%dimid
    outext = outext//'       dimid = '//temp//new_line('a')
    temp = self%varid
    outext = outext//'       varid = '//temp//new_line('a')
    outext = outext//'       units = '//self%units//new_line('a')
    temp = self%length
    outext = outext//'       lenght = '//temp
    call Log%put(outext,.false.)
    end subroutine printDimsHDF5

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho
    !> @brief
    !> corrects the time array to a more efficient format if needed.
    !> @param[in] timeComments, timeArray
    !---------------------------------------------------------------------------
    real(prec) function correcthdf5Time(timeArray)
    real(prec), dimension(:), intent(in) :: timeArray
    integer, dimension(6) :: date
    type(timedelta) :: dateOffset
    real(prec) :: offset
    type(string) :: isoDateStr
    type(datetime) :: hdf5Date

    isoDateStr = timeArray(1)//' '//timeArray(2)//' '//timeArray(3)//' '//timeArray(4)//' '//timeArray(5)//' '//timeArray(6)
    date = Utils%getDateFromISOString(isoDateStr)
    hdf5Date = Utils%getDateTimeFromDate(date)
    dateOffset = Globals%SimTime%StartDate - hdf5Date
    offset = -dateOffset%total_seconds()
    
    correcthdf5Time = offset
    end function correcthdf5Time


    end module hdf5Parser_mod