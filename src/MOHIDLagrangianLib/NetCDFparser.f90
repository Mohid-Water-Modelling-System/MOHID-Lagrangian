    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : NetCDFparser
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2019
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa Paz
    !
    ! DESCRIPTION:
    !> Module that defines a netcdf file model class, responsible for abstracting
    !> the parsing of netcdf files. Each object of this class is responsible for
    !> one file, effectivelly representing it.
    !> The object should return the necessary data fields with corresponding meta
    !> data, such as names and units. An variable library is consulted when
    !> required, because even NetCDF CF compliance allows ambiguity.
    !------------------------------------------------------------------------------

    module netcdfParser_mod

#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif

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
        logical :: reverse_axis = .false., negate = .false.
        logical :: reverse_data = .false.
    contains
    procedure :: print => printDimsNC
    end type dim_t

    type :: var_t
        type(string) :: name
        type(string) :: simName
        integer :: varid
        type (string) :: units
        integer :: ndims
        logical :: isDimensionVar
        integer, allocatable, dimension(:) :: dimids
        integer :: natts
        real(prec) :: offset, scale, fillvalue
    contains
    procedure :: print => printVarsNC
    end type var_t

    type :: ncfile_class !< A class that models a netcdf file
        type(string) :: filename   !< name of the file to read
        integer :: ncID             !< ID of the file
        integer :: nDims, nVars, nAtt, uDimID   !< number of dimensions, variables, attributes and dim IDs on the file
        type(dim_t), allocatable, dimension(:) :: dimData !< metadata from the dimensions on the file
        type(var_t), allocatable, dimension(:) :: varData   !< metadata from the variables on the file
        integer :: mDims
        logical :: grid_2D
        integer :: status
    contains
    procedure :: initialize => getFile
    procedure :: getVarDimensions_variable
    procedure :: getVar
    procedure :: finalize => closeFile
    procedure, private :: check
    procedure, private :: getNCid
    procedure, private :: getNCglobalMetadata
    procedure, private :: getNCDimMetadata
    procedure, private :: getNCVarMetadata
    procedure, private :: getDimByDimID
    procedure, private :: check2dGrid
    procedure :: print => printNcInfo
    end type ncfile_class
    
    type :: ncReader_class
    contains
    procedure :: getFullFile
    end type ncReader_class

    public :: ncfile_class, ncReader_class

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> instantiates and returns a background object with the data from a NC file
    !> @param[in] self, fileName, varList, syntecticVar
    !---------------------------------------------------------------------------
    type(background_class) function getFullFile(self, fileName, varList, syntecticVar)
    class(ncReader_class), intent(in) :: self
    type(string), intent(in) :: fileName
    type(string), dimension(:), intent(in) :: varList
    logical, dimension(:), intent(in) :: syntecticVar    
    type(ncfile_class) :: ncFile
    type(generic_field_class), allocatable, dimension(:) :: backgrounDims
    type(generic_field_class), allocatable, dimension(:) :: gfield
    type(string) :: name, units
    type(box) :: extents
    type(vector) :: pt
    real(prec), dimension(3,2) :: dimExtents
    integer :: i, realVarIdx
    type(string) :: outext
    real(prec) :: valorminimo

    allocate(gfield(size(syntecticVar)))
    realVarIdx = 0
    units = '-'
    outext = '->Reading '//fileName
    call Log%put(outext,.false.)

    call ncFile%initialize(fileName)
    do i=1, size(syntecticVar)
        if(.not.syntecticVar(i)) then !finding the first real variable to extract dimension arrays
            call ncFile%getVarDimensions_variable(varList(i), backgrounDims)
            !call ncFile%getVarDimensions(varList(i), backgrounDims)
            
            !write(*,*)"backgroundDims 1, 2, 3, 4 : ", backgrounDims(1)%Name, backgrounDims(2)%Name, backgrounDims(3)%Name
            if (allocated(backgrounDims)) then
                realVarIdx = i
                valorminimo = backgrounDims(3)%GetFieldMinBound()
                exit
            end if 
        end if
    end do
    if (realVarIdx /= 0) then
        do i=1, size(syntecticVar)
            !write(*,*)"Getting Var name = ", varList(i)
            if(.not.syntecticVar(i)) then !normal variable, put it on a generic field
                !write(*,*)"getting normal variable"
                call ncFile%getVar(varList(i), gfield(i))
            else                          !synthetic variable to be constructed based on the field of a normal variable
                call ncFile%getVar(varList(realVarIdx), gfield(i), .true., varList(i), units)
            end if
        end do
    end if
    
    call ncFile%finalize()
    
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
    valorminimo = backgrounDims(3)%getFieldMinBound()
    !write(*,*) "Sai do getFullFile"
    end function getFullFile

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Parse the netcdf file headers and assembles the file model
    !> @param[in] self, filename
    !---------------------------------------------------------------------------
    subroutine getFile(self, filename)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: filename
    self%filename = filename
    call self%getNCid()
    call self%getNCglobalMetadata()
    call self%getNCVarMetadata()
    call self%check2dGrid()
    call self%getNCDimMetadata()
    end subroutine getFile

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Opens the netcdf file, assigns the NCID and checks for errors
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNcid(self)
    class(ncfile_class), intent(inout) :: self
    self%status = NF90_open(trim(self%filename%chars()), NF90_NOWRITE, self%ncID)
    call self%check()
    end subroutine getNcid

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Inquires the nc file for global metadata - number  of dimensions,
    !> variables and attributes
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNCglobalMetadata(self)
    class(ncfile_class), intent(inout) :: self
    self%status = nf90_inquire(self%ncID, self%nDims, self%nVars, self%nAtt, self%uDimID)
    !write(*,*)"entrada getNCglobalMetadata"
    !write(*,*)"number dimensions = ", self%nDims
    call self%check()
    end subroutine getNCglobalMetadata

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Inquires the nc file for variable metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNCVarMetadata(self)
    class(ncfile_class), intent(inout) :: self
    integer :: i, ndims, nAtts, maxdims
    integer :: dimids(self%nDims)
    integer :: tempStatus
    character(CHAR_LEN) :: varName, units
    allocate(self%varData(self%nVars))
    maxdims=0
    do i=1, self%nVars
        self%status = nf90_inquire_variable(self%ncID, i, varName, ndims=ndims, dimids=dimids, nAtts=nAtts)
        call self%check()
        self%varData(i)%name = trim(varName)
        self%varData(i)%simName = Globals%Var%getVarSimName(self%varData(i)%name)
        self%varData(i)%varid = i
        self%varData(i)%ndims = ndims
        maxdims = max(maxdims, ndims)
        allocate(self%varData(i)%dimids(ndims))
        self%varData(i)%dimids = dimids(1:ndims)
        self%varData(i)%nAtts = nAtts
        tempStatus = nf90_get_att(self%ncID, i, 'units', units)
        if (tempStatus == -43) units = "not set"
        self%varData(i)%units = trim(units)
        tempStatus = nf90_get_att(self%ncID, i, "scale_factor", self%varData(i)%scale)
        if (tempStatus == -43) self%varData(i)%scale = 1.0
        tempStatus = nf90_get_att(self%ncID, i, "add_offset", self%varData(i)%offset)
        if (tempStatus == -43) self%varData(i)%offset = 0.0
        tempStatus = nf90_get_att(self%ncID, i, "_FillValue", self%varData(i)%fillvalue)
        if (tempStatus == -43) self%varData(i)%fillvalue = MV
    end do
    !self%mDims = maxdims
    end subroutine getNCVarMetadata

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Inquires the nc file for 2D grid coordinates
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine check2dGrid(self)
    class(ncfile_class), intent(inout) :: self
    integer :: i
    self%grid_2D = .false.
    do i=1, self%nVars
        if ((self%varData(i)%ndims == 2) .and. (Globals%Var%checkDimensionName(self%varData(i)%name))) then
            !Current variable is 2D and is a dimension variable
            self%grid_2D = .true.
        end if
    end do
    end subroutine check2dGrid
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC and Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Inquires the nc file for dimension metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNCDimMetadata(self)
    class(ncfile_class), intent(inout) :: self
    integer :: i, j, k, dimLength
    character(CHAR_LEN) :: dimName
    logical :: dimNameIsValid
    type(string) :: outext
    !Begin -----------------------------------------
    if (self%grid_2D) then
        !probably Maretec's netcdf containing 2D lat and lon
        !use max dims which is the maximum number of dimensions of the variables inside the nc.
        if (self%mDims < 5 .and. self%mDims > 2) then
            !3D + time
            !write(*,*)"entrada getNCDimMetadata"
            !write(*,*)"number dimensions = ", self%nDims
            !write(*,*)"number vars = ", self%nVars
            !data is 3D or 4D
            allocate(self%dimData(self%mDims))
            !k goes through max dimensions (3 or 4)
            k = 1
            !i goes through all file dimensions (in MOHID 6 are provided) so need to exclude those not wanted.
            do i = 1, self%nDims
                if (k > self%mDims) exit
                self%status = nf90_inquire_dimension(self%ncID, i, dimName, dimLength)
                call self%check()
                self%dimData(k)%name = trim(dimName)
                !write(*,*)"current dimension variable in file = ", trim(self%dimData(k)%name)
                dimNameIsValid = Globals%Var%checkDimensionName(self%dimData(k)%name)
                if (dimNameIsValid) then
                    !allocate new dimension to the dimData place holder
                    self%dimData(k)%simName = Globals%Var%getVarSimName(self%dimData(k)%name)
                    !write(*,*)"Sim name of variable in file = ", trim(self%dimData(k)%simName)
                    self%dimData(k)%length = dimLength
                    self%status = nf90_inq_dimid(self%ncID, self%dimData(k)%name%chars(), self%dimData(k)%dimid)
                    !write(*,*)"Dim data dimID = ", self%dimData(k)%dimid
                    call self%check()
                    do j=1, self%nVars
                        !write(*,*)"dim name = ", self%dimData(k)%name
                        !write(*,*)"vardata name = ", self%varData(j)%name
                        if (self%dimData(k)%name == self%varData(j)%name) then
                            self%dimData(k)%units = self%varData(j)%units
                            self%dimData(k)%varid = self%varData(j)%varid
                            self%varData(j)%isDimensionVar = dimNameIsValid
                            !write(*,*)"dim data ID = ", self%dimData(i)%varid
                            exit
                        end if
                    end do
                    k = k + 1
                end if
            end do
            
            if (k < self%mDims) then
                outext = '[NetCDFparser::getNCDimMetadata]: the 4D file '//trim(self%filename%chars())// ' has less thant 4 dimension variables. Stopping'
                call Log%put(outext)
                stop
            end if
        else
            outext = '[NetCDFparser::getNCDimMetadata]: ncfile '//trim(self%filename%chars())//' has less than 3 or more than 4 dimensions. Stopping'
            call Log%put(outext)
            stop
        end if
    
    else
        allocate(self%dimData(self%nDims))
        self%mDims = self%nDims
        !write(*,*)"entrada getNCDimMetadata"
        !write(*,*)"number dimensions = ", self%nDims
        !write(*,*)"number vars = ", self%nVars
        do i=1, self%nDims
            self%status = nf90_inquire_dimension(self%ncID, i, dimName, dimLength)
            call self%check()
            self%dimData(i)%name = trim(dimName)
            !write(*,*)"current dimension variable in file = ", self%dimData(i)%name
            self%dimData(i)%simName = Globals%Var%getVarSimName(self%dimData(i)%name)
            !write(*,*)"Sim name of variable in file = ", self%dimData(i)%simName
            self%dimData(i)%length = dimLength
            self%status = nf90_inq_dimid(self%ncID, self%dimData(i)%name%chars(), self%dimData(i)%dimid)
            !write(*,*)"Dim data dimID = ", self%dimData(i)%dimid
            call self%check()
            do j=1, self%nVars
                !write(*,*)"dim name = ", self%dimData(i)%name
                !write(*,*)"vardata name = ", self%varData(j)%name
                if (self%dimData(i)%name == self%varData(j)%name) then
                    self%dimData(i)%units = self%varData(j)%units
                    self%dimData(i)%varid = self%varData(j)%varid
                    !write(*,*)"dim data ID = ", self%dimData(i)%varid
                    exit
                end if
            end do
        end do    
    end if

    !write(*,*)"saida getNCDimMetadata"
    end subroutine getNCDimMetadata

    !---------------------------------------------------------------------------
    !> @author Joao Sobrinho - Colab Atlantic
    !> @brief
    !> Reads the dimension fields from the nc file for a given variable.
    !> returns an array of scalar 1D fields, each with a name, units and data (includes Lat, Lon and Depth)
    !> @param[in] self, varName, dimsArrays
    !---------------------------------------------------------------------------
    subroutine getVarDimensions_variable(self, varName, dimsArrays)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(generic_field_class), allocatable, dimension(:), intent(out) :: dimsArrays
    real(prec), allocatable, dimension(:) :: tempRealArray1D
    real(prec), allocatable, dimension(:,:) :: tempRealArray2D, auxRealField2D
    type(string) :: dimName, dimUnits
    integer, allocatable, dimension(:) :: varShape
    integer :: i, j, k, ndim, var, i1, j1
    type(dim_t) :: tempDim
    logical :: increase_flag, neg_flag, foundDimVar, found2dDimVar, LatOrLon
    !Begin-----------------------------------------------------------------------
    !write(*,*)"entrei getVarDimensions_variable"
    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName) then   !found the requested var
            !write(*,*)"var name = ", varName
            !write(*,*)"number of dims in varname = ", self%varData(i)%ndims
            allocate(dimsArrays(self%varData(i)%ndims)) !allocating dimension fields
            allocate(varShape(self%varData(i)%ndims))
            !write(*,*)"allocate done"
            !Get vector dimension vectors
            do ndim =1, self%varData(i)%ndims   !going trough all of the variable dimensions
                !write(*,*)"ndim =  ", ndim
                tempDim = self%getDimByDimID(self%varData(i)%dimids(ndim))
                varShape(ndim) = tempDim%length
            end do
            !write(*,*)"varshape done"
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                !write(*,*)"dims dimension = ", self%nDims
                foundDimVar = .false.
                do k=1, self%mDims  !going trough the dimensions of the file
                    !write(*,*)"enteri novo k "
                    !write(*,*)"dim ID da dimensao atual da variavel =  ", self%varData(i)%dimids(j)
                    !write(*,*)"dim ID do dimData =  ", self%dimData(k)%dimid
                    !write(*,*)"nome da dimensao (k) =  ", trim(self%dimData(k)%name)
                    if (self%varData(i)%dimids(j) == self%dimData(k)%dimid) then    !found a corresponding dimension between the variable and the file
                        found2dDimVar = .false.
                        !write(*,*)"enteri if self%varData(i)%dimids(j) == self%dimData(k)%dimid "
                        LatOrLon = Globals%Var%checkLatOrLon(self%dimData(k)%name) !Lat and lon have the same dimensions so dont need to check if it is lat or lon
                cd1:    if (self%grid_2D .and. LatOrLon) then
                    dvar:   do var = 1, self%nVars !going trough all file variables
                                LatOrLon = .false.
                                LatOrLon = Globals%Var%checkLatOrLon(self%varData(var)%name)
                                if ((LatOrLon) .and. (Globals%Var%getVarSimName(self%dimData(k)%name)==self%varData(var)%name)) then
                                    !Lat or Lon found, so allocate 2D field
                                    !write(*,*)"varShape(1) = ", varShape(1)
                                    !write(*,*)"varShape(2) = ", varShape(2)
                                    allocate(auxRealField2D(varShape(1), varShape(2))) !allocating a place to read the field data to
                                    allocate(tempRealArray2D(varShape(2), varShape(1))) !allocating a place to read the field data to
                                    dimName = self%varData(var)%simName
                                    !write(*,*)"simulation name of dimension = ", dimName
                                    dimUnits = self%varData(var)%units
                                    !write(*,*)"dimUnits of dimension = ", dimUnits
                                    !Here is where Lat, Lon and Depth are read from the netcdf
                                    self%status = nf90_get_var(self%ncID, self%varData(var)%varid, auxRealField2D)
                                    
                                    do i1=1,varShape(1)
                                    do j1=1,varShape(2)
                                        tempRealArray2D(j1,i1) = auxRealField2D(i1,j1)
                                    enddo
                                    enddo
                                    call self%check()
                                    found2dDimVar = .true.
                                    exit cd1
                                end if
                            end do dvar
                        end if cd1
                        
                        !write(*,*)"found2dDimVar = ", found2dDimVar
                        !for 1D time or depth
                        if (.not. found2dDimVar) then
                            allocate(tempRealArray1D(self%dimData(k)%length)) !allocating a place to read the field data to
                            dimName = self%dimData(k)%simName
                            !write(*,*)"Dim name = ", dimName
                            dimUnits = self%dimData(k)%units
                            !write(*,*)"dimUnits = ", dimUnits
                            !write(*,*)"self%dimData(k)%varid = ", self%dimData(k)%varid
                            !Here is where Lat, Lon and Depth are read from the netcdf
                            self%status = nf90_get_var(self%ncID, self%dimData(k)%varid, tempRealArray1D)
                            call self%check()
                        
                            !need to check for 'level' variable specific issues
    
                            !@Daniel
                            ! To interpolate on regular meshes, the search cell method uses the lower bound data and the 1 index
                            ! as a reference to locate the point inside the data cell in time, longitude and latitude axis.
                            ! To be consistent, the cell search in depth dimension must be done in the same way. It means:
                            ! The lower bound in depth axis (maximum depth) must be in the 1 index (lower bound) and growing 
                            ! up till the surface (minimum depth, upper bound) .
                            ! SO:
                            ! 0) index = (1,2,3,.....n),
                            !    axis = (-bottom,..., +surface)
                            ! To check this and adjust the data to this criteria, we need to check the two following conditions
                            if (dimName == Globals%Var%level) then
                                !1) The depth must increase. If it does not increase, must be reversed.
                                !2) The axis should be negative. If it is not negative, negate it.
                                increase_flag = all(tempRealArray1D(2:) >= tempRealArray1D(:size(tempRealArray1D)-1)) 
                                neg_flag = all(tempRealArray1D(:) <= 0)
    
                                if ((increase_flag .eqv. .true.) .and. (neg_flag .eqv. .true.))  then
                                    self%dimData(k)%reverse_data = .false.
                                    self%dimData(k)%reverse_axis = .false.
                                    self%dimData(k)%negate = .false.
                                else if ((increase_flag .eqv. .true.) .and. (neg_flag .eqv. .false.))  then
                                    self%dimData(k)%reverse_data = .true.
                                    self%dimData(k)%reverse_axis = .true.
                                    self%dimData(k)%negate = .true.
                                else if ((increase_flag .eqv. .false.) .and. (neg_flag .eqv. .false.)) then
                                    self%dimData(k)%reverse_data = .false.
                                    self%dimData(k)%reverse_axis = .false.
                                    self%dimData(k)%negate = .true.
                                else if ((increase_flag .eqv. .false.) .and. (neg_flag .eqv. .true.)) then
                                    self%dimData(k)%reverse_data = .true.
                                    self%dimData(k)%reverse_axis = .true.
                                    self%dimData(k)%negate = .false.
                                end if 
    
                                if (self%dimData(k)%reverse_axis .eqv. .true.) then
                                    tempRealArray1D = tempRealArray1D(size(tempRealArray1D):1:-1)
                                end if
    
                                if (self%dimData(k)%negate .eqv. .true.) then
                                    tempRealArray1D = -tempRealArray1D
                                end if
                        
                            end if
                            !need to check for 'time' variable specific issues
                            if (dimName == Globals%Var%time) then
                                call correctNCTime(dimUnits, tempRealArray1D)
                            end if
                        end if
                        
                        if (allocated(tempRealArray1D)) then
                            if (dimName == Globals%Var%lon) then
                                !write(*,*)"initializing dimArrays(2) with dimName = ", dimName
                                call dimsArrays(2)%initialize(dimName, dimUnits, tempRealArray1D)
                                !write(*,*)"dimData k = ", trim(Globals%Var%getVarSimName(self%dimData(k)%name))
                                dimsArrays(2)%name = Globals%Var%getVarSimName(self%dimData(k)%name)
                                !write(*,*)"nome gravado = ", j, trim(dimsArrays(j)%name)
                            elseif (dimName == Globals%Var%lat) then
                                !write(*,*)"initializing dimArrays(1) with dimName = ", dimName
                                call dimsArrays(1)%initialize(dimName, dimUnits, tempRealArray1D)
                                !write(*,*)"dimData k = ", trim(Globals%Var%getVarSimName(self%dimData(k)%name))
                                dimsArrays(2)%name = Globals%Var%getVarSimName(self%dimData(k)%name)
                                !write(*,*)"nome gravado = ", j, trim(dimsArrays(j)%name)
                            else
                                call dimsArrays(j)%initialize(dimName, dimUnits, tempRealArray1D)
                                dimsArrays(j)%name = Globals%Var%getVarSimName(self%dimData(k)%name)
                            endif
                            
                            deallocate(tempRealArray1D)
                            foundDimVar = .true.
                        elseif (allocated(tempRealArray2D)) then
                            call dimsArrays(j)%initialize(dimName, dimUnits, tempRealArray2D)
                            deallocate(tempRealArray2D)
                            foundDimVar = .true.
                        end if
                        !write(*,*)"dimData k = ", trim(Globals%Var%getVarSimName(self%dimData(k)%name))
                        !dimsArrays(j)%name = Globals%Var%getVarSimName(self%dimData(k)%name)
                        !write(*,*)"nome gravado = ", j, trim(dimsArrays(j)%name)
                    end if
                end do
            end do
        end if
    end do
    !write(*,*)"sai getVarDimensions_variable"
    end subroutine getVarDimensions_variable
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Reads the dimension fields from the nc file for a given variable.
    !> returns an array of scalar 1D fields, each with a name, units and data (includes Lat, Lon and Depth)
    !> @param[in] self, varName, dimsArrays
    !---------------------------------------------------------------------------
    subroutine getVarDimensions(self, varName, dimsArrays)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(scalar1d_field_class), allocatable, dimension(:), intent(out) :: dimsArrays
    real(prec), allocatable, dimension(:) :: tempRealArray, tempRealArrayDelta
    type(string) :: dimName, dimUnits
    integer :: i, j, k
    logical :: increase_flag, neg_flag 
    !Variable IDs (self%dimData(k)) : 1-time, 2-lat, 3-lon, 4-depth
    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName) then   !found the requested var
            !write(*,*)"dims var name = ", varName
            !write(*,*)"number of dims in vardata = ", self%varData(i)%ndims
            !write(*,*)"number of dimensions of file = ", self%nDims
            allocate(dimsArrays(self%varData(i)%ndims)) !allocating output fields
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                !write(*,*)"j = ", j
                do k=1, self%nDims  !going trough all available dimensions of the file
                    !write(*,*)"Dimension ID of current variable = ", self%varData(i)%dimids(j)
                    !write(*,*)"Current dimension ID = ", self%dimData(k)%dimid
                    !write(*,*)"Dim name = ", self%dimData(k)%name
                    if (self%varData(i)%dimids(j) == self%dimData(k)%dimid) then    !found a corresponding dimension between the variable and the file
                        allocate(tempRealArray(self%dimData(k)%length)) !allocating a place to read the field data to
                        !write(*,*)"Dim length = ", self%dimData(k)%length
                        dimName = self%dimData(k)%simName
                        dimUnits = self%dimData(k)%units
                        !write(*,*)"dimUnits = ", dimUnits
                       ! write(*,*)"variable ID of dimension = ", self%dimData(k)%varid
                        !Here is where Lat, Lon and Depth are read from the netcdf
                        self%status = nf90_get_var(self%ncID, self%dimData(k)%varid, tempRealArray)
                        call self%check()
                        !need to check for 'level' variable specific issues
    
                        !@Daniel
                        ! To interpolate on regular meshes, the search cell method uses the lower bound data and the 1 index
                        ! as a reference to locate the point inside the data cell in time, longitude and latitude axis.
                        ! To be consistent, the cell search in depth dimension must be done in the same way. It means:
                        ! The lower bound in depth axis (maximum depth) must be in the 1 index (lower bound) and growing 
                        ! up till the surface (minimum depth, upper bound) .
                        ! SO:
                        ! 0) index = (1,2,3,.....n),
                        !    axis = (-bottom,..., +surface)
                        ! To check this and adjust the data to this criteria, we need to check the two following conditions
                        if (dimName == Globals%Var%level) then
                            !1) The depth must increase. If it does not increase, must be reversed.
                            !2) The axis should be negative. If it is not negative, negate it.
                            increase_flag = all(tempRealArray(2:) >= tempRealArray(:size(tempRealArray)-1)) 
                            neg_flag = all(tempRealArray(:) <= 0)
    
                            if ((increase_flag .eqv. .true.) .and. (neg_flag .eqv. .true.))  then
                                self%dimData(k)%reverse_data = .false.
                                self%dimData(k)%reverse_axis = .false.
                                self%dimData(k)%negate = .false.
                            else if ((increase_flag .eqv. .true.) .and. (neg_flag .eqv. .false.))  then
                                self%dimData(k)%reverse_data = .true.
                                self%dimData(k)%reverse_axis = .true.
                                self%dimData(k)%negate = .true.
                            else if ((increase_flag .eqv. .false.) .and. (neg_flag .eqv. .false.)) then
                                self%dimData(k)%reverse_data = .false.
                                self%dimData(k)%reverse_axis = .false.
                                self%dimData(k)%negate = .true.
                            else if ((increase_flag .eqv. .false.) .and. (neg_flag .eqv. .true.)) then
                                self%dimData(k)%reverse_data = .true.
                                self%dimData(k)%reverse_axis = .true.
                                self%dimData(k)%negate = .false.
                            end if 
    
                            if (self%dimData(k)%reverse_axis .eqv. .true.) then
                                tempRealArray = tempRealArray(size(tempRealArray):1:-1)
                            end if
    
                            if (self%dimData(k)%negate .eqv. .true.) then
                                tempRealArray = -tempRealArray
                            end if
                        
                        end if
                        !need to check for 'time' variable specific issues
                        if (dimName == Globals%Var%time) then
                            call correctNCTime(dimUnits, tempRealArray)
                        end if
                        call dimsArrays(j)%initialize(dimName, dimUnits, 1, tempRealArray)
                        if (allocated(tempRealArray)) deallocate(tempRealArray)
                        if (allocated(tempRealArrayDelta)) deallocate(tempRealArrayDelta)
                    end if
                end do
            end do
        end if
    end do
    !write(*,*)"sai getvardimentions"
    end subroutine getVarDimensions

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
	!> @author Mohsen Shabani CRETUS - GFNL- 2025.11.12 | Email:shabani.mohsen@outlook.com
    !> @brief
    !> Reads the fields from the nc file for a given variable.
    !> returns a generic field, with a name, units and data
    !> @param[in] self, varName, varField, binaryVar, altName, altUnits
    !---------------------------------------------------------------------------
    subroutine getVar(self, varName, varField, binaryVar, altName, altUnits)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(generic_field_class), intent(out) :: varField
    logical, optional, intent(in) :: binaryVar
    type(string), optional, intent(in) :: altName, altUnits
    logical :: bVar
    real(prec), allocatable, dimension(:,:) :: tempRealField2D, auxRealField2D
    real(prec), allocatable, dimension(:,:,:) :: tempRealField3D, auxRealField3D
    real(prec), allocatable, dimension(:,:,:,:) :: tempRealField4D, auxRealField4D
    type(string) :: dimName, varUnits
    integer :: i, j, k, id_dim, t, indx, i1, j1, j2
    type(dim_t) :: tempDim, uDim
    integer, allocatable, dimension(:) :: varShape, u_Shape
    type(string) :: outext
    logical variable_u_is4D
    
    bVar= .false.
    if(present(binaryVar)) bVar = binaryVar
    variable_u_is4D = .false.

    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName ) then   !found the requested var
            !write(*,*)"Getting Var for simulation name = ", self%varData(i)%simName
            !write(*,*)"Getting Var name = ", varName
            allocate(varShape(self%varData(i)%ndims))
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                tempDim = self%getDimByDimID(self%varData(i)%dimids(j))
                varShape(j) = tempDim%length
            end do
            if(self%varData(i)%ndims == 3) then !3D variable
                allocate(auxRealField3D(varShape(1),varShape(2),varShape(3)))
                allocate(tempRealField3D(varShape(2),varShape(1),varShape(3))) !Write it as in MOHID hdf5s
                self%status = nf90_get_var(self%ncID, self%varData(i)%varid, auxRealField3D)
                !self%status = nf90_get_var(self%ncID, self%varData(i)%varid, tempRealField3D)
                do k=1,varShape(3)
                do i1=1,varShape(1)
                do j1=1,varShape(2)
                    tempRealField3D(j1,i1,k) = auxRealField3D(i1,j1,k)
                enddo
                enddo
                enddo
                
                call self%check()
                if (.not.bVar) then
                    where (tempRealField3D /= self%varData(i)%fillvalue)
                        tempRealField3D = tempRealField3D*self%varData(i)%scale + self%varData(i)%offset
                    elsewhere (tempRealField3D == self%varData(i)%fillvalue)
                        tempRealField3D = 0.0
                    end where
                else
                    if (self%varData(i)%fillvalue == MV) then
                        outext = '[NetCDFParser::getVar]:WARNING - variables without _fillvalue, you might have some problems in a few moments. Masks will not work properly (beaching, land exclusion,...)'
                    call Log%put(outext)
                    end if
                    where (tempRealField3D /= self%varData(i)%fillvalue) tempRealField3D = Globals%Mask%waterVal
                    where (tempRealField3D == self%varData(i)%fillvalue) tempRealField3D = Globals%Mask%landVal
                end if
                do id_dim=1,3 !reverting fields to have 'natural' coordinate-field storage
                    if (self%dimData(id_dim)%reverse_data) then
                        if (self%dimData(id_dim)%simName == Globals%Var%lon) then
                            tempRealField3D = tempRealField3D(varshape(1):1:-1,:,:)
                        else if (self%dimData(id_dim)%simName == Globals%Var%lat) then
                            tempRealField3D = tempRealField3D(:,varshape(2):1:-1,:)
                        else if (self%dimData(id_dim)%simName == Globals%Var%level) then
                            tempRealField3D = tempRealField3D(:,:,varshape(3):1:-1)
                        end if
                    end if
                end do
                
                !write(*,*)"Size temReakField3D = ", size(tempRealField3D,1), size(tempRealField3D,2), size(tempRealField3D,3)
                !write(*,*)"Value in temReakField4D(10,10,1) = ", tempRealField3D(10,10,1) 

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
                allocate(auxRealField4D(varShape(1),varShape(2),varShape(3),varShape(4)))
                allocate(tempRealField4D(varShape(2),varShape(1),varShape(3),varShape(4)))
                self%status = nf90_get_var(self%ncID, self%varData(i)%varid, auxRealField4D)
                !self%status = nf90_get_var(self%ncID, self%varData(i)%varid, tempRealField4D)
                do t=1,varShape(4)
                do k=1,varShape(3)
                do i1=1,varShape(1)
                do j1=1,varShape(2)
                    tempRealField4D(j1,i1,k,t) = auxRealField4D(i1,j1,k,t)
                enddo
                enddo
                enddo
                enddo
                !write(*,*)"Size temReakField4D = ", size(tempRealField4D,1), size(tempRealField4D,2), size(tempRealField4D,3), size(tempRealField4D,4)
                !write(*,*)"Value in temReakField4D(10,10,surface,1) = ", tempRealField4D(10,10,size(tempRealField4D,3),1) 
                call self%check()
                if (.not.bVar) then
                    if ((varName == Globals%Var%temp) .and. (Globals%simdefs%Temperature_add_offset /= 0)) then
                        where (tempRealField4D /= self%varData(i)%fillvalue)
                            tempRealField4D = tempRealField4D + Globals%simdefs%Temperature_add_offset
                        elsewhere (tempRealField4D == self%varData(i)%fillvalue)
                            tempRealField4D = 0.0
                        end where 
                    else
                        where (tempRealField4D /= self%varData(i)%fillvalue)
                            tempRealField4D = tempRealField4D*self%varData(i)%scale + self%varData(i)%offset
                        elsewhere (tempRealField4D == self%varData(i)%fillvalue)
                            tempRealField4D = 0.0
                        end where
                    end if
                else
                    if (self%varData(i)%fillvalue == MV) then
                        outext = '[NetCDFParser::getVar]:WARNING - variables without _fillvalue, you might have some problems in a few moments. Masks will not work properly (beaching, land exclusion,...)'
                    call Log%put(outext)
                    end if
                    where (tempRealField4D /= self%varData(i)%fillvalue) tempRealField4D = Globals%Mask%waterVal
                    where (tempRealField4D == self%varData(i)%fillvalue) tempRealField4D = Globals%Mask%landVal
                end if
                do id_dim=1,4 !reverting fields to have 'natural' coordinate-field storage
                    if (self%dimData(id_dim)%reverse_data) then
                        if (self%dimData(id_dim)%simName == Globals%Var%lon) then
                            tempRealField4D = tempRealField4D(varshape(1):1:-1,:,:,:)
                        elseif (self%dimData(id_dim)%simName == Globals%Var%lat) then
                            tempRealField4D = tempRealField4D(:,varshape(2):1:-1,:,:)
                        elseif (self%dimData(id_dim)%simName == Globals%Var%level) then
                            tempRealField4D = tempRealField4D(:,:,varshape(3):1:-1,:)
                        elseif (self%dimData(id_dim)%simName == Globals%Var%time) then
                            tempRealField4D = tempRealField4D(:,:,:,varshape(4):1:-1)
                        end if
                    end if
                end do
                
                if (.not.bVar) then
                    call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
                else
                    dimName = varName
                    if(present(altName)) dimName = altName
                    varUnits = self%varData(i)%units
                    if(present(altUnits)) varUnits = altUnits
                    call varField%initialize(dimName, varUnits, tempRealField4D)
                end if
            elseif(self%varData(i)%ndims == 2) then !2D variable, for now only bathymetry  !for the rugosity and D50 that has same structure like bathymetry
                if ((self%varData(i)%simName == Globals%Var%bathymetry) .or. (self%varData(i)%simName == Globals%Var%rugosityVar) .or. (self%varData(i)%simName == Globals%Var%D50Var) ) then
                    allocate(auxRealField2D(varShape(1),varShape(2)))
                    allocate(tempRealField2D(varShape(2),varShape(1)))
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
                                allocate(tempRealField4D(varShape(2),varShape(1), u_Shape(3), u_Shape(4)))
                                exit do1
                            else
                                allocate(tempRealField3D(varShape(2),varShape(1), u_Shape(3)))
                                exit do1
                            end if
                            
                        end if
                    end do do1
                    
                    self%status = nf90_get_var(self%ncID, self%varData(i)%varid, auxRealField2D)
                    do i1=1,varShape(1)
                    do j1=1,varShape(2)
                        tempRealField2D(j1,i1) = auxRealField2D(i1,j1)
                    enddo
                    enddo
                    call self%check()
                    if (.not.bVar) then
                        where (tempRealField2D /= self%varData(i)%fillvalue)
                            tempRealField2D = tempRealField2D*self%varData(i)%scale + self%varData(i)%offset
                        elsewhere (tempRealField2D == self%varData(i)%fillvalue)
                            tempRealField2D = 0.0
                        end where
                        if (variable_u_is4D) then
                            !For bathymetry, converts the 2D input field into a 4D field to be consistent with velocity matrixes
                            do t=1,size(tempRealField4D,4)
                                do k=1, size(tempRealField4D,3)
									if (self%varData(i)%simName == Globals%Var%bathymetry) tempRealField4D(:,:,k,t) = -tempRealField2D(:,:)
									if (self%varData(i)%simName == Globals%Var%rugosityVar) tempRealField4D(:,:,k,t) = +tempRealField2D(:,:)  
									if (self%varData(i)%simName == Globals%Var%D50Var) tempRealField4D(:,:,k,t) = +tempRealField2D(:,:)  
                                end do
                            end do
                            call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
                        else
                            !This needs to be completed .... will give errors
                            call varField%initialize(varName, self%varData(i)%units, tempRealField3D)
                        end if
                        
                    else
                        outext = '[NetCDFparser::getVar]: Variable '//varName//' is synthetic and cannot have 2D dimensionality. Stopping'
                        call Log%put(outext)
                        stop
                    end if
                else
                    outext = '[NetCDFparser::getVar]: Variable '//varName//' is 2D and not bathymetry-rugosity-D50, so it is not supported. Stopping'
                    call Log%put(outext)
                    stop
                end if
            else
                outext = '[NetCDFparser::getVar]: Variable '//varName//' has a non-supported dimensionality. Stopping'
                call Log%put(outext)
                stop
            end if
        end if
    end do

    end subroutine getVar

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> returns a dimension field metadata structure, given a dimID
    !> @param[in] self, dimID
    !---------------------------------------------------------------------------
    type(dim_t) function getDimByDimID(self, dimID)
    class(ncfile_class), intent(inout) :: self
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
        outext = '[NetCDFparser::getDimByDimID]: dimension with ID='//outext//' not found. Stopping'
        call Log%put(outext)
        stop
    end if

    end function getDimByDimID

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Close the netcdf file
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine closeFile(self)
    class(ncfile_class),intent(inout) :: self
    self%status = NF90_close(self%ncID)
    call self%check()
    end subroutine closeFile

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Debug the netcdf error after a netcdf command
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine check(self)
    class(ncfile_class), intent(inout) :: self
    type(string) :: outext
    if(self%status /= NF90_noerr) then
        outext = '[NetCDFparser::check]: '//trim(NF90_strerror(self%status))//', stoping'
        call Log%put(outext)
        stop
    end if
    end subroutine check

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> prints most of the file model metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printNcInfo(self)
    class(ncfile_class), intent(inout) :: self
    type(string) :: outext, temp
    outext = '-->Reading NetCDF file '//self%filename%chars()//new_line('a')
    temp = self%ncid
    outext = outext//'       file ID = '//temp//new_line('a')
    temp = self%nDims
    outext = outext//'       Number of dimensions = '//temp//new_line('a')
    temp = self%nVars
    outext = outext//'       Number of variable fields = '//temp
    call Log%put(outext,.false.)
    end subroutine printNcInfo

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> print variable metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printVarsNC(self)
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
    end subroutine printVarsNC

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> print dimensions metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine printDimsNC(self)
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
    end subroutine printDimsNC

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> corrects the time array to a more efficient format if needed.
    !> assumes netcdf comment as 'seconds since 1981-01-01 00:00:00'
    !> @param[in] timeComments, timeArray
    !---------------------------------------------------------------------------
    subroutine correctNCTime(timeComments, timeArray)
    type(string), intent(inout) :: timeComments
    real(prec), dimension(:), intent(inout) :: timeArray
    type(string), allocatable :: dc(:), dates(:), hours(:)
    type(string) :: isoDateStr
    integer, dimension(6) :: date
    type(datetime) :: NCDate
    type(timedelta) :: dateOffset
    type(string) :: outext
    real(prec) :: scale, offset

    call timeComments%split(tokens=dc, sep=' ')
    if (size(dc) == 4) then
        scale = 1.0
        if (dc(1) == 'seconds') scale = 1.0
        if (dc(1) == 'hours')   scale = 3600.0
        if (dc(1) == 'days')    scale = 3600.0*24.0
        if (dc(1) == 'months')  scale = 3600.0*24.0*30.0 !really hope no one gets such a brilliant idea as to use this as a time unit
        if (dc(1) == 'years')   scale = 3600.0*24.0*30.0*12.0 !or this
        call dc(3)%split(tokens=dates, sep='-')
        call dc(4)%split(tokens=hours, sep=':')
        isoDateStr = dates(1)//' '//dates(2)//' '//dates(3)//' '//hours(1)//' '//hours(2)//' '//hours(3)
        date = Utils%getDateFromISOString(isoDateStr)
        NCDate = Utils%getDateTimeFromDate(date)
        dateOffset = Globals%SimTime%StartDate - NCDate
        offset = -dateOffset%total_seconds()

        timeArray = timeArray*scale + offset
        timeComments = 'seconds since '//Globals%SimTime%StartDate%isoformat(' ')        
    
    elseif (size(dc) == 3) then
        scale = 1.0
        if (dc(1) == 'seconds') scale = 1.0
        if (dc(1) == 'hours')   scale = 3600.0
        if (dc(1) == 'days')    scale = 3600.0*24.0
        if (dc(1) == 'months')  scale = 3600.0*24.0*30.0 !really hope no one gets such a brilliant idea as to use this as a time unit
        if (dc(1) == 'years')   scale = 3600.0*24.0*30.0*12.0 !or this
        call dc(3)%split(tokens=dates, sep='-')
        isoDateStr = dates(1)//' '//dates(2)//' '//dates(3)//' '//'00'//' '//'00'//' '//'00'
        date = Utils%getDateFromISOString(isoDateStr)
        NCDate = Utils%getDateTimeFromDate(date)
        dateOffset = Globals%SimTime%StartDate - NCDate
        offset = -dateOffset%total_seconds()
        
        timeArray = timeArray*scale + offset
        timeComments = 'seconds since '//Globals%SimTime%StartDate%isoformat(' ')
        
    else       

        outext = '[NetCDF parser::correctNCTime]:WARNING - Time units may not be in the format *seconds since 1981-01-01 00:00:00*, you might have some problems in a few moments...'
        call Log%put(outext)
    end if

    end subroutine correctNCTime


    end module netcdfParser_mod