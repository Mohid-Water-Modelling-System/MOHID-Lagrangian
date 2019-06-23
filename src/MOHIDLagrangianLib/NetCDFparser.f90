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
    use background_mod

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
        integer :: status
    contains
    procedure :: initialize => getFile
    procedure :: getVarDimensions
    procedure :: getVar
    procedure :: finalize => closeFile
    procedure, private :: check
    procedure, private :: getNCid
    procedure, private :: getNCglobalMetadata
    procedure, private :: getNCDimMetadata
    procedure, private :: getNCVarMetadata
    procedure, private :: getDimByDimID
    procedure :: print => printNcInfo
    end type ncfile_class

    public :: ncfile_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Parse the netcdf file headers and assembles the file model
    !> @param[in] self, filename
    !---------------------------------------------------------------------------
    subroutine getFile(self, filename)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: filename
    integer :: i
    self%filename = filename
    call self%getNCid()
    call self%getNCglobalMetadata()
    call self%getNCVarMetadata()
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
    integer :: i, j, ndims, nAtts
    integer :: dimids(self%nDims)
    integer :: tempStatus
    character(CHAR_LEN) :: varName, units
    allocate(self%varData(self%nVars))
    do i=1, self%nVars
        self%status = nf90_inquire_variable(self%ncID, i, varName, ndims=ndims, dimids=dimids, nAtts=nAtts)
        call self%check()
        self%varData(i)%name = trim(varName)
        self%varData(i)%simName = Globals%Var%getVarSimName(self%varData(i)%name)
        self%varData(i)%varid = i
        self%varData(i)%ndims = ndims
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
        if (tempStatus == -43) self%varData(i)%fillvalue = 0.0
    end do
    end subroutine getNCVarMetadata

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Inquires the nc file for dimension metadata
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNCDimMetadata(self)
    class(ncfile_class), intent(inout) :: self
    integer :: i, j, dimLength
    character(CHAR_LEN) :: dimName
    allocate(self%dimData(self%nDims))
    do i=1, self%nDims
        self%status = nf90_inquire_dimension(self%ncID, i, dimName, dimLength)
        call self%check()
        self%dimData(i)%name = trim(dimName)
        self%dimData(i)%simName = Globals%Var%getVarSimName(self%dimData(i)%name)
        self%dimData(i)%length = dimLength
        self%status = nf90_inq_dimid(self%ncID, self%dimData(i)%name%chars(), self%dimData(i)%dimid)
        call self%check()
        do j=1, self%nVars
            if (self%dimData(i)%name == self%varData(j)%name) then
                self%dimData(i)%units = self%varData(j)%units
                self%dimData(i)%varid = self%varData(j)%varid
                exit
            end if
        end do
    end do
    end subroutine getNCDimMetadata

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Reads the dimension fields from the nc file for a given variable.
    !> returns an array of scalar 1D fields, each with a name, units and data
    !> @param[in] self, varName, dimsArrays
    !---------------------------------------------------------------------------
    subroutine getVarDimensions(self, varName, dimsArrays)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(scalar1d_field_class), allocatable, dimension(:), intent(out) :: dimsArrays
    real(prec), allocatable, dimension(:) :: tempRealArray, tempRealArrayDelta
    type(string) :: dimName, dimUnits
    integer :: i, j, k, l

    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName) then   !found the requested var
            allocate(dimsArrays(self%varData(i)%ndims)) !allocating output fields
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                do k=1, self%nDims  !going trough all available dimensions of the file
                    if (self%varData(i)%dimids(j) == self%dimData(k)%dimid) then    !found a corresponding dimension between the variable and the file
                        allocate(tempRealArray(self%dimData(k)%length)) !allocating a place to read the field data to
                        dimName = self%dimData(k)%simName
                        dimUnits = self%dimData(k)%units
                        self%status = nf90_get_var(self%ncID, self%dimData(k)%varid, tempRealArray)
                        call self%check()
                        !need to check for 'level' variable specific issues

                        !@Daniel
                        ! We seek for: 0) index=(0, ...n), axis=(+surface=0 .... -bottom)
                        ! Three situations:
                        !     index=(0,1,2,    ...  n-1,n)  
                        ! 1)  axis=(-bottom, ... ,-surface)   --> Transform: reverse
                        ! 2)  axis=(+bottom, ... , +bottom)  --> Transform: reverse and negate
                        ! 3)  axis=(+surface, ... , +bottom)  --> Transform: negate
                     
                        if (dimName == Globals%Var%level) then
                            if ((tempRealArray(1) <= 0) .and. (tempRealArray(1)) > tempRealArray(size(tempRealArray))) then
                                self%dimData(k)%reverse_axis = .false.
                                self%dimData(k)%negate = .false.
                            elseif ((tempRealArray(1) <= 0) .and. (tempRealArray(1)) < tempRealArray(size(tempRealArray))) then
                                self%dimData(k)%reverse_axis = .true.
                                self%dimData(k)%negate = .false.
                                self%dimData(k)%reverse_axis = .true.
                                !print*, '[NetCDFparser::warning]:', 'The axis',k,'has wrong directon. Correcting...'
                            elseif ((tempRealArray(1) >= 0) .and. (tempRealArray(1) > tempRealArray(size(tempRealArray)))) then
                                self%dimData(k)%reverse_axis = .false.
                                self%dimData(k)%negate = .true.
                                self%dimData(k)%reverse_data = .false.
                                !print*, '[NetCDFparser::warning]:', 'The axis',k,'has wrong sing/direction. Correcting...'  
                            elseif ((tempRealArray(1) >= 0) .and. (tempRealArray(1) < tempRealArray(size(tempRealArray)))) then
                                self%dimData(k)%reverse_axis = .true.
                                self%dimData(k)%negate = .true.
                                self%dimData(k)%reverse_data = .true.
                                !print*, '[NetCDFparser::warning]:', 'The axis',k,'has wrong sign. Correcting...'
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

    end subroutine getVarDimensions

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Reads the fields from the nc file for a given variable.
    !> returns a generic field, with a name, units and data
    !> @param[in] self, varName, varField
    !---------------------------------------------------------------------------
    subroutine getVar(self, varName, varField)
    class(ncfile_class), intent(inout) :: self
    type(string), intent(in) :: varName
    type(generic_field_class), intent(out) :: varField
    real(prec), allocatable, dimension(:) :: tempRealField1D
    real(prec), allocatable, dimension(:,:,:) :: tempRealField3D
    real(prec), allocatable, dimension(:,:,:,:) :: tempRealField4D
    type(string) :: dimName, varUnits
    integer :: i, j, k, id_dim, first,last
    type(dim_t) :: tempDim
    integer, allocatable, dimension(:) :: varShape
    type(string) :: outext

    do i=1, self%nVars !going trough all variables
        if (self%varData(i)%simName == varName ) then   !found the requested var
            allocate(varShape(self%varData(i)%ndims))
            do j=1, self%varData(i)%ndims   !going trough all of the variable dimensions
                tempDim = self%getDimByDimID(self%varData(i)%dimids(j))
                varShape(j) = tempDim%length
            end do
            if(self%varData(i)%ndims == 3) then !3D variable
                allocate(tempRealField3D(varShape(1),varShape(2),varShape(3)))
                self%status = nf90_get_var(self%ncID, self%varData(i)%varid, tempRealField3D)
                call self%check()
                where (tempRealField3D == self%varData(i)%fillvalue) tempRealField3D = 0.0
                tempRealField3D = tempRealField3D*self%varData(i)%scale + self%varData(i)%offset ! scale + offset transform
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
                call varField%initialize(varName, self%varData(i)%units, tempRealField3D)
            else if(self%varData(i)%ndims == 4) then !4D variable
                allocate(tempRealField4D(varShape(1),varShape(2),varShape(3),varShape(4)))
                self%status = nf90_get_var(self%ncID, self%varData(i)%varid, tempRealField4D)
                call self%check()                
                where (tempRealField4D == self%varData(i)%fillvalue) tempRealField4D = 0.0
                tempRealField4D = tempRealField4D*self%varData(i)%scale + self%varData(i)%offset
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
                call varField%initialize(varName, self%varData(i)%units, tempRealField4D)
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
    integer :: i
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
        NCDate = datetime(date(1),date(2),date(3),date(4),date(5),date(6))
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