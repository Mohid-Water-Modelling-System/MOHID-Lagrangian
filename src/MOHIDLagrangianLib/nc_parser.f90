    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : common_modules
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Daniel Garaboa Paz
    !
    ! DESCRIPTION:
    !> Module read netcdf in a easy way.
    !------------------------------------------------------------------------------

    module nc_parser

#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif

    use common_modules
    use background_mod
    use field_types_mod
    use xmlparser_mod
    use FoX_dom
    use background_mod

    implicit none

    type :: field
        real(prec),dimension(:),allocatable :: s1d
        real(prec),dimension(:,:),allocatable :: s2d
        real(prec),dimension(:,:,:),allocatable :: s3d
        real(prec),dimension(:,:,:,:),allocatable :: s4d
    end type field

    type :: dims
        real(prec),dimension(:),allocatable :: values
    end type dims

    type :: nc_class
        type(string) :: file_name, model_input_name
        type(string) :: varname, units
        integer :: ncid, varid
        integer :: ndims
        integer,dimension(:),allocatable :: var_dims
        type(field) :: variable
        type(dims),dimension(:),allocatable :: dims
        integer :: status
    contains
    procedure :: initNcLibHeaders
    procedure :: getNcid
    procedure :: getDimsNumber
    procedure :: getDimsShape
    procedure :: getDimsValues
    procedure :: getVarUnits
    procedure :: getVarName
    procedure :: getVarData
    procedure :: closeNcid
    procedure :: transferToGenericField
    procedure :: printNcInfo
    procedure :: check
    procedure :: ncToField
    end type nc_class

    contains

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Init the netcdf header form XML_INPUT
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine initNcLibHeaders(self, xmlfilename)
    implicit none
    class(nc_class),intent(inout) :: self
    type(string),dimension(2) :: model
    type(string),dimension(5) :: stdvarname
    type(string),intent(in) :: xmlfilename
    type(Node),pointer :: xmlmodel
    type(Node),pointer :: xmlvars
    type(string) :: tag, tag_name,  tag_value, varname
    integer :: flag, model_id, var_id, var_index, model_index

    model(1) = 'MOHID'
    model(2) = 'ROMS'

    stdvarname(1) = 'eastward_sea_water_velocity'
    stdvarname(2) = 'northward_sea_water_velocity'
    stdvarname(3) = 'upward_sea_water_velocity'
    stdvarname(4) = 'sea_water_temperature'
    stdvarname(5) = 'sea_water_salinity'

    call XMLReader%getFile(xmlmodel, xmlfilename)
    print*, 'I am here'
    tag = 'MOHID'
    call XMLReader%gotoNode(xmlvars, xmlvars, tag)
    print*,'but not here'

    do model_index = 1, size(model)
        print*,'Checking if the ouput is an input from', model(model_index)%chars()

        call XMLReader%gotonode(xmlvars, xmlmodel, model(model_index))
        tag_value = 'value'

        do var_index = 1, size(stdvarname)

            call XMLReader%getNodeAttribute(xmlmodel, stdvarname(var_index), tag_value, varname)
            self%status = NF90_inq_varid(self%ncid, varname%chars(), var_id)

            if (self%status == 0) then
                flag = flag + 1
            else
                stop
            end if

            if (flag == 5) then
                print*, 'Hey, this is a', model(model_index)%chars(), 'input. That means, the code is working'
            endif
        enddo

    enddo

    end subroutine initNcLibHeaders

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Open the netcdf file
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getNcid(self)
    class(nc_class),intent(inout) :: self
    self%status = NF90_open(trim(self%file_name%chars()), NF90_NOWRITE, self%ncid)
    call self%check()
    end subroutine getNcid

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> get the dimensions of the variable
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getDimsNumber(self)
    class(nc_class),intent(inout) :: self
    self%status = NF90_inq_varid(self%ncid, self%varname%chars(), self%varid)
    call self%check()
    self%status = nf90_inquire_variable(self%ncid, self%varid, ndims = self%ndims)
    call self%check()
    end subroutine

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> get the dimensions of the variable
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getDimsShape(self)
    class(nc_class),intent(inout) :: self
    integer,dimension(:),allocatable :: dimids
    integer :: include_parents
    integer :: i

    allocate(dimids(self%ndims))
    self%status = nf90_inquire_variable(self%ncid, self%varid, dimids = dimids)
    call self%check()

    allocate(self%var_dims(self%ndims))
    do i=1,self%ndims
        self%status = nf90_inquire_dimension(self%ncid, dimids(i), len=self%var_dims(i))
        call self%check()
    end do

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !>  Get the data extents
    !> @param[in] self
    !> WARNING: Incomplete command, working progress to prepare it to use with
    !> the background in order to sent the data to it!
    !---------------------------------------------------------------------------
    subroutine getDimsValues(self)
    class(nc_class),intent(inout) :: self
    integer,dimension(:),allocatable :: dimids
    integer :: i

    self%status = NF90_inq_varid(self%ncid, self%varname%chars(), self%varid)
    call self%check()

    allocate(dimids(self%ndims))

    self%status = nf90_inquire_variable(self%ncid, self%varid, dimids = dimids)
    call self%check()

    allocate(self%dims(self%ndims))

    do i=1,self%ndims
        allocate(self%dims(i)%values(self%var_dims(i)))
        self%status = NF90_GET_VAR(self%ncid, dimids(i), self%dims(i)%values)
    end do

    end subroutine

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Get the units of the variable from the attributte.
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getVarUnits(self)
    class(nc_class),intent(inout) :: self
    character(len=:),allocatable :: units
    integer :: units_len, attnum
    self%status = NF90_inquire_attribute(self%ncid, self%varid, 'units', len=units_len, attnum=attnum)
    call self%check()
    allocate(character(len=units_len) :: units)
    self%status = nf90_get_att(self%ncid, self%varid, 'units', units)
    self%units=units
    call self%check()
    end subroutine getVarUnits

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Get the standard name of the variable
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getVarName(self)
    class(nc_class),intent(inout) :: self
    character(len=:),allocatable :: varname
    integer :: name_len, attnum
    self%status = NF90_inquire_attribute(self%ncid, self%varid, 'standard_name', len=name_len, attnum=attnum)
    call self%check()
    allocate(character(len=name_len) ::varname)
    self%status = nf90_get_att(self%ncid, self%varid, 'standard_name', varname)
    ! for a strange reason the netcdf output does not admit the string as argument.
    self%varname=varname
    call self%check()
    end subroutine getVarName

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> get the field data and added it to the field class.
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine getVarData(self)
    class(nc_class),intent(inout) :: self

    if (self%ndims == 1) then
        allocate(self%variable%s1d(self%var_dims(1)))
        self%status = NF90_GET_VAR(self%ncid, self%varid, self%variable%s1d)
        call self%check()

    else if (self%ndims == 2) then
        allocate(self%variable%s2d(self%var_dims(1),self%var_dims(2)))
        self%status = NF90_GET_VAR(self%ncid, self%varid, self%variable%s2d)
        call self%check()

    else if (self%ndims == 3) then
        allocate(self%variable%s3d(self%var_dims(1),self%var_dims(2),self%var_dims(3)))
        self%status = NF90_GET_VAR(self%ncid, self%varid, self%variable%s3d)
        call self%check()

    else if(self%ndims == 4) then
        allocate(self%variable%s4d(self%var_dims(1),self%var_dims(2),self%var_dims(3),self%var_dims(4)))

        self%status = NF90_GET_VAR(self%ncid, self%varid, self%variable%s4d)
        call self%check()

    end if
    end subroutine getVarData

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Close the netcdf file
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine closeNcid(self)
    class(nc_class),intent(inout) :: self
    self%status = NF90_close(self%ncid)
    call self%check()
    end subroutine closeNcid

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Debug the netcdf error after a netcdf command
    !> @param[in] self
    !
    !---------------------------------------------------------------------------
    subroutine check(self)
    class(nc_class), intent(inout) :: self
    if(self%status /= NF90_noerr) then
        print *, trim(NF90_strerror(self%status))
        stop "Stopped"
    end if
    end subroutine check

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Debug the netcdf error after a netcdf command
    !> @param[in] self
    !> The generic initialize methods does not work
    !---------------------------------------------------------------------------
    subroutine transferToGenericField(self,gfield)
    class(nc_class), intent(inout) :: self
    type(generic_field_class),intent(out) :: gfield    
    if (self%ndims == 1) then
        call gfield%initialize(self%varname, self%units, self%variable%s1d)    
    else if (self%ndims == 2) then
        call gfield%initialize(self%varname, self%units, self%variable%s2d)    
    else if (self%ndims == 3) then
        call gfield%initialize(self%varname, self%units, self%variable%s3d)    
    else if(self%ndims == 4) then
        call gfield%initialize(self%varname, self%units, self%variable%s4d)
    end if
    end subroutine transferToGenericField

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Read a field and transfer it to a gfield
    !> @param[in] self
    !> The generic initialize methods does not work
    !---------------------------------------------------------------------------
    subroutine ncToField(self,gfield)
    class(nc_class), intent(inout) :: self
    type(generic_field_class),intent(out) :: gfield
    type(string) :: xmlfile
    xmlfile = './nc_institution_library.xml'
    call self%getNcid()
    call self%initNcLibHeaders(xmlfile)
    call self%getDimsNumber()
    call self%getDimsShape()
    call self%getDimsValues()
    call self%getVarUnits()
    call self%getVarName()
    call self%getVarData()
    call self%closeNcid()
    call self%TransferToGenericField(gfield)
    end subroutine ncToField

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa Paz - USC
    !> @brief
    !> Read a field and transfer it to a gfield
    !> @param[in] self
    !> print the main nc information to chekc everytthing is fine
    !---------------------------------------------------------------------------
    subroutine printNcInfo(self)
    class(nc_class), intent(inout) :: self

    print*, 'Id Nc opened', self%file_name%chars()
    print*, 'Id Nc opened', self%ncid
    print*, 'Number of dimensions', self%ndims
    print*, 'Dimensions size', self%var_dims
    print*, 'Field readed name', self%varname%chars()
    end subroutine printNcInfo

    end module nc_parser