    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : simulation_testmaker
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : November 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Defines synthetic hydrodynamic fields to test Tracer behavior.
    !------------------------------------------------------------------------------

    module simulation_testmaker_mod

    use common_modules
    use background_mod
    use field_types_mod

    use nc_parser
    use ncparser_mod

    implicit none
    private

    type :: testmaker_class
        integer :: TestCode = 1
    contains
    procedure :: initialize => initTestMaker
    procedure, private :: makeTaylorGreen
    procedure, private :: makeConstantVel
    procedure, private :: makeRealVel
    end type testmaker_class

    type(testmaker_class) :: TestMaker

    !Public access vars
    public :: TestMaker

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Initializes the TestMaker
    !---------------------------------------------------------------------------
    subroutine initTestMaker(self, testcode, testbox, testbackground)
    class(testmaker_class), intent(inout) :: self
    integer, intent(in) :: testcode
    type(box), intent(in) :: testbox
    type(background_class), intent(inout) :: testbackground
    integer :: resolution
    type(string) :: outext
    resolution = 20
    self%TestCode = testcode
    if (self%TestCode == 1) then
        !make Taylor-Green vortices test
        outext = 'Creating Taylor-Green vortices, stand by...'
        call Log%put(outext)
        call self%makeTaylorGreen(resolution, testbox, testbackground)
    else if (self%TestCode == 2) then
        !make constant velocity test
        outext = 'Creating constant velocity fields, stand by...'
        call Log%put(outext)
        call self%makeConstantVel(resolution, testbox, testbackground)
    else if  (self%TestCode == 3) then
        !make real (from file) velocity test
        outext = 'Reading velocity fields from a netcdf file, stand by...'
        call Log%put(outext)
        call self%makeRealVel(resolution, testbox, testbackground)
    end if
    end subroutine initTestMaker

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Fills a domain's Background with taylor green vortices
    !---------------------------------------------------------------------------
    subroutine makeTaylorGreen(self, res, testbox, testbackground)
    class(testmaker_class), intent(inout) :: self
    integer, intent(in) :: res
    type(box), intent(in) :: testbox
    type(background_class), intent(inout) :: testbackground
    type(scalar1d_field_class), allocatable, dimension(:) :: testbackgroundims
    type(generic_field_class) :: gfield1, gfield2, gfield3
    type(string) :: name, units
    real(prec), allocatable, dimension(:) :: sx, sy, sz, t
    real(prec), allocatable, dimension(:,:,:,:) :: vx, vy, vz
    integer :: i, j, m
    real(prec) :: pi = 4*atan(1.0)
    real(prec) :: sst

    allocate(sx(res), sy(res), sz(res), t(res))
    allocate(vx(res,res,res,res))
    allocate(vy(res,res,res,res))
    allocate(vz(res,res,res,res))
    !create scaled coordinate arrays
    do i=1, res
        sx(i) = (i-1)*2.0*pi/(res-1)
        sz(i) = i-1
    end do
    sy = sx
    t = sz
    !create the velocity fields in time and space
    do m=1, res
        do i=1, res
            do j=1, res
                sst = 1.0 - m/res
                vx(i,j,:,m) = cos(sx(i))*sin(sy(j))*sst
                vy(i,j,:,m) = -sin(sx(i))*cos(sy(j))*sst
            end do
        end do
    end do
    vz = 0.0
    !put the data on generic fields
    units = 'm/s'
    call gfield1%initialize(Globals%Var%u, units, vx)
    call gfield2%initialize(Globals%Var%v, units, vy)
    call gfield3%initialize(Globals%Var%w, units, vz)
    !create the dimensions for the background
    allocate(testbackgroundims(4))
    units = 'deg'
    sx = sx*testbox%size%x/(2.0*pi) + testbox%pt%x
    call testbackgroundims(1)%initialize(Globals%Var%lon, units, 1, sx)
    units = 'deg'
    sy = sy*testbox%size%y/(2.0*pi) + testbox%pt%y
    call testbackgroundims(2)%initialize(Globals%Var%lat, units, 1, sy)
    units = 'm'
    sz = sz*testbox%size%z/(res-1) + testbox%pt%z
    call testbackgroundims(3)%initialize(Globals%Var%depth, units, 1, sz)
    units = 'seg'
    t = t*Globals%Parameters%TimeMax/(res-1)
    call testbackgroundims(4)%initialize(Globals%Var%time, units, 1, t)
    !construct background
    name = 'Taylor-Green test'
    testbackground = Background(1, name, testbox, testbackgroundims)
    call testbackground%add(gfield1)
    call testbackground%add(gfield2)
    call testbackground%add(gfield3)

    end subroutine makeTaylorGreen

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Fills a domain's Background with a constant velocity field
    !---------------------------------------------------------------------------
    subroutine makeConstantVel(self, res, testbox, testbackground)
    class(testmaker_class), intent(inout) :: self
    integer, intent(in) :: res
    type(box), intent(in) :: testbox
    type(background_class), intent(inout) :: testbackground
    type(scalar1d_field_class), allocatable, dimension(:) :: testbackgroundims
    type(generic_field_class) :: gfield1, gfield2, gfield3
    type(string) :: name, units
    real(prec), allocatable, dimension(:) :: sx, sy, sz, t
    real(prec), allocatable, dimension(:,:,:,:) :: vx, vy, vz
    integer :: i, j, k, m
    real(prec) :: pi = 4*atan(1.0)
    real(prec) :: sst

    allocate(sx(res), sy(res), sz(res), t(res))
    allocate(vx(res,res,res,res))
    allocate(vy(res,res,res,res))
    allocate(vz(res,res,res,res))
    !create scaled coordinate arrays
    do i=1, res
        sx(i) = i-1
    end do
    sy = sx
    sz = sx
    t = sz
    !create the velocity fields in time and space
    vx = 1.0
    vy = 2.0
    vz = 0.0
    !put the data on generic fields
    units = 'm/s'
    call gfield1%initialize(Globals%Var%u, units, vx)
    call gfield2%initialize(Globals%Var%v, units, vy)
    call gfield3%initialize(Globals%Var%w, units, vz)
    !create the dimensions for the background
    allocate(testbackgroundims(4))
    units = 'deg'
    sx = sx*testbox%size%x/(res-1) + testbox%pt%x
    call testbackgroundims(1)%initialize(Globals%Var%lon, units, 1, sx)
    units = 'deg'
    sy = sy*testbox%size%y/(res-1) + testbox%pt%y
    call testbackgroundims(2)%initialize(Globals%Var%lat, units, 1, sy)
    units = 'm'
    sz = sz*testbox%size%z/(res-1) + testbox%pt%z
    call testbackgroundims(3)%initialize(Globals%Var%depth, units, 1, sz)
    units = 'seg'
    t = t*Globals%Parameters%TimeMax/(res-1)
    call testbackgroundims(4)%initialize(Globals%Var%time, units, 1, t)
    !construct background
    name = 'Constant velocity test'
    testbackground = Background(1, name, testbox, testbackgroundims)
    call testbackground%add(gfield1)
    call testbackground%add(gfield2)
    call testbackground%add(gfield3)

    end subroutine makeConstantVel

    !---------------------------------------------------------------------------
    !> @author Daniel Garaboa - USC
    !> @brief
    !> Fills a domain's with sfc velocity field
    !---------------------------------------------------------------------------
    subroutine makeRealVel(self, res, testbox, testbackground)
    class(testmaker_class), intent(inout) :: self
    integer, intent(in) :: res
    type(box), intent(in) :: testbox
    type(background_class), intent(inout) :: testbackground
    type(scalar1d_field_class), allocatable, dimension(:) :: testbackgroundims
    type(generic_field_class) :: gfield1, gfield2,gfield3
    type(nc_class),dimension(3) :: nc_input_file
    type(string) :: name, units
    
    type(ncfile_class) :: ncFile
    type(string) :: ncFileName
    type(string) :: temp
    integer :: i
    type(scalar1d_field_class), allocatable, dimension(:) :: testbackgroundims_alt
    
    ncFileName = 'MOHID_Vigo_20180904_0000.nc4'
    call ncFile%initialize(ncFileName)
    temp = 'u'
    call ncFile%getVarDimensions(temp, testbackgroundims_alt)
    call ncFile%getVar(temp, gfield1)
    do i=1, size(testbackgroundims_alt)
        call testbackgroundims_alt(i)%print()
    end do
    call gfield1%print()

    nc_input_file(1)%file_name = 'MOHID_Vigo_20180904_0000.nc4'
    nc_input_file(2)%file_name = 'MOHID_Vigo_20180904_0000.nc4'
    nc_input_file(3)%file_name = 'MOHID_Vigo_20180904_0000.nc4'
    nc_input_file(1)%varname = 'u'
    nc_input_file(2)%varname = 'v'
    nc_input_file(3)%varname = 'w'

    call nc_input_file(1)%ncToField(gfield1)
    call nc_input_file(2)%ncToField(gfield2)
    call nc_input_file(3)%ncToField(gfield3)

    call nc_input_file(1)%printNcInfo()
    allocate(testbackgroundims(3))

    !create the dimensions for the background
    name = 'lon'
    units = 'deg'
    call testbackgroundims(1)%initialize(name, units, 1, nc_input_file(1)%dims(1)%values)
    name = 'lat'
    units = 'deg'
    call testbackgroundims(2)%initialize(name, units, 1, nc_input_file(1)%dims(2)%values)
    name = 'depth'
    units = 'm'
    call testbackgroundims(3)%initialize(name, units, 1, nc_input_file(1)%dims(3)%values)
    name = 'time'
    units = 's'
    call testbackgroundims(4)%initialize(name, units, 1, nc_input_file(1)%dims(4)%values)
    print*, nc_input_file(3)%dims(4)%values
    !construct background
    name = 'Real surface velocity field'
    testbackground = Background(1, name, testbox, testbackgroundims)
    call testbackground%add(gfield1)
    call testbackground%add(gfield2)
    call testbackground%add(gfield3)

    end subroutine makeRealVel


    end module simulation_testmaker_mod