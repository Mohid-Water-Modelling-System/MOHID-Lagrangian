    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : source_identity
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : March 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a source class and related methods.
    !------------------------------------------------------------------------------

    module sources_mod

    use common_modules
    use simulationGlobals_mod
    use csv_module
    use csvParser_mod
    use MTimeSeriesParser_mod

    implicit none
    private

    type :: source_par                          !<Type - parameters of a source object
        type(string) :: name                    !< source name
        integer :: id                           !< unique source identification (integer)
        real(prec) :: emitting_rate             !< Emitting rate of the Source (Hz)
        logical :: emitting_fixed_rate          !< Type of emitter rate: true-fixed rate(Hz); false-variable(from file)
        type(string) :: rate_file               !< File name of the emission rate data (csv, .dat)
        real(prec), dimension(:), allocatable :: variable_rate !< Emission rate read from the rate_file of the Source - interpolated on a regular time series spaced dt
        real(prec) :: rateScale                 !< Scaling for the variable emission rate
        logical :: fixed_position               !< Type of position mode: true-fixed position; false-variable position(from file)
        type(string) :: position_file           !< File name of the position data (csv, .dat)
        type(vector), dimension(:), allocatable :: variable_position !< Position read from the position_file of the Source - interpolated on a regular time series spaced dt
        logical, dimension(:), allocatable :: activeTime      !< array of logicals that maps active state for every dt
        type(string) :: source_geometry         !< Source type : 'point', 'line', 'sphere', 'box'
        class(shape), allocatable :: geometry   !< Source geometry
    end type source_par

    type :: source_prop                         !<Type - material properties of a source object
        type(string) :: propertyType            !< source property type (plastic, paper, fish, etc)
        type(string) :: propertySubType         !< source property name
        logical :: particulate                  !< true for a Source that emitts particulate tracers (a concentration of particles)
        real(prec) :: radius                    !< radius of the emitted Tracers (size of the particle if not particulate, volume of the Tracer if particulate)
        real(prec) :: volume                    !< volume of the emitted particles
        real(prec) :: area                      !< surface area  of the emitted particles
        real(prec) :: density                   !< density of the Tracers
        type(string), dimension(:), allocatable :: propName !< name of a given property
        real(prec), dimension(:), allocatable :: propValue  !< value of a given property
    end type source_prop

    type :: source_state             !<Type - state variables of a source object
        real(prec) :: age                   ! time variables
        logical :: active                   !< active switch
        real(prec) :: emission_stack        !< number of emissions on the stack for the current time step
        type(vector) :: pos                 !< Position of the source baricenter (m)
        type(vector) :: vel                 !< Velocity of the source (m s-1)
    end type source_state

    type :: source_stats             !<Type - statistical variables of a source object
        integer :: particles_emitted        !< Number of emitted particles by this source
        integer :: ns                       !< Number of sampling steps
    end type source_stats

    type :: source_stencil         !<Type - holder for the tracer creation stencil of the source
        type(vector) :: dp                  !< resolution in x, y and z directions
        integer :: np                       !< Number of tracers by emission
        integer :: total_np                 !< Total number of tracers that this source will generate
        type(vector), allocatable, dimension(:) :: ptlist !<list of points (coordinates), relative to the source geometry point, to be generated at every emission.
    end type source_stencil

    type :: source_class           !<Type - The source class
        type(source_par)   :: par           !<To access parameters
        type(source_prop)  :: prop          !<To access Tracer properties
        type(source_state) :: now           !<To access state variables
        type(source_stencil) :: stencil     !<To acess stencil variables
        type(source_stats) :: stats         !<To access statistics
    contains
    procedure :: initialize => initializeSource
    procedure :: isParticulate
    procedure :: setPropertyNumber
    procedure :: setPropertyBaseAtribute
    procedure :: setPropertyAtribute
    procedure :: check
    procedure :: setVariableRate
    procedure :: setVariablePosition
    procedure, private :: getVariableRate
    procedure, private :: getVariablePosition
    procedure, private :: setActiveTimes
    procedure, private :: setotalnp
    procedure, private :: linkProperty
    procedure :: print => printSource
    end type source_class

    type :: sourceArray_class
        type(source_class), allocatable, dimension(:) :: src
    contains
    procedure :: initialize => initSources
    procedure :: finalize => killSources
    procedure :: setPropertyNames
    end type sourceArray_class

    !Simulation variables
    type(sourceArray_class) :: tempSources !< Temporary Source array, used exclusively for building the case from a description file

    !Public access vars
    public :: tempSources, sourceArray_class, source_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source allocation routine - allocates sources objects
    !> @param[in] self,nsources
    !---------------------------------------------------------------------------
    subroutine initSources(self,nsources)
    implicit none
    class(sourceArray_class), intent(inout) :: self
    integer, intent(in) :: nsources
    integer err
    type(string) :: outext, temp
    allocate(self%src(nsources), stat=err)
    if(err/=0)then
        outext='[initSources]: Cannot allocate Sources, stoping'
        call Log%put(outext)
        stop
    else
        temp = nsources
        outext = 'Allocated '// temp // ' Sources.'
        call Log%put(outext)
    endif
    end subroutine initSources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source group destructor - deallocates sources objects
    !---------------------------------------------------------------------------
    subroutine killSources(self)
    implicit none
    class(sourceArray_class), intent(inout) :: self
    integer err
    type(string) :: outext
    if (allocated(self%src)) deallocate(self%src, stat=err)
    if(err/=0)then
        outext='[killSources]: Cannot deallocate Sources, stoping'
        call Log%put(outext)
        stop
    endif
    end subroutine killSources

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source property setting procedure - initializes Source variables
    !> @param[in] src, ptype, pname
    !---------------------------------------------------------------------------
    subroutine linkProperty(src, ptype, pname)
    implicit none
    class(source_class), intent(inout) :: src
    type(string), intent(in) :: ptype
    type(string), intent(in) :: pname
    src%prop%propertyType = ptype
    src%prop%propertySubType = pname
    end subroutine linkProperty

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source property setting routine, calls source by id to set its properties
    !> @param[in] self,srcid_str,ptype,pname
    !---------------------------------------------------------------------------
    subroutine setPropertyNames(self, srcid_str, ptype, pname)
    implicit none
    class(sourceArray_class), intent(inout) :: self
    type(string), intent(in) :: srcid_str      !<Source id tag
    type(string), intent(in) :: ptype          !<Property type to set
    type(string), intent(in) :: pname          !<Property name to set
    integer :: srcid
    type(string) :: outext, temp
    integer :: i
    logical :: notlinked
    srcid = srcid_str%to_number(kind=1_I1P)
    notlinked = .true.  !assuming not linked
    do i=1, size(self%src)
        if (self%src(i)%par%id == srcid) then ! found the correct source to link to
            call self%src(i)%linkProperty(ptype,pname) ! calling Source method to link property
            temp = self%src(i)%par%id
            outext='      Source id = '// temp // ', '// self%src(i)%par%name //' is of type '// self%src(i)%prop%propertyType //', with property name ' // self%src(i)%prop%propertySubType
            call Log%put(outext,.false.)
            notlinked = .false. ! we linked it
            exit
        endif
    enddo
    if (notlinked) then ! property has no corresponding Source
        temp = srcid
        outext='      Source id = '// temp //' not listed, property '// pname //', of type ' // ptype // ' not linked, ignoring'
        call Log%put(outext,.false.)
    endif
    end subroutine setPropertyNames

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Allocates the space for the properties
    !> @param[in] src, nProps
    !---------------------------------------------------------------------------
    subroutine setPropertyNumber(src, nProps)
    class(source_class), intent(inout) :: src
    integer, intent(in) :: nProps
    if(allocated(src%prop%propName)) deallocate(src%prop%propName)
    if(allocated(src%prop%propValue)) deallocate(src%prop%propValue)
    allocate(src%prop%propValue(nProps))
    allocate(src%prop%propName(nProps))
    end subroutine setPropertyNumber

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source property atribute setting proceadure - initializes Source variables
    !> @param[in] src, i, pName, pValue
    !---------------------------------------------------------------------------
    subroutine setPropertyAtribute(src, i, pName, pValue)
    class(source_class), intent(inout) :: src
    type(string), intent(in) :: pName
    type(string), intent(in) :: pValue
    integer, intent(in) :: i
    type(string) :: outext
    src%prop%propName(i) = pName
    src%prop%propValue(i) = pValue%to_number(kind=1._R8P)
    end subroutine setPropertyAtribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source property atribute setting proceadure - initializes Source variables
    !> @param[in] src, pname, pvalue
    !---------------------------------------------------------------------------
    subroutine setPropertyBaseAtribute(src, pname, pvalue)
    class(source_class), intent(inout) :: src
    type(string), intent(in) :: pname
    type(string), intent(in) :: pvalue
    type(string) :: outext
    select case (pname%chars())
    case ('particulate')
        if (pvalue%to_number(kind=1_I1P) == 1) then
            src%prop%particulate = .true.
        end if
    case ('radius')
        src%prop%radius = pvalue%to_number(kind=1._R8P)
    case ('volume')
        src%prop%volume = pvalue%to_number(kind=1._R8P)
    case ('area')
        src%prop%area = pvalue%to_number(kind=1._R8P)
    case ('density')
        src%prop%density = pvalue%to_number(kind=1._R8P)
        case default
        outext='[Sources::setPropertyBaseAtribute]: unexpected atribute '//pname//' for property '//src%prop%propertySubType//', ignoring'
        call Log%put(outext)
    end select
    end subroutine setPropertyBaseAtribute

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that checks for the consistency of the Source properties.
    !---------------------------------------------------------------------------
    subroutine check(self)
    class(source_class), intent(in) :: self
    type(string) :: outext, temp(2)
    logical :: failed
    failed = .false.
    temp(1) = self%par%id
    if (self%prop%radius == MV) then
        failed = .true.
        temp(2) = 'radius'
    elseif (self%prop%volume == MV) then
        failed = .true.
        temp(2) = 'volume'
    elseif (self%prop%area == MV) then
        failed = .true.
        temp(2) = 'area'
    elseif (self%prop%density == MV) then
        failed = .true.
        temp(2) = 'density'
    end if
    if (failed) then
        outext = 'Property '//temp(2)//' from Source id = '//temp(1)//' is not set in the material property library file, stoping'
        call Log%put(outext)
        stop
    end if
    end subroutine check

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that reads a file with variable emission rate data, allocates
    !> the space, and stores the data in the Source.
    !> @param[in] self, filename, rateScale
    !---------------------------------------------------------------------------
    subroutine getVariableRate(self, filename, rateScale)
    class(source_class), intent(inout) :: self
    type(string), intent(in) :: filename
    real(prec), intent(in) :: rateScale
    type(csv_file) :: rateFile
    character(len=30), dimension(:), allocatable :: header
    real(prec), dimension(:), allocatable :: time, rate, stime, srate
    integer :: i, index
    real(prec) :: weight
    type(csvparser_class) :: CSVReader
    type(mTimeSeriesParser_class) :: MTSReader
    type(string), dimension(:), allocatable :: varList
    type(string) :: outext

    if (filename%extension() == '.csv') then
        call CSVReader%getFile(filename, 1)
        call CSVReader%getDataByLabel(Globals%Var%time, time)
        call CSVReader%getDataByLabel(Globals%Var%rate, rate)
    else if (filename%extension() == '.dat') then
        allocate(varList(2))
        varList(1) = Globals%Var%time
        varList(2) = Globals%Var%rate
        call MTSReader%getFile(filename, varList)
        call MTSReader%getDataByLabel(Globals%Var%time, time)
        call MTSReader%getDataByLabel(Globals%Var%rate, rate)
    end if
    
    rate = rate*rateScale
    
    !interpolating the data to a regular dt spaced array and storing the data in the Source
    allocate(stime(int(min(Globals%Parameters%TimeMax,maxval(time))/Globals%SimDefs%dt)+1))
    allocate(srate(size(stime)))
    do i=1, size(stime)
        stime(i)=Globals%SimDefs%dt*(i-1)
    end do
    do i=1, size(stime)
        if (stime(i)<time(1)) then
            srate(i) = 0.0
        else
            do index=1, size(time)-1
                if (stime(i)>=time(index)) then
                    weight = (stime(i)-time(index))/(time(index+1)-time(index))
                    srate(i) = (1.0-weight)*rate(index) + weight*rate(index+1)
                end if
            end do
        end if
    end do
    allocate(self%par%variable_rate(size(srate)))
    self%par%variable_rate = srate
    end subroutine getVariableRate

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the rate for the current time step, in case of a variable
    !> rate.
    !> @param[in] self, index
    !---------------------------------------------------------------------------
    subroutine setVariableRate(self, index)
    class(source_class), intent(inout) :: self
    integer :: index
    self%par%emitting_rate = self%par%variable_rate(max(1,min(index,size(self%par%variable_rate))))
    end subroutine setVariableRate
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that reads a file with variable position data, allocates
    !> the space, and stores the data in the Source.
    !> @param[in] self, filename
    !---------------------------------------------------------------------------
    subroutine getVariablePosition(self, filename)
    class(source_class), intent(inout) :: self
    type(string), intent(in) :: filename
    type(csv_file) :: csvFile
    character(len=30), dimension(:), allocatable :: header
    real(prec), dimension(:), allocatable :: time, stime
    real(prec), dimension(:), allocatable :: xx, yy, zz
    real(prec), dimension(:,:), allocatable :: posi, sposi
    integer :: i, index
    real(prec) :: weight
    type(csvparser_class) :: CSVReader
    type(mTimeSeriesParser_class) :: MTSReader
    type(string), dimension(:), allocatable :: varList
    type(string) :: outext

    if (filename%extension() == '.csv') then
        call CSVReader%getFile(filename, 1)
        call CSVReader%getDataByLabel(Globals%Var%time, time)
        call CSVReader%getDataByLabel(Globals%Var%lon, xx)
        call CSVReader%getDataByLabel(Globals%Var%lat, yy)
        call CSVReader%getDataByLabel(Globals%Var%level, zz)
    else if (filename%extension() == '.dat') then
        allocate(varList(4))
        varList(1) = Globals%Var%time
        varList(2) = Globals%Var%lon
        varList(3) = Globals%Var%lat
        varList(4) = Globals%Var%level
        call MTSReader%getFile(filename, varList)
        call MTSReader%getDataByLabel(Globals%Var%time, time)
        call MTSReader%getDataByLabel(Globals%Var%lon, xx)
        call MTSReader%getDataByLabel(Globals%Var%lat, yy)
        call MTSReader%getDataByLabel(Globals%Var%level, zz)
    end if
    
    allocate(posi(size(time),3))
    posi(:,1) = xx
    posi(:,2) = yy
    posi(:,3) = zz
    
    !interpolating the data to a regular dt spaced array and storing the data in the Source
    allocate(stime(int(min(Globals%Parameters%TimeMax,maxval(time))/Globals%SimDefs%dt)+1))
    do i=1, size(stime)
        stime(i)=Globals%SimDefs%dt*(i-1)
    end do
    allocate(sposi(size(stime),3))
    do i=1, size(stime)
        if (stime(i)<time(1)) then
            sposi(i,1) = self%par%geometry%pt%x
            sposi(i,2) = self%par%geometry%pt%y
            sposi(i,3) = self%par%geometry%pt%z
        else
            do index=1, size(time)-1
                if (stime(i)>=time(index)) then
                    weight = (stime(i)-time(index))/(time(index+1)-time(index))
                    sposi(i,1) = (1.0-weight)*posi(index,1) + weight*posi(index+1,1)
                    sposi(i,2) = (1.0-weight)*posi(index,2) + weight*posi(index+1,2)
                    sposi(i,3) = (1.0-weight)*posi(index,3) + weight*posi(index+1,3)
                end if
            end do
        end if
    end do
    allocate(self%par%variable_position(size(stime)))
    do i=1, size(stime)
        self%par%variable_position(i) = sposi(i,1)*ex + sposi(i,2)*ey + sposi(i,3)*ez
    end do
    end subroutine getVariablePosition
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets the position for the current time step, in case of a variable
    !> position.
    !> @param[in] self, index
    !---------------------------------------------------------------------------
    subroutine setVariablePosition(self, index)
    class(source_class), intent(inout) :: self
    integer :: index
    self%now%pos = self%par%variable_position(max(1,min(index,size(self%par%variable_position))))
    end subroutine setVariablePosition

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that sets active times logical array of the Source.
    !> @param[in] self, activeTimes
    !---------------------------------------------------------------------------
    subroutine setActiveTimes(self, activeTimes)
    class(source_class), intent(inout) :: self
    real(prec), dimension(:,:), intent(in) :: activeTimes
    integer :: i, j
    real(prec) :: time
    allocate(self%par%activeTime(int(Globals%Parameters%TimeMax/Globals%SimDefs%dt)+1))
    if (size(activeTimes, 1) < 1) then
        self%par%activeTime = .true.
        return
    end if
    self%par%activeTime = .false.
    do i=1, size(self%par%activeTime)
        time = i*Globals%SimDefs%dt
        do j=1, size(activeTimes, 1)
            if (time >= activeTimes(j,1)) then
                if (time <= activeTimes(j,2)) then
                    self%par%activeTime(i) = .true.
                end if
            end if
        end do
    end do
    end subroutine setActiveTimes

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source inititialization proceadure - initializes Source variables
    !> @param[in] src, id, name, emitting_rate, emitting_fixed_rate, rate_file, rateScale, posi_fixed, posi_file, activeTimes, source_geometry, shapetype, res
    !---------------------------------------------------------------------------
    subroutine initializeSource(src, id, name, emitting_rate, emitting_fixed_rate, rate_file, rateScale, posi_fixed, posi_file, activeTimes, source_geometry, shapetype, res)
    class(source_class) :: src
    integer, intent(in) :: id
    type(string), intent(in) :: name
    real(prec), intent(in) :: emitting_rate
    logical, intent(in) :: emitting_fixed_rate
    type(string), intent(in) :: rate_file
    real(prec), intent(in) :: rateScale
    logical, intent(in) :: posi_fixed
    type(string), intent(in) :: posi_file
    real(prec), dimension(:,:), intent(in) :: activeTimes
    type(string), intent(in) :: source_geometry
    class(shape), intent(in) :: shapetype
    type(vector), intent(in) :: res
    integer :: sizem, i
    type(string) :: outext
    integer :: err
    !Setting parameters
    src%par%id=id
    src%par%name=name
    src%par%emitting_rate=emitting_rate
    src%par%emitting_fixed_rate = emitting_fixed_rate
    src%par%rate_file = rate_file
    if (.not.emitting_fixed_rate) then
        call src%getVariableRate(src%par%rate_file, rateScale)
        src%par%emitting_rate=sum(src%par%variable_rate)/size(src%par%variable_rate)
    end if
    src%par%rateScale = rateScale
    !Setting possible variable position
    src%par%fixed_position = posi_fixed
    src%par%position_file = posi_file
    if (.not.posi_fixed) then
        call src%getVariablePosition(posi_file)
    end if
    call src%setActiveTimes(activeTimes)
    src%par%source_geometry=source_geometry
    allocate(src%par%geometry, source=shapetype)
    !Setting properties
    src%prop%propertyType = "base" ! pure Lagrangian trackers by default
    src%prop%propertySubType = "base"
    src%prop%particulate = .false.
    src%prop%radius = MV
    src%prop%volume = MV
    src%prop%area = MV
    src%prop%density = MV
    !Setting state variables
    src%now%age=0.0
    src%now%active=.false. !disabled by default
    src%now%emission_stack = 1
    src%now%pos=src%par%geometry%pt !coords of the Source (meaning depends on the geometry type!)
    !setting statistical samplers
    src%stats%particles_emitted=0
    src%stats%ns=0
    !setting stencil variables
    src%stencil%dp = Globals%SimDefs%Dp
    if (res%x > 0.0) src%stencil%dp = res !the source has a custom resolution
    allocate(src%stencil%ptlist, source = Geometry%getFillPoints(src%par%geometry, src%stencil%dp))
    src%stencil%np = size(src%stencil%ptlist)
    call src%setotalnp()
    src%stencil%ptlist = Utils%m2geo(src%stencil%ptlist)

    sizem = sizeof(src)
    call SimMemory%addsource(sizem)
    call src%print()

    end subroutine initializeSource

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Returns particulate status of this Source, i.e, true if the emitted
    !> Tracers are actually a collection of particles with an evolving
    !> concentration
    !---------------------------------------------------------------------------
    logical function isParticulate(self)
    class(source_class) :: self
    isParticulate = self%prop%particulate
    end function isParticulate

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that sets the total number of tracers a source will potentially create
    !> \f${NP}_{total}^{source-i}=int((T_{end}^{source-i}-T_{start}^{source-i})/(Dt/{Rate}^{source-i})*{NP}_{emission}^{source-i})\f$
    !---------------------------------------------------------------------------
    subroutine setotalnp(self)
    implicit none
    class(source_class), intent(inout) :: self
    self%stencil%total_np=int(count(self%par%activeTime)*Globals%SimDefs%dt*self%par%emitting_rate*self%stencil%np)
    end subroutine setotalnp

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> source print routine - prints a source info on console/log
    !> @param[in] src
    !---------------------------------------------------------------------------
    subroutine printSource(src)
    implicit none
    class(source_class) :: src
    type(string) :: outext
    type(string) :: temp_str(3)
    temp_str(1)=src%par%id
    outext = '-->Source '//src%par%name//new_line('a')//&
        '       Id = '//temp_str(1)//new_line('a')//&
        '       Geometry type is '//src%par%source_geometry//new_line('a')
    temp_str(1)=src%now%pos%x
    temp_str(2)=src%now%pos%y
    temp_str(3)=src%now%pos%z
    outext = outext//'       Initially at coordinates'//new_line('a')//&
        '       '//temp_str(1)//' '//temp_str(2)//' '//temp_str(3)//new_line('a')
    temp_str(1)=src%stencil%np
    if (src%par%emitting_fixed_rate) then
        temp_str(2)=src%par%emitting_rate
        outext = outext//'       Emitting '//temp_str(1)//' tracers at '//temp_str(2)//' Hz'
    else
        outext = outext//'       Emitting '//temp_str(1)//' tracers at a rate defined in '//src%par%rate_file//new_line('a')
        temp_str(2) = src%par%rateScale
        outext = outext//'       scaled by '//temp_str(2)
    end if    
    call Log%put(outext,.false.)
    end subroutine printSource

    end module sources_mod
