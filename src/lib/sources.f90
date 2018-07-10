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
    use simulation_globals_mod

    implicit none
    private

    type :: source_par               !<Type - parameters of a source object
        integer :: id                       !< unique source identification (integer)
        integer :: emitting_rate    !< Emitting rate of the source (Hz)
        real(prec_time) :: startime         !< time to start emitting tracers
        real(prec_time) :: stoptime         !< time to stop emitting tracers
        type(string) :: name                !< source name
        type(string) :: property_type       !< source property type (plastic, paper, fish, etc)
        type(string) :: property_name       !< source property name
        type(string) :: source_geometry     !< Source type : 'point', 'line', 'sphere', 'box'
        class(shape), allocatable :: geometry   !< Source geometry
    end type

    type :: source_state             !<Type - state variables of a source object
        real(prec_time) :: age              ! time variables
        logical :: active                   !< active switch
        integer :: emission_stride          !< Number of time steps to wait until next emission
        type(vector) :: pos                 !< Position of the source baricenter (m)
        type(vector) :: vel                 !< Velocity of the source (m s-1)
        real(prec) :: depth                 !< Depth of the source baricenter (m)
        real(prec) :: T                     !< Temperature of the source (Celcius)
    end type

    type :: source_stats             !<Type - statistical variables of a source object
        ! All stats variables at writing precision (prec_wrt)
        ! Avegarge variable is computed by Accumulated_var / ns
        integer :: particles_emitted        !< Number of emitted particles by this source
        real(prec_wrt) :: acc_T             !< Accumulated temperature of the tracer (Celcius)
        integer :: ns                       !< Number of sampling steps
    end type

    type :: source_stencil         !<Type - holder for the tracer creation stencil of the source
        integer :: np                       !< Number of tracers by emission
        integer :: total_np                 !< Total number of tracers that this source will generate
        type(vector), allocatable, dimension(:) :: ptlist !<list of points (coordinates), relative to the source geometry point, to be generated at every emission
    end type

    type :: source_class           !<Type - The source class
        type(source_par)   :: par           !<To access parameters
        type(source_state) :: now           !<To access state variables
        type(source_stencil) :: stencil     !<To acess stencil variables
        type(source_stats) :: stats         !<To access statistics
    contains
    procedure :: initialize => initializeSource
    procedure :: setotalnp
    procedure :: linkproperty
    procedure :: print => printSource
    end type

    type :: source_group_class
        type(source_class), allocatable, dimension(:) :: src
    contains
    procedure :: initialize => initSources
    procedure :: finalize => killSources
    procedure :: setProps
    end type

    !Simulation variables
    type(source_group_class) :: tempSources

    !Public access vars
    public :: tempSources, source_group_class, source_class
    !Public access procedures

    contains

    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source allocation routine - allocates sources objects
    !
    !> @param[in] nsources
    !---------------------------------------------------------------------------
    subroutine initSources(self,nsources)
    implicit none
    class(source_group_class), intent(inout) :: self
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source group destructor - deallocates sources objects
    !---------------------------------------------------------------------------
    subroutine killSources(self)
    implicit none
    class(source_group_class), intent(inout) :: self
    integer err
    type(string) :: outext
    if (ALLOCATED(self%src)) deallocate(self%src, stat=err)    
    if(err/=0)then
        outext='[killSources]: Cannot deallocate Sources, stoping'
        call Log%put(outext)
        stop    
    endif
    end subroutine killSources


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source property setting routine, calls source by id to set its properties
    !
    !> @param[in] srcid,ptype,pname
    !---------------------------------------------------------------------------
    subroutine setProps(self,srcid_str,ptype,pname)
    implicit none
    class(source_group_class), intent(inout) :: self
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
            call self%src(i)%linkproperty(ptype,pname) ! calling Source method to link property
            temp = self%src(i)%par%id
            outext='      Source id = '// temp // ', '// self%src(i)%par%name //' is of type '// self%src(i)%par%property_type //', with property name ' // self%src(i)%par%property_name
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
    end subroutine

    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source inititialization proceadure - initializes Source variables
    !
    !> @param[in] src, id, name, emitting_rate, source_geometry
    !---------------------------------------------------------------------------
    subroutine initializeSource(src,id,name,emitting_rate,start,finish,source_geometry,shapetype)
    implicit none
    class(source_class) :: src
    integer, intent(in) :: id
    type(string), intent(in) :: name
    real(prec), intent(in) :: emitting_rate
    real(prec), intent(in) :: start
    real(prec), intent(in) :: finish
    type(string), intent(in) :: source_geometry
    class(shape), intent(in) :: shapetype

    integer :: sizem, i
    type(string) :: outext
    integer :: err

    !Setting parameters
    src%par%id=id
    src%par%emitting_rate=emitting_rate
    src%par%startime=start
    src%par%stoptime=finish
    src%par%name=name
    src%par%source_geometry=source_geometry
    allocate(src%par%geometry, source=shapetype)
    src%par%property_type = "base" ! pure Lagrangian trackers by default
    src%par%property_name = "base"
    !Setting state variables
    src%now%age=0.0
    src%now%active=.false. !disabled by default
    src%now%emission_stride=1 !first time-step once active the Source emitts
    src%now%pos=src%par%geometry%pt !coords of the Source (meaning depends on the geometry type!)
    !setting statistical samplers
    src%stats%particles_emitted=0
    src%stats%acc_T=0.0
    src%stats%ns=0
    !setting stencil variables
    src%stencil%np = Geometry%fillsize(src%par%geometry)
    call src%setotalnp()
    allocate(src%stencil%ptlist(src%stencil%np), stat=err)
    if(err/=0)then
        outext='Cannot allocate point list for Source '// src%par%name //', stoping'
        call Log%put(outext)
        stop
    endif
    call Geometry%fill(src%par%geometry, src%stencil%np, src%stencil%ptlist)

    sizem = sizeof(src)
    call SimMemory%addsource(sizem)
    call src%print()

    !DBG
    !do i=1, src%stencil%np
    !print*, src%stencil%ptlist(i)
    !end do
    end subroutine initializeSource


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source property setting proceadure - initializes Source variables
    !
    !> @param[in] src,ptype,pname
    !---------------------------------------------------------------------------
    subroutine linkproperty(src,ptype,pname)
    implicit none
    class(source_class) :: src
    type(string), intent(in) :: ptype
    type(string), intent(in) :: pname
    src%par%property_type = ptype
    src%par%property_name = pname
    end subroutine linkproperty

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that sets the total number of tracers a source will potentially create

    !---------------------------------------------------------------------------
    subroutine setotalnp(self)
        implicit none
        class(source_class), intent(inout) :: self
        !< computing the total as \f${NP}_{total}^{source-i}=int((T_{end}^{source-i}-T_{start}^{source-i})/(Dt/{Rate}^{source-i})*{NP}_{emission}^{source-i})\f$
        self%stencil%total_np=int((self%par%stoptime-self%par%startime)/(Globals%SimDefs%dt)/self%par%emitting_rate*self%stencil%np)
    end subroutine setotalnp

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source print routine - prints a source info on console/log
    !
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
    temp_str(1)=src%par%emitting_rate
    temp_str(2)=src%stencil%np
    temp_str(3)=src%stencil%total_np
    outext = outext//'       Emitting '//temp_str(2)//' tracers at every '//temp_str(1)//' time-steps'//new_line('a')
    outext = outext//'       For an estimated total of '//temp_str(3)//' tracers' //new_line('a')
    temp_str(1)=src%par%startime
    temp_str(2)=src%par%stoptime
    outext = outext//'       Active from '//temp_str(1)//' to '//temp_str(2)//' seconds'

    call Log%put(outext,.false.)

    end subroutine printSource

    end module sources_mod
