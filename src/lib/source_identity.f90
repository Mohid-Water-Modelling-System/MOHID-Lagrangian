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

    module source_identity

    use commom_modules
    use simulation_globals

    implicit none
    private

    type source_par               !<Type - parameters of a source object
        integer :: id                       !< unique source identification (integer)
        real(prec_time) :: emitting_rate    !< Emitting rate of the source (Hz)
        real(prec_time) :: starttime        !< time to start emitting tracers
        real(prec_time) :: stoptime         !< time to stop emitting tracers
        type(string) :: name                !< source name
        type(string) :: property_name       !< source property name
        type(string) :: source_geometry     !< Source type : 'point', 'line', 'sphere', 'box'
        class(shape), allocatable :: geometry   !< Source geometry
    end type

    type source_state             !<Type - state variables of a source object
        real(prec_time) :: age              ! time variables
        logical :: active                   !< active switch
        type(vector) :: pos                 !< Position of the source baricenter (m)
        type(vector) :: vel                 !< Velocity of the source (m s-1)
        real(prec) :: depth                 !< Depth of the source baricenter (m)
        real(prec) :: T                     !< Temperature of the source (Celcius)
    end type

    type source_stats             !<Type - statistical variables of a source object
        ! All stats variables at writing precision (prec_wrt)
        ! Avegarge variable is computed by Accumulated_var / ns
        integer :: particles_emitted        !< Number of emitted particles by this source
        real(prec_wrt) :: acc_T             !< Accumulated temperature of the tracer (Celcius)
        integer :: ns                       !< Number of sampling steps
    end type
    
    type source_stencil         !<Type - holder for the tracer creation stencil of the source
        integer :: np                       !<number of tracers by emission
        integer :: total_np                 !< Total number of tracers that this source will generate
        type(vector), allocatable :: ptlist !<list of points (coordinates), relative to the source geometry point, to be generated at every emission
    end type

    type source_class           !<Type - The source class
        type(source_par)   :: par           !<To access parameters
        type(source_state) :: now           !<To access state variables
        type(source_stencil) :: stencil     !<To acess stencil variables
        type(source_stats) :: stats         !<To access statistics
    contains
    procedure :: initialize
    procedure :: printout
    end type

    !Simulation variables
    type(source_class), allocatable, dimension(:) :: Source


    !Public access vars
    public :: Source
    !Public access procedures
    public :: allocSources

    contains
    
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Private source stencil initialization routine
    !
    !> @param[in] np, total_np, plist, geometry
    !---------------------------------------------------------------------------
    subroutine setnpstencil(np,total_np,geometry)
    implicit none
    integer, intent(inout) :: np
    integer, intent(inout) :: total_np
    class(shape), intent(in) :: geometry
    type(vector) :: temp
    type(string) :: outext
    
    select type (geometry)
    type is (shape)
        !nothing to do
    class is (box)
        np = int(geometry%size%x/SimDefs%Dp)*int(geometry%size%y/SimDefs%Dp)*int(geometry%size%z/SimDefs%Dp)        
    class is (point)
        np = 1
    class is (line)
        temp = geometry%pt-geometry%last !simply the difference
        np = int(temp%normL2()/SimDefs%Dp)
    class is (sphere)
        
        class default
        outext='setnpstencil: unexpected type for geometry object!'
        call ToLog(outext)
        stop
    end select
    
    end subroutine
    
    
    !!---------------------------------------------------------------------------
    !!> @Ricardo Birjukovs Canelas - MARETEC
    !! Routine Author Name and Affiliation.
    !!
    !!> @brief
    !!> Private source stencil initialization and allocation routine
    !!
    !!> @param[in] np, total_np, plist, geometry
    !!---------------------------------------------------------------------------
    !subroutine setsourcestencil(np,total_np,plist,geometry)
    !implicit none
    !integer, intent(inout) :: np
    !integer, intent(inout) :: total_np
    !
    !end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source allocation routine - allocates the sources objects
    !
    !> @param[in] nsources
    !---------------------------------------------------------------------------
    subroutine allocSources(nsources)
    implicit none
    integer, intent(in) :: nsources
    integer err
    type(string) :: outext
    allocate(Source(nsources), stat=err)
    if(err/=0)then
        outext='Cannot allocate Sources, stoping'
        call ToLog(outext)
        stop
    endif
    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source inititialization proceadure - initializes Source variables
    !
    !> @param[in] src, id, name, emitting_rate, source_geometry
    !---------------------------------------------------------------------------
    subroutine initialize(src,id,name,emitting_rate,start,finish,source_geometry,geometry)
    implicit none
    class(source_class) :: src
    integer, intent(in) :: id
    type(string), intent(in) :: name
    real(prec), intent(in) :: emitting_rate
    real(prec), intent(in) :: start
    real(prec), intent(in) :: finish
    type(string), intent(in) :: source_geometry
    class(shape), intent(in) :: geometry
    integer :: size

    !Setting parameters
    src%par%id=id
    src%par%emitting_rate=emitting_rate
    src%par%starttime=start
    src%par%stoptime=finish
    src%par%name=name
    src%par%source_geometry=source_geometry
    allocate(src%par%geometry, source=geometry)
    !Setting state variables
    src%now%age=0.0
    src%now%active=.false. !disabled by default
    src%now%pos=src%par%geometry%pt !coords of the Source (meaning depends on the geometry type!)
    !setting statistical samplers
    src%stats%particles_emitted=0
    src%stats%acc_T=0.0
    src%stats%ns=0
    !setting stencil variables
    !call setsourcestencil(src%stencil%np,src%stencil%total_np,src%stencil%plist,src%par%geometry)
    
    size = sizeof(src)
    call SimMemory%addsource(size)
    call src%printout()

    end subroutine
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> source print routine - prints a source info on console/log
    !
    !> @param[in] src
    !---------------------------------------------------------------------------
    subroutine printout(src)
    implicit none
    class(source_class) :: src

    type(string) :: outext
    type(string) :: temp_str(3)

    temp_str(1)=src%par%id
    outext = '-->Source '//src%par%name//' allocated'//new_line('a')//&
        '       Id = '//temp_str(1)//new_line('a')//&
        '       Geometry type is '//src%par%source_geometry//new_line('a')
    temp_str(1)=src%now%pos%x
    temp_str(2)=src%now%pos%y
    temp_str(3)=src%now%pos%z
    outext = outext//'       Initially at coordinates'//new_line('a')//&
        '       '//temp_str(1)//new_line('a')//&
        '       '//temp_str(2)//new_line('a')//&
        '       '//temp_str(3)//new_line('a')
    temp_str(1)=src%par%emitting_rate
    outext = outext//'       Emitting rate of '//temp_str(1)//' Hz'//new_line('a')
    temp_str(1)=src%par%starttime
    temp_str(2)=src%par%stoptime
    outext = outext//'       Active from '//temp_str(1)//' to'//temp_str(2)//' seconds'

    call ToLog(outext,.false.)

    end subroutine

    end module source_identity
