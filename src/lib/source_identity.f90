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

    implicit none
    private

    type source_par_class               !<Type - parameters of a source object
        integer :: id                       !< unique source identification (integer)
        real(prec) :: emitting_rate         !< Emitting rate of the source (Hz)
        type(string) :: name                !< source name
        type(string) :: property_name       !< source property name
        type(string) :: source_geometry     !< Source type : 'point', 'line', 'sphere', 'box'
        class(shape), allocatable :: geometry   !< Source geometry
    end type

    type source_state_class             !<Type - state variables of a source object
        real(prec_time) :: age              ! time variables
        logical :: active                   !< active switch
        type(vector) :: pos                 !< Position of the source baricenter (m)
        type(vector) :: vel                 !< Velocity of the source (m s-1)
        real(prec) :: depth                 !< Depth of the source baricenter (m)
        real(prec) :: T                     !< Temperature of the source (Celcius)
    end type

    type source_stats_class             !<Type - statistical variables of a source object
        ! All stats variables at writing precision (prec_wrt)
        ! Avegarge variable is computed by Accumulated_var / ns
        integer ::particles_emitted         !< Number of emitted particles by this source
        real(prec_wrt) :: acc_T             !< Accumulated temperature of the tracer (Celcius)
        integer :: ns                       !< Number of sampling steps
    end type

    type source_class                   !<Type - The source class
        type(source_par_class)   :: par     !<To access parameters
        type(source_state_class) :: now     !<To access state variables
        type(source_stats_class) :: stats   !<To access statistics
    end type

    !Simulation variables
    type(source_class), allocatable, dimension(:) :: Source

    
    !Public access vars
    public :: Source
    !Public access procedures
    public :: initSource, allocSources

    contains
    
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
        if(err==0)then
            !print*, nsources, " sources allocated"
        else
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
    !> source inititialization routine - Generates a source and initializes its variables
    !
    !> @param[out] source
    !> @param[in] num, id, name, emitting_rate, source_geometry
    !---------------------------------------------------------------------------
    subroutine initSource(num,id,name,emitting_rate,source_geometry,geometry)
    implicit none
    integer, intent(in) :: num
    integer, intent(in) :: id
    type(string), intent(in) :: name
    real(prec), intent(in) :: emitting_rate
    type(string), intent(in) :: source_geometry
    class(shape), intent(in) :: geometry
    
    !Setting parameters
    Source(num)%par%id=id
    Source(num)%par%emitting_rate=emitting_rate
    Source(num)%par%name=name
    Source(num)%par%source_geometry=source_geometry
    allocate(Source(num)%par%geometry, source=geometry)
    !Setting state variables
    Source(num)%now%age=0.0
    Source(num)%now%active=.false. !disabled by default    
    Source(num)%now%pos=Source(num)%par%geometry%pt !coords of the Source (meaning depends on the geometry type!)
    !setting statistical samplers
    Source(num)%stats%particles_emitted=0
    Source(num)%stats%acc_T=0.0
    Source(num)%stats%ns=0  
    
    call printsource(Source(num))
    
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
    subroutine printsource(src)
    implicit none
    type(source_class), intent(in) :: src
    
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
    outext = outext//'       Emitting rate of '//temp_str(1)//' Hz'
    
    call ToLog(outext,.false.)    
    
    end subroutine    

    end module source_identity
