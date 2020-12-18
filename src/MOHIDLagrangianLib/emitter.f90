    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : emitter
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : April 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines an emitter class and related methods. This module is
    !> responsible for building a potential tracer list based on the availble
    !> sources and calling their initializers.
    !------------------------------------------------------------------------------

    module emitter_mod

    use sources_mod
    use tracers_mod
    use tracerList_mod
    use sourcesList_mod
    use common_modules

    implicit none
    private

    type :: emitter_class       !< Emitter class
        integer :: emitted      !< number of Tracers this Emitter has created
    contains
    procedure :: initialize => initializeEmitter
    procedure :: emitt
    procedure :: emitt_src
    procedure :: tracerMaker
    end type emitter_class

    !Public access vars
    public :: emitter_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that initializes an emmiter class object. Sets default values
    !---------------------------------------------------------------------------
    subroutine initializeEmitter(self)
    implicit none
    class(emitter_class), intent(inout) :: self
    self%emitted = 0
    end subroutine initializeEmitter

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that emitts the Tracers, based on the Sources on the list of the
    !> Emitter's Block
    !> @param[in] self, srclist, trclist
    !---------------------------------------------------------------------------
    subroutine emitt(self, srclist, trclist)
    implicit none
    class(emitter_class), intent(inout) :: self          !> the Emmiter from the Block where the Source is
    class(sourceList_class), intent(inout)  :: srclist   !>the Source that will emitt new Tracers
    class(tracerList_class), intent(inout)  :: trclist   !>the Tracer list from the Block where the Source is
    class(*), pointer :: aSource
    type(string) :: outext
    integer :: i
    logical :: reset_stack

    reset_stack = .false.
    call srclist%reset()                   ! reset list iterator
    do while(srclist%moreValues())         ! loop while there are values
        aSource => srclist%currentValue()  ! get current value
        select type(aSource)
        class is (source_class)
            if (.not.aSource%par%fixed_position) call aSource%setVariablePosition(Globals%Sim%getnumdt())
            if (aSource%now%active) then
                if (.not.aSource%par%emitting_fixed_rate) call aSource%setVariableRate(Globals%Sim%getnumdt())                
                aSource%now%emission_stack = aSource%now%emission_stack + aSource%par%emitting_rate*Globals%SimDefs%dt  !adding to the emission stack
                do i=1, floor(aSource%now%emission_stack)
                    call self%emitt_src(aSource, trclist)
                    reset_stack = .true.                    
                end do 
                if (reset_stack) then
                    aSource%now%emission_stack = 0 !reseting for the next time step              
                end if
            end if
            class default
            outext = '[Emitter] Unexepected type of content, not a Source'
            call Log%put(outext)
            stop
        end select
        call srclist%next()            ! increment the list iterator
    end do
    call srclist%reset()               ! reset list iterator

    end subroutine emitt

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that emitts the Tracers, given a particular Source
    !> @param[in] self, src, trclist
    !---------------------------------------------------------------------------
    subroutine emitt_src(self, src, trclist)
    implicit none
    class(emitter_class), intent(inout) :: self !> the Emmiter from the Block where the Source is
    class(source_class), intent(inout)  :: src  !>the Source that will emitt new Tracers
    class(tracerList_class), intent(inout)   :: trclist  !>the Tracer list from the Block where the Source is
    integer i
    class(*), allocatable :: newtrc
    do i=1, src%stencil%np
        !PARALLEL The calls inside this routine MUST be atomic in order to get the correct sequencial Tracer Id
        call self%tracerMaker(newtrc, src, i)
        call trclist%add(newtrc)
    end do
    self%emitted = self%emitted + src%stencil%np
    src%stats%particles_emitted = src%stats%particles_emitted + src%stencil%np
    end subroutine emitt_src

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that calls the corresponding Tracer constructor, depending on the
    !> requested type from the emitting Source
    !> @param[in] self, trc, src, p
    !---------------------------------------------------------------------------
    subroutine tracerMaker(self, trc, src, p)
    implicit none
    class(emitter_class), intent(in) :: self
    class(*), allocatable, intent(out) :: trc
    class(source_class), intent(in) :: src
    integer, intent(in) :: p
    type(string) :: outext, temp

    !PARALLEL Globals%Sim%getnumTracer() MUST be atomic in order to get the correct sequencial Tracer Id
    select case (src%prop%propertyType%chars())
    case ('base')
        allocate(trc, source = Tracer(Globals%Sim%getnumTracer(), src, Globals%SimTime%CurrTime, p)) !Beacause ifort 2017 is not F2008 compliant...
        !trc = Tracer(1, src, Globals%Time%CurrTime, p) !Otherwise instinsic allocation would be enough and more readable, like this. Compiles fine in GFortran
    case ('paper')
        allocate(trc, source = paperTracer(Globals%Sim%getnumTracer(), src, Globals%SimTime%CurrTime, p))
    case ('plastic')
        allocate(trc, source = plasticTracer(Globals%Sim%getnumTracer(), src, Globals%SimTime%CurrTime, p))
    case ('coliform')
        allocate(trc, source = coliformTracer(Globals%Sim%getnumTracer(), src, Globals%SimTime%CurrTime, p))
        case default
        outext='[Emitter::tracerMaker]: unexpected type for Tracer object: '//src%prop%propertyType
        call Log%put(outext)
        stop
    end select

    end subroutine tracerMaker


    end module emitter_mod
