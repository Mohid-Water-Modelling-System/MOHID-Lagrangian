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
    use tracer_list_mod
    use sources_list_mod
    use common_modules

    implicit none
    private

    type :: emitter_class       !< Emitter class
        integer :: emitted      !< number of Tracers this Emitter has created
        integer :: emittable    !< number of Tracers this Emitter should create throughout the simulation
    contains
    procedure :: initialize => initializeEmitter
    procedure :: addSource
    procedure :: removeSource
    procedure :: emitt
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
    self%emittable = 0
    end subroutine initializeEmitter

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to compute the total emittable particles per source and allocate
    !> that space in the Blocks Tracer array
    !> @param[in] self, src
    !---------------------------------------------------------------------------
    subroutine addSource(self, src)
    implicit none
    class(emitter_class), intent(inout) :: self
    class(source_class),intent(in) :: src
    self%emittable = self%emittable + src%stencil%total_np
    end subroutine addSource

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method to remove from the total emittable particles count a Source
    !> @param[in] self, src
    !---------------------------------------------------------------------------
    subroutine removeSource(self, src)
    implicit none
    class(emitter_class), intent(inout) :: self
    class(source_class),intent(in) :: src
    self%emittable = self%emittable - src%stencil%total_np
    end subroutine removeSource

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that emitts the Tracers, based on the Sources on this Block Emitter
    !> @param[in] self, src, trclist
    !---------------------------------------------------------------------------
    subroutine emitt(self, src, trclist)
    implicit none
    class(emitter_class), intent(inout) :: self !> the Emmiter from the Block where the Source is
    class(source_class), intent(inout)  :: src  !>the Source that will emitt new Tracers
    class(tracerList_class), intent(inout)   :: trclist  !>the Tracer list from the Block where the Source is
    integer i
    class(*), allocatable :: newtrc
    do i=1, src%stencil%np
        self%emitted = self%emitted + 1
        !PARALLEL The calls inside this routine MUST be atomic in order to get the correct sequencial Tracer Id
        call self%tracerMaker(newtrc, src, i)
        call trclist%add(newtrc)
    end do
    src%stats%particles_emitted = src%stats%particles_emitted + src%stencil%np
    end subroutine emitt

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
    select case (src%prop%property_type%chars())
    case ('base')
        allocate(trc, source = Tracer(Globals%Sim%getnumTracer(), src, Globals%SimTime, p)) !Beacause ifort is not F2008 compliant.
        !trc = Tracer(1, src, Globals%SimTime, p) !Otherwise instinsic allocation would be enough and more readable, like this. Compiles fine in GFortran
    case ('paper')
        allocate(trc, source = paperTracer(Globals%Sim%getnumTracer(), src, Globals%SimTime, p))
    case ('plastic')
        allocate(trc, source = Tracer(Globals%Sim%getnumTracer(), src, Globals%SimTime, p))
        case default
        outext='[Emitter::tracerMaker]: unexpected type for Tracer object: '//src%prop%property_type
        call Log%put(outext)
        stop
    end select

    end subroutine tracerMaker


    end module emitter_mod
