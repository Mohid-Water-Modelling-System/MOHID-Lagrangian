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

    use common_modules
    use sources_mod
    use tracers_mod
    use tracer_array_mod
    use sources_array_mod    

    implicit none
    private

    type :: emitter_class
        integer :: id
        integer :: emitted
        integer :: emittable
    contains
    procedure :: initialize => initializeEmitter
    procedure :: addSource
    procedure :: removeSource
    procedure :: emitt
    procedure :: tracerMaker
    end type

    !Public access vars
    public :: emitter_class

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that initializes an emmiter class object. Sets default values
    !
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine initializeEmitter(self, id)
    implicit none
    class(emitter_class), intent(inout) :: self
    integer, intent(in) :: id
    self%id = id
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
    !> this method returns a resized Tracer array if needed to the corresponding
    !> Block
    !> @param[in] self, src, trc
    !---------------------------------------------------------------------------
    subroutine emitt(self, src, trcarr)
    implicit none
    class(emitter_class), intent(inout) :: self !> the Emmiter from the Block where the Source is
    class(source_class), intent(inout)  :: src  !>the Source that will emitt new Tracers
    class(TracerArray), intent(inout)   :: trcarr  !>the Tracer array from the Block where the Source is
    integer err, i
    type(string) :: outext, temp(2)
    integer :: allocstride = 3
    class(*), allocatable :: newtrc

    if (self%emittable <= 0) then
        !nothing to do as we have no Sources or no emittable Tracers
        temp(1) = self%id
        temp(2) = src%par%id
        outext='-->Source '//temp(2)//' trying to emitt Tracers from an exausted Emitter '//temp(1)
        call Log%put(outext,.false.)
    else
        !check if the Block Tracer Array has enough free places for this emission
        if (src%stencil%np > (trcarr%getLength() - trcarr%lastActive)) then
            call trcarr%resize(trcarr%getLength() + allocstride*src%stencil%np, initvalue = dummyTracer) !resizing the Block Tracer array to accomodate more emissions
        end if
        !there is space to emmitt the Tracers
        do i=1, src%stencil%np
            self%emitted = self%emitted + 1 
            self%emittable = self%emittable - 1
            trcarr%lastActive = trcarr%lastActive + 1
            trcarr%numActive = trcarr%numActive + 1
            !PARALLEL The calls inside this routine MUST be atomic in order to get the correct sequencial Tracer Id
            call self%tracerMaker(newtrc, src, i)
            call trcarr%put(trcarr%lastActive, newtrc)
        end do
        src%stats%particles_emitted = src%stats%particles_emitted + src%stencil%np
    endif

    end subroutine
    
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> method that calls the corresponding Tracer constructor, depending on the
    !> requested type from the emitting Source
    !> @param[in] sself, trc, src, p
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
