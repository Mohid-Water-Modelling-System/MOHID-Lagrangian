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

    use commom_modules
    use sources_mod
    use tracers_mod

    implicit none
    private

    type :: emitter_class
        integer :: emitted
        integer :: emittable
    contains
    procedure :: initialize => initializeEmitter
    procedure :: addSource
    procedure :: alloctracers
    procedure :: initracers
    !procedure :: activecheck
    end type

    !Public access vars
    public :: emitter_class

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that calls the tracer initialization from the emmiter object
    !
    !> @param[in] self, src
    !---------------------------------------------------------------------------
    subroutine initracers(self, srcs)
    implicit none
    class(emitter_class), intent(inout) :: self
    class(source_class), dimension(:), intent(inout) :: srcs
    integer num_emiss, i, j, k, p
    type(string) :: outext, temp(4)
    integer :: sizem

    p=0
    do i=1, size(srcs)
        num_emiss = srcs(i)%stencil%total_np/size(srcs(i)%stencil%ptlist)
        do j=1, num_emiss
            do k=1, size(srcs(i)%stencil%ptlist)
                p=p+1
                call Tracer(p)%initialize(p, srcs(i)%par%id, Globals%SimTime, srcs(i)%stencil%ptlist(k))
            enddo
        enddo
    enddo
    sizem = sizeof(Tracer)
    call SimMemory%addtracer(sizem)

    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that allocates the tracers respective to a given source
    !
    !> @param[in] self, src
    !---------------------------------------------------------------------------
    subroutine alloctracers(self, src)
    implicit none
    class(emitter_class), intent(inout) :: self
    class(source_class), intent(inout) :: src
    integer err
    type(string) :: outext, temp

    if (self%emittable .le. 0) then
        outext='[Emitter::alloctracers]: No Tracers will be simulated, stoping'
        call Log%put(outext)
        stop
    else
        allocate(Tracer(self%emittable), stat=err)
        if(err/=0)then
            outext='[Emitter::alloctracers]: Cannot allocate Tracers, stoping'
            call Log%put(outext)
            stop
        endif
    endif

    temp = size(Tracer)
    outext='Allocated '// temp // ' Tracers.'
    call Log%put(outext)
    !receiving Sources as argument so latter we can differentiate between tracer types

    end subroutine

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that initializes an emmiter class object. Sets default values
    !
    !> @param[in] self
    !---------------------------------------------------------------------------
    subroutine initializeEmitter(self)
    implicit none
    class(emitter_class), intent(inout) :: self
    self%emitted = 0
    self%emittable = 0
    end subroutine initializeEmitter
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to compute the total emittable particles per source and allocate 
    !> them
    !
    !> @param[in] self, src
    !---------------------------------------------------------------------------
    subroutine addSource(self, src)
    implicit none
    class(emitter_class), intent(inout) :: self
    class(source_class),intent(inout) :: src
    integer :: i
    
    call setotalnp(src) !finding the total tracers this Source will pass the emmiter
    self%emittable = self%emittable + src%stencil%total_np
        !print*, srcs(i)%stencil%total_np
   
    !allocating and initializing the tracers by the emitter, for all sources
    call self%alloctracers(src)
    !call self%initracers(srcs)

    end subroutine addSource

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> private routine that returns the total number of tracers an input
    !> source will potentially create
    !
    !> @param[in] src
    !---------------------------------------------------------------------------
    subroutine setotalnp(src)
    implicit none
    class(source_class), intent(inout) :: src
    !> \f${NP}_{total}^{source-i}=(T_{end}^{source-i}-T_{start}^{source-i})*{Rate}^{source-i}*{NP}_{emission}^{source-i}\f$
    src%stencil%total_np=(src%par%stoptime-src%par%startime)*src%par%emitting_rate*src%stencil%np
    end subroutine

  end module emitter_mod
