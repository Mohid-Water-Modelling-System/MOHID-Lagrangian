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
        integer :: emitted
        integer :: emittable
    contains
    procedure :: initialize => initializeEmitter
    procedure :: addSource
    procedure :: removeSource
    procedure :: emitt
    !procedure :: initracers
    !procedure :: activecheck
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
    subroutine initializeEmitter(self)
    implicit none
    class(emitter_class), intent(inout) :: self
    self%emitted = 0
    self%emittable = 0
    end subroutine initializeEmitter
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to compute the total emittable particles per source and allocate 
    !> that space in the Blocks Tracer array
    !
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method to remove from the total emittable particles count a Source
    !
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
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> method that emitts the Tracers, based on the Sources on this Block Emitter
    !> this method returns a resized Tracer array if needed to the corresponding
    !> Block
    !
    !> @param[in] self, src, trc
    !---------------------------------------------------------------------------
    subroutine emitt(self, src, trcarr)
    implicit none
    class(emitter_class), intent(inout) :: self !> the Emmiter from the Block where the Source is
    class(source_class), intent(inout)  :: src  !>the Source that will emitt new Tracers
    class(TracerArray), intent(inout)   :: trcarr  !>the Tracer array from the Block where the Source is
    integer err, i
    type(string) :: outext, temp
    integer :: allocstride = 3
    class(*), pointer :: newtrc

    if (self%emittable <= 0) then
        !nothing to do as we have no Sources or no emittable Tracers
    else
        !check if the Block Tracer Array has enough free places for this emission
        if (src%stencil%np >= (trcarr%getLength() - trcarr%lastActive)) then
             print*, 'Required space  = ', src%stencil%np
             print*, 'Active tracers  = ', trcarr%lastActive
             print*, 'Array length    = ', trcarr%getLength()
             print*, 'Initial length  = ', src%stencil%total_np
             print*, 'Available space = ', trcarr%getLength() - trcarr%lastActive
            call trcarr%resize(trcarr%getLength() + allocstride*src%stencil%np, initvalue = dummyTracer) !resizing the Block Tracer array to accomodate more emissions
        end if
        !there is space to emmitt the Tracers
        do i=1, src%stencil%np
            trcarr%lastActive = trcarr%lastActive + 1 !will need to change to paralelize
            trcarr%numActive = trcarr%numActive + 1
            !newtrc = tracer(src,i) !calling the constructor for a tracer

        end do
    endif

    !temp = size(Tracer)
    !outext='Allocated '// temp // ' Tracers.'
    !call Log%put(outext)
    !!receiving Sources as argument so latter we can differentiate between tracer types

    end subroutine
    
    
    !!---------------------------------------------------------------------------
    !!> @author Ricardo Birjukovs Canelas - MARETEC
    !! Routine Author Name and Affiliation.
    !!
    !!> @brief
    !!> method that calls the tracer initialization from the emmiter object
    !!
    !!> @param[in] self, src
    !!---------------------------------------------------------------------------
    !subroutine initracers(self, srcs)
    !implicit none
    !class(emitter_class), intent(inout) :: self
    !class(source_class), dimension(:), intent(inout) :: srcs
    !integer num_emiss, i, j, k, p
    !type(string) :: outext, temp(4)
    !integer :: sizem
    !
    !p=0
    !do i=1, size(srcs)
    !    num_emiss = srcs(i)%stencil%total_np/size(srcs(i)%stencil%ptlist)
    !    do j=1, num_emiss
    !        do k=1, size(srcs(i)%stencil%ptlist)
    !            p=p+1
    !            call Tracer(p)%initialize(p, srcs(i)%par%id, Globals%SimTime, srcs(i)%stencil%ptlist(k))
    !        enddo
    !    enddo
    !enddo
    !sizem = sizeof(Tracer)
    !call SimMemory%addtracer(sizem)
    !
    !end subroutine
    

  end module emitter_mod
