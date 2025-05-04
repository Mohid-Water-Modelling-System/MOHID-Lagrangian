    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : stateVector_mod
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the State Vector class and its methods. This class
    !> defines a state vector that encodes the state of a type of Tracers, 
    !> allowing for trivial application of kernels and subsequent integration. 
    !> These must be exported back into the objects from this class.
    !------------------------------------------------------------------------------

    module stateVector_mod

    use tracers_mod
    use tracerList_mod
    use common_modules

    implicit none
    private

    type :: trcPtr_class                   !< tracer pointer class, because foooooortraaaaaaan
        class(tracer_class), pointer :: ptr => null() !< the actual pointer
    end type trcPtr_class

    type :: stateVector_class
        integer :: ttype
        type(trcPtr_class), allocatable, dimension(:) :: trc   !< pointer to the Tracer
        real(prec), allocatable, dimension(:,:) :: state
        type(string), allocatable, dimension(:) :: varName
        integer, allocatable, dimension(:) :: source, id
        real(prec) , allocatable, dimension(:) :: landIntMask
        real(prec) , allocatable, dimension(:) :: resolution
        logical, allocatable, dimension(:) :: active
        integer :: idx
    contains
    procedure :: toTracers
    procedure :: copyState
    procedure :: finalize => cleanState
    end type stateVector_class

    public :: stateVector_class

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Copies a State Vector to another.
    !---------------------------------------------------------------------------
    subroutine copyState(self, newsv)
    class(stateVector_class), intent(in) :: self
    type(stateVector_class), intent(out) :: newsv
    
    newsv%ttype = self%ttype
    allocate(newsv%state(size(self%state,1), size(self%state,2)))
    newsv%state = self%state    
    allocate(newsv%landIntMask(size(self%landIntMask)))
    newsv%landIntMask = self%landIntMask
    allocate(newsv%resolution(size(self%resolution)))
    newsv%resolution = self%resolution
    allocate(newsv%active(size(self%active)))
    newsv%active = self%active
    allocate(newsv%varName(size(self%varName)))
    newsv%varName = self%varName
    allocate(newsv%source(size(self%source)))
    newsv%source = self%source  
    !maybe no need to copy source, id and tracer pointer
    
    end subroutine copyState

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Sends the data on the State Vector to the Tracer objects.
    !---------------------------------------------------------------------------
    subroutine toTracers(self)
    class(stateVector_class), intent(in) :: self
    integer :: i
    class(tracer_class), pointer :: aTracer
    type(string) :: outext
    if(allocated(self%active)) then
        do i=1, size(self%active)
            if (associated(self%trc(i)%ptr)) then
                aTracer => self%trc(i)%ptr
                if (self%active(i)) then
                    aTracer%now%active = .true.
                    call aTracer%setStateArray(self%state(i,:))
                end if
            else
                outext = '[stateVector::toTracers]: pointer to Tracer not associated, stoping'
                call Log%put(outext)
                stop
            end if
        end do
    end if
    end subroutine toTracers

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Cleans the State Vector object
    !---------------------------------------------------------------------------
    subroutine cleanState(self)
    class(stateVector_class), intent(inout) :: self
    self%ttype = MV
    if (allocated(self%trc)) deallocate(self%trc)
    if (allocated(self%state)) deallocate(self%state)
    if (allocated(self%landIntMask)) deallocate(self%landIntMask)
    if (allocated(self%resolution)) deallocate(self%resolution)
    if (allocated(self%active)) deallocate(self%active)
    if (allocated(self%source)) deallocate(self%source)
    if (allocated(self%id)) deallocate(self%id)
    if (allocated(self%varName)) deallocate(self%varName)
    self%idx = 1
    end subroutine cleanState

    end module stateVector_mod