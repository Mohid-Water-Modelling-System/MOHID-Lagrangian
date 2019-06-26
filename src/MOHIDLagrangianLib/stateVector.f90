    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : AoT
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the Arrays of Tracers class and its methods. This class
    !> defines a collection of id, xyz, uvw, .. arrays that allow for easy and
    !> efficient manipulation of the Tracer objects. These must be exported
    !> into the objects from this class
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
        integer, allocatable, dimension(:) :: landMask, landIntMask, source, id
        logical, allocatable, dimension(:) :: active
        integer :: idx
    contains
    procedure :: toTracers
    procedure :: finalize => cleanState
    end type stateVector_class

    public :: stateVector_class

    contains
    
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
    if (allocated(self%landMask)) deallocate(self%landMask)
    if (allocated(self%landIntMask)) deallocate(self%landIntMask)
    if (allocated(self%active)) deallocate(self%active)
    if (allocated(self%source)) deallocate(self%source)
    if (allocated(self%id)) deallocate(self%id)
    self%idx = 1    
    end subroutine cleanState
    
    end module stateVector_mod