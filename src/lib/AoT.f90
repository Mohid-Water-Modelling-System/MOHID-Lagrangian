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

    module AoT_mod

    use tracers_mod
    use tracer_array_mod
    use tracer_list_mod
    use common_modules

    implicit none
    private

    type :: aot_class !< Arrays of Tracers class
        integer, allocatable, dimension(:) :: id        !< Id of the Tracer
        integer, allocatable, dimension(:) :: index     !< index of the Tracer
        real(prec), allocatable, dimension(:) :: x,y,z  !< coordinates of the Tracer
        real(prec), allocatable, dimension(:) :: u,v,w  !< velocities of the Tracer
    contains
    !sort
    procedure :: Clean
    procedure :: print => print_AoT
    end type aot_class

    interface AoT !> Constructor
        procedure constructor
    end interface

    public :: aot_class, AoT
    
    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Constructor for AoT object with data from a tracerList_class object
    !> @parm[in] trclist
    !---------------------------------------------------------------------------
    function constructor(trclist)
    implicit none
    type(aot_class) :: constructor
    class(tracerList_class), intent(in) :: trclist
    integer :: nt, i
    class(*), pointer :: aTracer
    type(string) :: outext
    !allocating the necessary space
    nt = trclist%getSize()
    allocate(constructor%id(nt))
    !allocate(constructor%index(nt))
    allocate(constructor%x(nt))
    allocate(constructor%y(nt))
    allocate(constructor%z(nt))
    allocate(constructor%u(nt))
    allocate(constructor%v(nt))
    allocate(constructor%w(nt))
    nt=1
    call trclist%reset()               ! reset list iterator
    do while(trclist%moreValues())     ! loop while there are values
        aTracer => trclist%currentValue() ! get current value
        select type(aTracer)
        class is (tracer_class)
            if (aTracer%now%active) then
                constructor%id(nt) = aTracer%par%id
                !constructor%index(nt) = i
                constructor%x(nt) = aTracer%now%pos%x
                constructor%y(nt) = aTracer%now%pos%y
                constructor%z(nt) = aTracer%now%pos%z
                constructor%u(nt) = aTracer%now%vel%x
                constructor%v(nt) = aTracer%now%vel%y
                constructor%w(nt) = aTracer%now%vel%z
                nt= nt + 1
            end if
            class default
            outext = '[AoT::Constructor]: Unexepected type of content, not a Tracer'
            call Log%put(outext)
            stop
        end select
        call trclist%next()            ! increment the list iterator
    end do
    call trclist%reset()               ! reset list iterator
    
    end function constructor
    
     !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Destructor for AoT object, deallocates all contents
    !---------------------------------------------------------------------------
    subroutine Clean(self)
    implicit none
    class(aot_class), intent(inout) :: self
    deallocate(self%id)
    !deallocate(self%index)
    deallocate(self%x)
    deallocate(self%y)
    deallocate(self%z)
    deallocate(self%u)
    deallocate(self%v)
    deallocate(self%w)
    end subroutine Clean

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the elements of the array
    !---------------------------------------------------------------------------
    subroutine print_AoT(self)
    class(aot_class), intent(in) :: self
    type(string) :: outext, t(7)
    integer :: i
    do i=1, size(self%id)
        t(1) = self%id(i)
        t(2) = self%x(i)
        t(3) = self%y(i)
        t(4) = self%z(i)
        outext = 'Tracer['//t(1)//']::xyz('//t(2)//','//t(3)//','//t(4)//')'
        call Log%put(outext,.false.)
    end do
    end subroutine print_AoT


    end module AoT_mod
