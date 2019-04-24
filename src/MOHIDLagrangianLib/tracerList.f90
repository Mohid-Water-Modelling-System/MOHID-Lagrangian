    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : tracer_list
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the tracer linked list class and its methods. This class
    !> defines a double linked list to store any variable type, but with specific
    !> methods with type guards for Tracer objects. The class allows for insertion,
    !> deletion and iteration of the desired contents
    !------------------------------------------------------------------------------

    module tracerList_mod

    use tracers_mod
    use abstract_LinkedList_mod
    use common_modules

    implicit none
    private

    type, extends(linkedlist) :: tracerList_class !< List of Tracers class
        ! integer, private :: numActive = 0
    contains
    ! procedure :: getNumActive
    ! procedure :: lowerNumActive
    ! procedure :: addTrc
    procedure :: print => print_tracerList
    procedure :: printCurrent => print_tracerListCurrent
    end type tracerList_class

    public :: tracerList_class

    contains

    ! !---------------------------------------------------------------------------
    ! !> @author Ricardo Birjukovs Canelas - MARETEC
    ! !> @brief
    ! !> Method that returns the number of active Tracers on this list
    ! !---------------------------------------------------------------------------
    ! integer function getNumActive(this)
    ! class(tracerList_class), intent(in) :: this
    ! getNumActive = this%numActive
    ! end function getNumActive

    ! !---------------------------------------------------------------------------
    ! !> @author Ricardo Birjukovs Canelas - MARETEC
    ! !> @brief
    ! !> Method to increment the number of active Tracers on this list
    ! !---------------------------------------------------------------------------
    ! subroutine lowerNumActive(this)
    ! class(tracerList_class), intent(inout) :: this
    ! this%numActive = this%numActive - 1
    ! end subroutine lowerNumActive
    
    ! !---------------------------------------------------------------------------
    ! !> @author Ricardo Birjukovs Canelas - MARETEC
    ! !> @brief
    ! !> Method that adds a tracer to the list
    ! !---------------------------------------------------------------------------
    ! subroutine addTrc(this, value)
    ! class(tracerList_class), intent(inout) :: this
    ! class(*), intent(in) :: value
    ! this%numActive = this%numActive + 1
    ! call this%add(value)
    ! end subroutine addTrc

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the links of the list
    !---------------------------------------------------------------------------
    subroutine print_tracerList(this)
    class(tracerList_class), intent(in) :: this
    class(*), pointer :: curr
    type(string) :: outext
    outext = this%getSize()
    outext = 'Tracer list has '// outext// ' Tracers'
    call Log%put(outext)
    call this%reset()               ! reset list iterator
    do while(this%moreValues())     ! loop while there are values to print
        call this%printCurrent()
        call this%next()            ! increment the list iterator
    end do
    call this%reset()               ! reset list iterator
    end subroutine print_tracerList

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the current link of the list
    !---------------------------------------------------------------------------
    subroutine print_tracerListCurrent(this)
    class(tracerList_class), intent(in) :: this
    class(*), pointer :: curr
    type(string) :: outext
    curr => this%currentValue() ! get current value
    select type(curr)
    class is (tracer_class)
        call curr%print()
    class default
        outext = '[tracerList_class::print] Unexepected type of content, not a Tracer'
        call Log%put(outext)
        stop
    end select
    end subroutine print_tracerListCurrent


    end module tracerList_mod
