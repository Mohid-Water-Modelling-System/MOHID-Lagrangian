    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : sources_list
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : August 2018
    ! REVISION      : Canelas 0.1
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module to hold the Sources linked list class and its methods. This class
    !> defines a double linked list to store any variable type, but with specific
    !> methods with type guards for Source objects. The class allows for insertion,
    !> deletion and iteration of the desired contents
    !------------------------------------------------------------------------------

    module sources_list_mod

    use sources_mod
    use abstract_LinkedList_mod
    use common_modules

    implicit none
    private

    type, extends(linkedlist) :: sourceList_class !< List of Source class
    contains
    procedure :: print => print_sourceList
    procedure :: printCurrent => print_sourceListCurrent
    end type sourceList_class

    public :: sourceList_class

    contains

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the links of the list
    !---------------------------------------------------------------------------
    subroutine print_sourceList(this)
    class(sourceList_class), intent(in) :: this
    class(*), pointer :: curr
    call this%reset()               ! reset list iterator
    do while(this%moreValues())     ! loop while there are values to print
        call this%printCurrent()
        call this%next()            ! increment the list iterator
    end do
    call this%reset()               ! reset list iterator
    end subroutine print_sourceList

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints the current link of the list
    !---------------------------------------------------------------------------
    subroutine print_sourceListCurrent(this)
    class(sourceList_class), intent(in) :: this
    class(*), pointer :: curr
    type(string) :: outext
    curr => this%currentValue() ! get current value
    select type(curr)
    class is (source_class)
        call curr%print()
        class default
        outext = '[sourceList_class::print] Unexepected type of content, not a Source'
        call Log%put(outext)
        stop
    end select
    end subroutine print_sourceListCurrent


    end module sources_list_mod
