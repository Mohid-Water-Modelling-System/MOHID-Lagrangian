    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : container_linked_list
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a container linked list class, that encapsulates a
    !> link class and provides a compreensive API. This list is a container,
    !> it should work as a polymorphic array of objects.
    !------------------------------------------------------------------------------
    module container_linked_list

      use container_link_class

      private
      public :: container_list

      type container_list
         private
         class(pointer_link),pointer :: firstLink => null() !> first link in list
         class(pointer_link),pointer :: lastLink => null()  !> last link in list
       contains
         procedure :: add    !> add class(*) to linked list
         procedure :: firstValue  !> return value associated with firstLink
         procedure :: isEmpty     !> return true if list is empty
      end type container_list

    contains

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> This is the generic entry point for our list. It adds a link that uses a
      !> a polymorphic pointer so it can be used for any target
      !
      !> @param[value]
      !---------------------------------------------------------------------------
      subroutine add(this, value)
        class(pointer_list) :: this
        class(*) :: value
        class(pointer_link), pointer :: newLink

        if (.not. associated(this%firstLink)) then
           this%firstLink => pointer_link(value, this%firstLink)
           this%lastLink => this%firstLink
        else
           newLink => pointer_link(value, this%lastLink%nextLink())
           call this%lastLink%setNextLink(newLink)
           this%lastLink => newLink
        end if

      end subroutine add

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Returns the first link of the list.
      !---------------------------------------------------------------------------
      function firstValue(this)
        class(pointer_list) :: this
        class(*), pointer :: firstValue

        firstValue => this%firstLink%getValue()

      end function firstValue

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Returns a logical with the list condition (empty(true) or filled(false)).
      !---------------------------------------------------------------------------
      function isEmpty(this)
        class(pointer_list) :: this
        logical isEmpty

        if (associated(this%firstLink)) then
           isEmpty = .false.
        else
           isEmpty = .true.
        endif
      end function isEmpty

    end module container_linked_list
