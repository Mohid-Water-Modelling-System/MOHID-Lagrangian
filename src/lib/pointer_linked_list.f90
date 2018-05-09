    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : pointer_linked_list
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a pointer based linked list class, that encapsulates a
    !> link class and provides a compreensive API. This list is not a container,
    !> just a pointer list.
    !------------------------------------------------------------------------------
    module pointer_linked_list

      use pointer_link_list

      implicit none
      private

      type pointer_list
         private
         class(pointer_link),pointer :: firstLink => null() !> first link in list
         class(pointer_link),pointer :: lastLink => null()  !> last link in list
       contains
         procedure :: add    !> add class(*) to linked list
         procedure :: firstValue  !> return value associated with firstLink
         procedure :: isEmpty     !> return true if list is empty
      end type pointer_list

      !public access types
      public :: pointer_list

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

    end module pointer_linked_list


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    ! Routine Author Name and Affiliation.
    !
    !> @brief
    !> Module that defines a pointer based link. Not a container, just a pointer
    !---------------------------------------------------------------------------
    module pointer_link_list

      implicit none
      private
      
      type pointer_link
         private
         class(*), pointer :: value => null() !> value stored in link
         type(pointer_link), pointer :: next => null()!> next link in list
         contains
         procedure :: getValue    !> return value pointer
         procedure :: nextLink    !> return next pointer
         procedure :: setNextLink !> set next pointer
      end type pointer_link

      interface pointer_link
        procedure constructor !> construct/initialize a link
      end interface

      !public access types
      public :: pointer_link

    contains

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Returns pointer to next link.
      !---------------------------------------------------------------------------
      function nextLink(this)
      class(pointer_link) :: this
      class(pointer_link), pointer :: nextLink
        nextLink => this%next
      end function nextLink

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Sets the next pointer link.
      !---------------------------------------------------------------------------
      subroutine setNextLink(this,next)
      class(pointer_link) :: this
      class(pointer_link), pointer :: next
         this%next => next
      end subroutine setNextLink

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Returns the value the link points to.
      !---------------------------------------------------------------------------
      function getValue(this)
      class(pointer_link) :: this
      class(*), pointer :: getValue
      getValue => this%value
      end function getValue

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Pointer link constructor. Used as pointer_link(value_to_point, next_link)
      !---------------------------------------------------------------------------
      function constructor(value, next)
        class(pointer_link), pointer :: constructor
        class(*) :: value
        class(pointer_link), pointer :: next
        allocate(constructor)
        constructor%next => next
        allocate(constructor%value, source=value)
      end function constructor

    end module pointer_link_list
