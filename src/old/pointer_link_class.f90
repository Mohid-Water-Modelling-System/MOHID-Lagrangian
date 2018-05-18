    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : pointer_link_class
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a pointer based link class, to be encapsulated in a
    !> linked list class. This list is not a container,
    !> just a pointer list.
    !------------------------------------------------------------------------------
    module pointer_link_class

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

    end module pointer_link_class
