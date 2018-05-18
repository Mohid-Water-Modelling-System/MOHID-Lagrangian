    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Lagrangian Tracer
    ! MODULE        : container_link_class
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a container link class, that will be encapsulataded
    !> by a linked list class. This list is a container,
    !> it should work as a polymorphic array of objects.
    !------------------------------------------------------------------------------
    module container_link_class
      private
      public :: container_link
      
      type container_link
         private
         class(*), allocatable :: value !> value stored in link
         type(container_link), pointer :: next => null()!> next link in list
         contains
         procedure :: getValue    !> return value
         procedure :: nextLink    !> return next pointer
         procedure :: setNextLink !> set next link
      end type container_link

      interface container_link
        procedure constructor !> construct/initialize a link
      end interface

    contains

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Returns next link.
      !---------------------------------------------------------------------------
      function nextLink(this)
      class(container_link) :: this
      type(container_link) :: nextLink
        nextLink = this%next
      end function nextLink

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Sets the next link.
      !---------------------------------------------------------------------------
       subroutine setNextLink(this,next)
       class(container_link) :: this
       class(container_link), pointer :: next
          this%next => next
       end subroutine setNextLink

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Returns the value the link contains.
      !---------------------------------------------------------------------------
       function getValue(this)
       class(container_link) :: this
       class(*), allocatable :: getValue
         allocate(getValue, source=this%value)
       end function getValue

      !---------------------------------------------------------------------------
      !> @Ricardo Birjukovs Canelas - MARETEC
      ! Routine Author Name and Affiliation.
      !
      !> @brief
      !> Container link constructor. Used as container_link(value_to_store, next_link)
      !---------------------------------------------------------------------------
      function constructor(value, next)
        type(container_link) :: constructor
        class(*), intent(in) :: value
        class(container_link), pointer :: next
        constructor%next => next
        allocate(constructor%value, source=value)
      end function constructor

    end module container_link_class
