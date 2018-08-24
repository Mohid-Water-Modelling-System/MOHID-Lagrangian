    !MARETEC - Research Centre for Marine, Environment and Technology
    !Copyright (C) 2018  Ricardo Birjukovs Canelas
    !
    !This program is free software: you can redistribute it and/or modify
    !it under the terms of the GNU General Public License as published by
    !the Free Software Foundation, either version 3 of the License, or
    !(at your option) any later version.
    !
    !This program is distributed in the hope that it will be useful,
    !but WITHOUT ANY WARRANTY; without even the implied warranty of
    !MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !GNU General Public License for more details.
    !
    !You should have received a copy of the GNU General Public License
    !along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    !------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
    !------------------------------------------------------------------------------
    !
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : May 2018
    ! REVISION      : Canelas 0.2
    !> @author
    !> Ricardo Birjukovs Canelas
    !
    ! DESCRIPTION:
    !> Module that defines a link based on an unlimited polymorphic container class
    !------------------------------------------------------------------------------

    module link_mod

    use container_mod

    implicit none
    private
    !Public access vars
    public :: link

    type, extends(container) :: link !< link class based on a container 
        private !contents are only accessible trough the methods, no direct access is allowed
        integer, public :: key
        type(link), pointer :: next => null() !< pointer to a next link
        type(link), pointer :: previous => null() !< pointer to a previous link
    contains
    procedure :: get => getValue        !< returns stored content
    procedure :: nextLink               !< gets the next link
    procedure :: previousLink           !< gets the previous link
    procedure :: setNextLink            !< sets the next link pointer
    procedure :: setPreviousLink        !< sets the previous link pointer
    procedure :: removeLink
    end type link

    interface link
        procedure constructor !< construct/initialize a container
    end interface

    contains
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a pointer to the values stored in the container in
    !> this link
    !---------------------------------------------------------------------------
    function getValue(this)
    class(link) :: this
    class(*), pointer :: getValue
    getValue => this%getContent()
    end function getValue

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a pointer to the next link in a list
    !---------------------------------------------------------------------------
    function nextLink(this)
    class(link) :: this
    class(link), pointer :: nextLink
    nextLink => this%next
    end function nextLink
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a pointer to the previous link in a list
    !---------------------------------------------------------------------------
    function previousLink(this)
    class(link) :: this
    class(link), pointer :: previousLink
    previousLink => this%previous
    end function previousLink
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to set the next link in a list.
    !---------------------------------------------------------------------------
    subroutine setNextLink(this,next)
    class(link) :: this
    class(link), pointer :: next
    this%next => next
    end subroutine setNextLink
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to set the previous link in a list.
    !---------------------------------------------------------------------------
    subroutine setPreviousLink(this,prev)
    class(link) :: this
    class(link), pointer :: prev
    this%previous => prev
    end subroutine setPreviousLink
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to remove a link in a list.
    !---------------------------------------------------------------------------
    subroutine removeLink(this)
    class(link), intent(inout) :: this
    call this%deleteContent()
    end subroutine removeLink

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Link constructor, can be used with the 'link' name since it was defined 
    !> as such in an interface declaration
    !> @param[to_store, prev, next]
    !---------------------------------------------------------------------------
    function constructor(to_store, prev, next, key)
    class(link), pointer :: constructor
    class(*), intent(in) :: to_store
    class(link), pointer, intent(in) :: prev
    class(link), pointer, intent(in) :: next
    integer, intent(in), optional :: key
    allocate(constructor)
    call constructor%setPreviousLink(prev)
    call constructor%setNextLink(next)
    call constructor%storeContent(to_store)
    if (present(key)) then
        constructor%key = key
    end if
    end function constructor

    end module link_mod
