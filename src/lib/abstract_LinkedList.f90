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
    !> Module that defines an unlimited polymorphic container list class and related
    !> methods. A container is a fundamental entity allowing to build data
    !> structures such as lists and arrays.
    !> This is an abstract type, so a derived type must be defined for any
    !> specific contents that may be required. Those derived types should provide
    !> type-specific methods that require type-guards, such as printing.
    !------------------------------------------------------------------------------

    module abstract_LinkedList_mod

    use link_mod

    private
    !Public access vars
    public :: linkedlist

    type, abstract :: linkedlist !< Abstract Linked List based on an unlimited polymorphic link class
        private
        class(link), pointer :: firstLink => null()   !< First link in List
        class(link), pointer :: lastLink => null()    !< Last link in List
        class(link), pointer :: currLink => null()    !> list iterator
        integer :: numLinks = 0
    contains
    procedure, non_overridable :: addValue    !< stores a value on the list
    procedure, non_overridable :: getValue    !< get nth value in list
    procedure, non_overridable :: removeCurrent !< Method that removes the current link from a list
    procedure, non_overridable :: remove      !< Method that removes the nth link from a list
    procedure, non_overridable :: getFirst    !< returns the fist link of the list
    procedure, non_overridable :: getLast     !< returns the last link of the list
    procedure, non_overridable :: getSize     !< returns the size of the list
    procedure, non_overridable :: reset       !< reset list iterator
    procedure, non_overridable :: next        !< iterate to next value in list
    procedure, non_overridable :: previous    !< iterate to previous value in list
    procedure, non_overridable :: currentValue!< get current value in list
    procedure, non_overridable :: moreValues  !< more values to iterate?
    generic :: add => addValue
    end type linkedlist

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that stores a value on a new link
    !> @param[this, value, key]
    !---------------------------------------------------------------------------
    subroutine addValue(this, value, key)
    class(linkedlist) :: this
    class(*), intent(in) :: value
    integer, intent(in), optional :: key
    class(link), pointer :: newLink
    if (.not. associated(this%firstLink)) then
        if (present(key)) then
            this%firstLink => link(value, this%firstLink, this%firstLink, key)
        else
            this%firstLink => link(value, this%firstLink, this%firstLink)
        end if
        this%lastLink => this%firstLink
    else
        if (present(key)) then
            newLink => link(value, this%lastLink, this%lastLink%nextLink(), key)
        else
            newLink => link(value, this%lastLink, this%lastLink%nextLink())
        end if
        call this%lastLink%setNextLink(newLink)
        this%lastLink => newLink
    end if
    this%numLinks = this%numLinks + 1
    end subroutine addValue
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that removes a link from the list
    !---------------------------------------------------------------------------    
    subroutine removeCurrent(this)
    class(linkedlist), intent(inout) :: this    
    class(link), pointer :: previouslink
    class(link), pointer :: nextlink 
    
    previouslink => this%currLink%previousLink()
    nextlink => this%currLink%nextLink()
    
    if (associated(this%currLink,this%firstLink)) then !This is the first link
        this%firstLink => nextlink
    end if
    if (associated(previouslink)) then
        call previouslink%setNextLink(nextlink)
    end if
    if (associated(nextlink)) then
        call nextlink%setPreviousLink(previouslink)
    end if
    
    call this%currLink%removeLink()
    deallocate(this%currLink)    
    this%currLink => nextlink
    this%numLinks = this%numLinks - 1
    
    end subroutine removeCurrent
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that removes the nth link from a list
    !---------------------------------------------------------------------------    
    subroutine remove(this, n)
    class(linkedlist), intent(inout) :: this
    integer, intent(in) :: n
    class(link), pointer :: previouslink
    class(link), pointer :: nextlink
    integer :: i
    if (associated(this%firstLink)) then
        if (this%numLinks>=n) then
            call this%reset()
            do i=1, n-1    !iterating trough the list until the desired position
                call this%next()
            end do
            if (this%moreValues()) then
                call this%removeCurrent()
            end if
        else
            stop '[LinkedList::remove]: index out of bounds'
        end if
    end if
    end subroutine remove

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the first link of the list
    !---------------------------------------------------------------------------
    function getFirst(this) result(firstlink)
    class(linkedlist) :: this
    class(link), pointer :: firstlink
    firstlink => this%firstLink
    end function getFirst

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the last link of the list
    !---------------------------------------------------------------------------
    function getLast(this) result(lastLink)
    class(linkedlist) :: this
    class(link), pointer :: lastLink
    lastLink => this%lastLink
    end function getLast

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the size (number of links) of a list
    !---------------------------------------------------------------------------
    pure integer function getSize(this)
    class(linkedlist), intent(in) :: this
    getSize = this%numLinks
    end function getSize

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the value of the nth link of a list
    !---------------------------------------------------------------------------
    function getValue(this,n) result(res)
    class(linkedlist), intent(in) :: this
    integer, intent(in) :: n
    class(*), pointer :: res
    integer :: i
    type(link),pointer :: alink
    if (associated(this%firstLink)) then
        if (this%numLinks>=n) then
            call this%reset()
            do i=1, n-1    !iterating trough the list until the desired position            
                call this%next()
            end do
            if (this%moreValues()) then
            res => this%currLink%get()
            else
                stop '[LinkedList::getValue]: link non-existant, something went wrong!'
            end if
            call this%reset()
        else
            stop '[LinkedList::getValue]: index out of bounds'
        end if
    else
        stop '[LinkedList::getValue]: list is empty'
    end if
    end function getValue

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the value of the current link
    !---------------------------------------------------------------------------
    function currentValue(this)
    class(linkedlist) :: this
    class(*), pointer :: currentValue
    currentValue => this%currLink%get()
    end function CurrentValue

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the next link in the list
    !---------------------------------------------------------------------------
    subroutine next(this)
    class(linkedlist) :: this
    this%currLink => this%currLink%nextLink()
    end subroutine next
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the previous link in the list
    !---------------------------------------------------------------------------
    subroutine previous(this)
    class(linkedlist) :: this
    this%currLink => this%currLink%previousLink()
    end subroutine previous

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a logical with signaling if the current link is ok
    !---------------------------------------------------------------------------
    pure logical function moreValues(this)
    class(linkedlist), intent(in) :: this
    moreValues = associated(this%currLink)
    end function moreValues

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that resets the list iterator
    !---------------------------------------------------------------------------
    subroutine reset(this)
    class(linkedlist) :: this
    this%currLink => this%firstLink
    end subroutine reset

    end module abstract_LinkedList_mod
