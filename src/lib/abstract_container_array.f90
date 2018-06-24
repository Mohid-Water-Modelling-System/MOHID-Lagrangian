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
    !> Module that defines an unlimited polymorphic container class and related
    !> methods. A container is a fundamental entity allowing to build data
    !> structures such as lists and arrays.
    !> This is an abstract type, so a derived type must be defined for any
    !> specific contents that may be required. Those derived types should provide
    !> type-specific methods that require type-guards, such as printing.
    !------------------------------------------------------------------------------

    module abstract_container_array_mod

    use container_mod

    private
    !Public access vars
    public :: container_array

    type, abstract :: container_array !< Abstract unlimited polymorphic array class
        private
        class(container), allocatable, dimension(:) :: contents !< Allocatable unlimited polymorphic container array
        integer :: length !< Lenght of the array, for easy access
    contains
    procedure :: resize => resizeArray  !< Grows (adds empty space) or shrinks (discards the last entries) of the array
    procedure :: init => initArray !< Allocates the container array. Deallocates if already allocated
    procedure, non_overridable :: get => getValue !< returns the requested entry (pointer)
    procedure, non_overridable :: put => putValue !< stores a value on the requested index
    procedure, non_overridable :: getLength !< returns the length of the array
    procedure, non_overridable :: getMemSize !< returns the memory occupied by the contents, in bytes (1E-6 mb)
    end type container_array

    logical, target :: MC = .true.

    contains


    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the requested entry (pointer)
    !> @param[in] this, index
    !---------------------------------------------------------------------------
    function getValue(this, index)
    class(container_array), intent(in) :: this
    integer, intent(in) :: index
    class(*), pointer :: getValue
    if (index .le. this%getLength()) then
        getValue => this%contents(index)%getContent()
    else
        stop '[getValue]: index out of bounds'
    endif
    end function getValue

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that stores a value on the requested index
    !> @param[in] this, index, value
    !---------------------------------------------------------------------------
    subroutine putValue(this, index, value)
    class(container_array), intent(inout) :: this
    integer, intent(in) :: index
    class(*), intent(in) :: value
    if (index .le. this%getLength()) then
        call this%contents(index)%storeContent(value)
    else
        stop '[putValue]: index out of bounds'
    endif
    end subroutine putValue

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the length of the array
    !> @param[in] this
    !---------------------------------------------------------------------------
    function getLength(this)
    class(container_array), intent(in) :: this
    integer :: getLength
    getLength = this%length
    end function getLength

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that grows (adds empty space) or shrinks (discards the
    !> last entries) of the array. Use sparsely as this might get expensive
    !> for large array operations. Should think of a way to use move_alloc()
    !> @param[in] this, newsize
    !---------------------------------------------------------------------------
    subroutine resizeArray(this, newsize, initvalue)
    class(container_array), intent(inout) :: this
    integer, intent(in) :: newsize
    class(*), target, optional, intent(in) :: initvalue
    class(*), pointer :: value
    integer :: i, tocopy
    type(container), allocatable, dimension(:) :: temp
    value => MC
    tocopy=min(this%getLength(),newsize)
    allocate(temp(newsize))
    do i=1, tocopy
        call temp(i)%storeContent(this%get(i))
    enddo
    if (present(initvalue)) then
        value => initvalue !pointing to the initial value to store in every entry
    end if
    do i= tocopy+1, newsize
        call temp(i)%storeContent(value)
    end do
    call this%init(newsize,source=temp)
    end subroutine resizeArray

    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that allocates the container array. Deallocates if already allocated
    !> @param[in] this, entries, tocopy
    !---------------------------------------------------------------------------
    subroutine initArray(this, entries, source, initvalue)
    class(container_array), intent(inout) :: this
    integer, intent(in) :: entries
    type(container), dimension(:), optional, intent(in) :: source
    class(*), target, optional, intent(in) :: initvalue
    class(*), pointer :: value
    integer :: i
    value => MC
    if (allocated(this%contents)) then
        deallocate(this%contents)
    end if
    if (.not.present(source)) then
        allocate(this%contents(entries))
        this%length=entries
        if (present(initvalue)) then
            value => initvalue !pointing to the initial value to store in every entry
        end if
        do i=1, size(this%contents)
            call this%put(i,value)
        end do
    else if (present(source)) then !using sourced allocation
        allocate(this%contents, source=source)
        this%length=size(source)
    endif
    end subroutine initArray
    
    !---------------------------------------------------------------------------
    !> @author Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns the memory used by the array in bytes
    !> @param[in] this
    !---------------------------------------------------------------------------
    function getMemSize(this)
    class(container_array), intent(in) :: this
    integer :: getMemSize
    getMemSize = sizeof(this%contents)
    end function getMemSize

    end module abstract_container_array_mod
