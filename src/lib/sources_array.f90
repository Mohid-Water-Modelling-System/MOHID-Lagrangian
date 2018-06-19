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


    module sources_array_mod

    use abstract_container_array_mod
    use sources_mod

    private
    public :: SourceArray

    type, extends(container_array) :: SourceArray
        integer :: usedLength
    contains
    procedure :: printArray => print_SourceArray
    procedure :: printElement => print_SourceArray_Element
    end type SourceArray

    contains

    subroutine print_SourceArray(this)
    class(SourceArray), intent(in) :: this
    class(*), pointer :: curr
    integer :: i
    do i=1, this%getLength()
        curr => this%get(i)
        select type(curr)
        type is (source_class)
            !call curr%print()
            class default
            stop '[print_SourceArray]: unexepected type of content: not a Source or derived type'
        end select
    end do
    end subroutine print_SourceArray

    subroutine print_SourceArray_Element(this,index)
    class(SourceArray), intent(in) :: this
    integer, intent(in) :: index
    class(*), pointer :: curr
    if (index .le. this%getLength()) then
        curr => this%get(index)
        select type(curr)
        type is (source_class)
            !call curr%print()
            class default
            stop '[print_SourceArray_Element]: unexepected type of content, not a Source or derived type'
        end select
    else
        stop '[print_SourceArray_Element]: index out of bounds'
    endif
    end subroutine print_SourceArray_Element

    end module sources_array_mod
