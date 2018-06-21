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


    module tracer_array_mod

    use abstract_container_array_mod
    use tracers_mod

    private
    public :: TracerArray

    type, extends(container_array) :: TracerArray
        integer :: usedLength
    contains
    procedure :: printArray => print_TracerArray
    procedure :: printElement => print_TracerArray_Element
    end type TracerArray

    contains
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the elements of the array
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine print_TracerArray(this)
    class(TracerArray), intent(in) :: this
    class(*), pointer :: curr
    integer :: i
    do i=1, this%usedLength
        curr => this%get(i)
        select type(curr)
        type is (tracer_class)
            !call curr%print()
        class is (paper_class)
            !call curr%print()
        class is (plastic_class)
            !call curr%print()
            class default
            stop '[print_TracerArray]: unexepected type of content: not a shape or derived type'
        end select
    end do
    end subroutine print_TracerArray

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints an element of the array
    !> @param[in] this, index
    !---------------------------------------------------------------------------
    subroutine print_TracerArray_Element(this,index)
    class(TracerArray), intent(in) :: this
    integer, intent(in) :: index
    class(*), pointer :: curr
    if (index .le. this%usedLength) then
        curr => this%get(index)
        select type(curr)
        type is (tracer_class)
            !call curr%print()
        class is (paper_class)
            !call curr%print()
        class is (plastic_class)
            !call curr%print()
            class default
            stop '[print_TracerArray_Element]: unexepected type of content, not a shape or derived type'
        end select
    else
        stop '[print_TracerArray_Element]: index out of bounds'
    endif
    end subroutine print_TracerArray_Element

    end module tracer_array_mod
