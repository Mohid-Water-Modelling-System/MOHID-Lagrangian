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
    use common_modules
    
    implicit none
    private

    type, extends(container_array) :: tracerarray_class
        integer :: numActive = 0  !> number of active Tracers in the array
        integer :: lastActive = 0 !> position of the last active Tracer on the array
    contains
    procedure :: printArray => print_TracerArray
    procedure :: printElement => print_TracerArray_Element
    procedure :: findLastActive
    end type tracerarray_class

    public :: tracerarray_class

    contains

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that traverses the array from the end, finding the last 
    !> active Tracer. Returns its index.
    !---------------------------------------------------------------------------
    integer function findLastActive(this)
    implicit none
    class(tracerarray_class), intent(in) :: this
    class(*), pointer :: curr
    integer :: i
    type(string) :: outext
    i=this%getLength()
    do while (i>0)
        curr => this%get(i)
        select type(curr)
        class is (tracer_class)
            if (curr%now%active) then
                findLastActive = i
                i = 0
            end if
            class default
            outext = '[tracerarray_class::findLastActive]: unexepected type of content not a Tracer'
            call Log%put(outext)
            stop
        end select
        i = i-1
    end do   
    end function findLastActive


    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that prints all the elements of the array
    !> @param[in] this
    !---------------------------------------------------------------------------
    subroutine print_TracerArray(this)
    implicit none
    class(tracerarray_class), intent(in) :: this
    class(*), pointer :: curr
    integer :: i
    do i=1, this%lastActive
        curr => this%get(i)
        select type(curr)
        class is (tracer_class)
            call curr%print()
            class default
            stop '[print_TracerArray]: unexepected type of content not a Tracer'
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
    implicit none
    class(tracerarray_class), intent(in) :: this
    integer, intent(in) :: index
    class(*), pointer :: curr
    if (index .le. this%lastActive) then
        curr => this%get(index)
        select type(curr)
        class is (tracer_class)
            call curr%print()
            class default
            stop '[print_TracerArray_Element]: unexepected type of content not a Tracer'
        end select
    else
        stop '[print_TracerArray_Element]: index out of bounds'
    endif
    end subroutine print_TracerArray_Element

    end module tracer_array_mod
