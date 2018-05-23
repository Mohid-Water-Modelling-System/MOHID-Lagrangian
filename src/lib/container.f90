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
    !> Module that defines an unlimited polymorphic container class and related methods. 
    !> A container is a fundamental entity allowing to build data structures such as lists
    !> and arrays.
    !------------------------------------------------------------------------------

    module container_mod

    implicit none
    private
    !Public access vars
    public :: container

    type container !< Unlimited polymorphic container class 
        private !contents are only accessible trough the methods, no direct access is allowed
        class(*), pointer :: value => null() !< value stored in container
    contains
    procedure :: getContent     !< returns stored content (pointer)
    procedure :: storeContent   !< stores the provided values (sourced allocation)
    procedure :: printContainer !< prints container contents (only primitive types implemented)
    end type container

    interface container
        procedure constructor !< construct/initialize a container
    end interface

    contains
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that returns a pointer to the values stored in the container
    !> @param[this]
    !---------------------------------------------------------------------------
    function getContent(this)
    class(container), intent(in) :: this
    class(*), pointer :: getContent
    getContent => this%value
    end function getContent

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method that stores the provided value in the container using 
    !> sourced allocation
    !> @param[this, to_store]
    !---------------------------------------------------------------------------
    subroutine storeContent(this,to_store)
    class(container), intent(inout) :: this
    class(*), intent(in) :: to_store
    allocate(this%value, source=to_store)
    end subroutine storeContent
    
    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Method to print the stored value. Only knows about instrinsic types, 
    !> ignores (but warns) if other types are passed.
    !> @param[this]
    !---------------------------------------------------------------------------
    subroutine printContainer(this)
    class(container), intent(in) :: this
    select type(v => this%value)
    type is (integer)
        print *, v
    type is (character(*))
        print *, v(1:1)
    type is (real)
        print *, v
        class default
        print*, "[printContainer]: don't know how to print this value, ignoring"
    end select
    end subroutine printContainer

    !---------------------------------------------------------------------------
    !> @Ricardo Birjukovs Canelas - MARETEC
    !> @brief
    !> Container constructor, can be used with the 'container' name since it is 
    !> defined as an interface
    !> @param[to_store]
    !---------------------------------------------------------------------------
    function constructor(to_store)
    class(container), pointer :: constructor
    class(*), intent(in) :: to_store
    allocate(constructor)
    allocate(constructor%value, source=to_store)
    end function constructor

    end module container_mod
