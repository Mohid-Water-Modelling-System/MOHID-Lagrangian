    !
    !     Copyright (c) 2015, NVIDIA CORPORATION.  All rights reserved.
    !
    ! NVIDIA CORPORATION and its licensors retain all intellectual property
    ! and proprietary rights in and to this software, related documentation
    ! and any modifications thereto.  Any use, reproduction, disclosure or
    ! distribution of this software and related documentation without an express
    ! license agreement from NVIDIA CORPORATION is strictly prohibited.
    !

    !          THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT
    !   WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT
    !   NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR
    !   FITNESS FOR A PARTICULAR PURPOSE.
    !

    !! link.f90

    module link_mod
    private
    public :: link
    type link
        private
        class(*), pointer :: value => null() ! value stored in link
        type(link), pointer :: next => null()! next link in list
    contains
    procedure :: getValue    ! return value pointer
    procedure :: printLinks  ! print linked list starting with this link
    procedure :: nextLink    ! return next pointer
    procedure :: setNextLink ! set next pointer
    end type link

    interface link
    procedure constructor ! construct/initialize a link
    end interface

    contains

    function nextLink(this)
    class(link) :: this
    class(link), pointer :: nextLink
    
    nextLink => this%next
    
    end function nextLink

    subroutine setNextLink(this,next)
    class(link) :: this
    class(link), pointer :: next
    
    this%next => next
    
    end subroutine setNextLink

    function getValue(this)
    class(link) :: this
    class(*), pointer :: getValue
    
    getValue => this%value
    
    end function getValue

    subroutine printLink(this)
    class(link) :: this

    select type(v => this%value)
    type is (integer)
        print *, v
    type is (character(*))
        print *, v(1:1)
    type is (real)
        print *, v
        class default
        stop 'printLink: unexepected type for link'
    end select

    end subroutine printLink

    subroutine printLinks(this)
    class(link) :: this
    class(link), pointer :: curr

    call printLink(this)
    curr => this%next
    do while(associated(curr))
        call printLink(curr)
        curr => curr%next
    end do

    end subroutine

    function constructor(value, next)
    class(link),pointer :: constructor
    class(*) :: value
    class(link), pointer :: next
    
    allocate(constructor)
    constructor%next => next
    allocate(constructor%value, source=value)
    
    end function constructor

    end module link_mod
