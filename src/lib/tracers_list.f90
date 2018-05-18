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

!! integerList.f90

module integer_list_mod
  use abstract_list_mod
  private
  public :: integerList
  type, extends(list) :: integerList
   contains
     !procedure :: addInteger                    ! add integer in list
     procedure :: current => currentInteger ! get integer pointed by iterator
     procedure :: printList => printIntegerList ! print the integer list
     !generic :: add => addInteger
  end type integerList

contains

  subroutine printIntegerList(this)
    class(integerList) :: this
    class(*), pointer :: curr

    call this%reset()
    do while(this%moreValues())
       curr => this%CurrentValue()
       select type(curr)
       type is (integer)
          print *, curr
       end select
       call this%next()
    end do
    call this%reset()
  end subroutine printIntegerList

  !subroutine addInteger(this, value)
  !  class(integerList) :: this
  !  integer value
  !  class(*), allocatable :: v
  !
  !  allocate(v,source=value)
  !  call this%list%add(v)
  !
  !end subroutine addInteger

  integer function currentInteger(this)
    class(integerList) :: this
    class(*), pointer :: v

    v => this%currentValue()
    select type(v)
    type is (integer)
      currentInteger = v
    end select
  end function currentInteger

end module integer_list_mod