!
! To compile:
!
!   Pre-processor now required:
!   /opt/local/bin/gfortran-mp-4.8 -cpp kinds.f90 cla.f90 cla_test.f90 -o cla_test
!
! To run:
!   ./cla_test -h                             : -h is a command line flag
!   ./cla_test --jello 999                    : --jello is long form key, matched with value 999
!   ./cla_test --jello=999                    : --jello is long form key, matched with value 999, alternate syntax
!   ./cla_test -j 999                         : -j is short form key, matched with value 999
!
! Copyright (c) 2015 Edward D. Zaron, Portland, Oregon, USA
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or
! substantial portions of the Software.
!
! Except as contained in this notice, the name(s) of the above copyright holders shall not
! be used in advertising or otherwise to promote the sale, use or other dealings in this
! Software without prior written authorization.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
! OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
!

  program cla_test

    use kinds
    use cla

    implicit none
    character(len=STRLEN)  :: key
    character(len=XSTRLEN) :: description
    integer(kind=int_kind) :: kkind
    character(len=STRLEN)  :: default

    integer(kind=int_kind) :: i
    integer(kind=8)        :: j, jello
    real(kind=real_kind)   :: x
    real(kind=4)           :: y4
    real(kind=8)           :: y8, arg4
    character(len=STRLEN)  :: a
    !character(len=XSTRLEN) :: b
    logical                :: l, ll, arg1, f, s
    character(len=STRLEN)  :: arg2, arg3

    ! Test using command line arguments:
    call cla_init
    
!    ! Test using values supplied in a string:
!    call cla_init(trim('f  edztest_pos2 edztest_pos3 9.626   -j 999 --integer   9 -f'))

    ! Register your keys for key/value pairs.
    ! EACH KEY MUST BEGIN WITH -
    ! (long form, maybe more readable)
    key         = '-i'
    description = 'an integer parameter'
    kkind       = cla_int
    default     = '0'
    call cla_register(key,'--integer',description,kkind,default)
    ! (compact form)
    call cla_register('-j','--jello', 'another int', cla_int    ,'72')
    call cla_register('-x','--xcrypt','real', cla_float  ,'3.14159')
    call cla_register('-y','-yellow', 'a real for testing r4, r8',   cla_float  ,'3.14159')
    call cla_register('-a','--author','str', cla_char, 'Ed Zaron')
    call cla_register('-f','--theflag','a flag', cla_flag,'f')
    call cla_register('-s','--second', 'another flag', cla_flag,'f')
    call cla_register('-l','--logical',  'bool', cla_logical,'t')

    ! Register positional parameters in the order they appear on the
    ! command line.
    call cla_posarg_register('arg1', 'first positional',  cla_logical,'t')    
    call cla_posarg_register('arg2', 'second positional',  cla_char,'Eli S')  
    call cla_posarg_register('arg3', 'third positional',  cla_char,'posthree') 
    call cla_posarg_register('arg4', 'fourth positional',  cla_float,'3.14159265358') 
    
    ! These can be called in any order:
    print *,'*********************************************************'
    print *,'Calling cla_help:'
    call cla_help("cla_test")
    print *,'-------> Returned from cla_help.'
    print *,'*********************************************************'

    print *,'*********************************************************'
    print *,'Calling cla_show:'
    call cla_show
    print *,'-------> Returned from cla_show.'
    print *,'*********************************************************'

    print *,'*********************************************************'
    print *,'Calling cla_validate:'
    call cla_validate("cla_test")
    print *,'-------> Returned from cla_validate.'
    print *,'*********************************************************'

    ! Test to see if the --jello key is present using short and long key:
    ll = cla_key_present('--jello')
    if (ll) then
       print *,' --jello is present'
    else
       print *,' --jello is not present'
    endif

    ll = cla_key_present('-j')
    if (ll) then
       print *,' -j is present'
    else
       print *,' -j is not present'
    endif    
    
    ! Test to see if the -f flag is present:
    ll = cla_key_present('-f')
    if (ll) then
       print *,' -f is present'
    else
       print *,' -f is not present'
    endif

    ll = cla_key_present('-s')
    if (ll) then
       print *,' -s is present'
    else
       print *,' -s is not present'
    endif

    print*,"***********************"
    ! Get values:
    call cla_get('-i',i)
    call cla_get('-j',j)
    call cla_get('--jello',jello)
    call cla_get('-x',x)
    call cla_get('-y',y4)
    call cla_get('-y',y8)
    call cla_get('-a',a)
!    call cla_get_xchar('-b',b) ! Looks like extended character strings not implemented.
    call cla_get('-l',l)
    call cla_get('-f',f)
    call cla_get('-s',s)
    print*,"***********************"
    ! Code crashes if no positional arguments are present --> fixed
    print*,"Get positional arguments"
    call cla_get('arg1',arg1)
    call cla_get('arg2',arg2)
    call cla_get('arg3',arg3)
    call cla_get('arg4',arg4)
    print*,"***********************"
    
    print *,'After processing command line arguments:'
    print *,' i     = ',i
    print *,' j     = ',j,' should equal jello:'
    print *,' jello = ',jello
    print *,' x     = ',x
    print *,' y (r4)= ',y4
    print *,' y (r8)= ',y8
    print *,' a     = ',trim(a)
    print *,' l     = ',l
    print *,' f     = ',f
    print *,' s     = ',s    
    print *,' arg1  = ',arg1
    print *,' arg2  = ',arg2
    print *,' arg3  = ',arg3
    print *,' arg4  = ',arg4    
  end program cla_test

