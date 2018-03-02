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

module cla

  use kinds
#ifdef f2003
  use, intrinsic :: iso_fortran_env, only : &
       input_unit  => stdin, &
       output_unit => stdout, &
       error_unit  => stderr
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif
  
  implicit none

  ! Command Line Arguments
  ! A key-value parser for the commandline

  integer(kind=int_kind), parameter :: &
       cla_int   = 1, &
       cla_float = 2, &
       cla_char  = 3, &
       cla_xchar = 4, & ! NOT IMPLEMENTED
       cla_logical=5, &
       cla_flag  = 6

  ! Optionally, these routines can parse a string contaning arguments, cla_cla,
  ! rather than the command line arguments.
  integer(kind=int_kind), parameter, private :: CLALEN=1024
  character(len=CLALEN), private  :: cla_cla
  integer(kind=int_kind), private :: cla_cla_len
  
  character(len=STRLEN), dimension(6) :: cla_kindstr
  character(len=STRLEN), private :: cla_empty
  character(len=STRLEN), dimension(6) :: cla_true_str
 
  type, private :: cla_t
     character(len=2)  :: key
     character(len=STRLEN)  :: longkey 
     character(len=XSTRLEN) :: description
     integer(kind=int_kind) :: kind
     character(len=STRLEN)  :: default
  end type cla_t

  type, private :: cla_posarg_t
     character(len=STRLEN)  :: key
     character(len=XSTRLEN) :: description
     integer(kind=int_kind) :: kind
     character(len=STRLEN)  :: default     
  end type cla_posarg_t
  
  
  type(cla_t), private, dimension(:), pointer :: cla_registry
  type(cla_posarg_t), private, dimension(:), pointer :: cla_posarg_registry
  
  
  integer(kind=int_kind), private :: cla_num
  integer(kind=int_kind), private :: cla_posarg_num  

  interface cla_init
     module procedure &
          cla_init_default, & ! no input parameters ==> read and parse the command line
          cla_init_str        ! string input parameter ==> read and parse the string instead of the command line
  end interface cla_init
     
  interface cla_get
     module procedure &
          cla_get_float_r4, &
          cla_get_float_r8, &
          cla_get_int_i4, &
          cla_get_int_i8, &
          cla_get_char, &
          cla_get_logical
  end interface

  contains

    integer function cla_command_argument_count()
      ! Define these wrappers for these intrinsic functions because we want to implement
      ! a cla parser that can work with plain strings, read from stdin for example,
      ! not just from the command line.
      ! The purpose is to easily adapt a program that runs using command line inputs into
      ! a program that can run with string input it eats from a pipe.
      implicit none
      CHARACTER(len=CLALEN) :: outs
      INTEGER               :: i, k, n
           
      if (cla_cla_len == 0) then
         cla_command_argument_count = command_argument_count()
      else
         ! Count number of spaces to get
         ! number of parameters minus 1:
         cla_command_argument_count = 0
         ! Now count arguments:
!         print *,'Counting arguments in :',trim(cla_cla)
         cla_cla_len = len_trim(cla_cla)
         if (cla_cla_len == 0) then
            return
         end if
         do i=1,cla_cla_len
            if (cla_cla(i:i) == ' ') then
               cla_command_argument_count = cla_command_argument_count + 1
            end if
         end do
         cla_command_argument_count = cla_command_argument_count + 1
!         print *,'cla_command_argument_count found ',cla_command_argument_count
      end if
    end function cla_command_argument_count

    subroutine cla_get_command_argument(i,arg)
      implicit none
      integer :: i,n,nm1,nn
      character(len=*) :: arg
      if (cla_cla_len == 0) then
         call get_command_argument(i,arg)
      else
         nm1 = 0
         nn = 0
         do n=1,i
            nm1 = nn + nm1
            nn = index(cla_cla(nm1+1:),' ')
         end do
         arg(1:nn) = cla_cla(nm1+1:nm1+nn-1)
         arg(nn:) = ' '
!         print *,'trim(cla_cla)=',trim(cla_cla)
!         print *,'get_command_argument got : arg=',i,' val=',trim(arg)
      end if
    end subroutine cla_get_command_argument

    subroutine cla_message(message)
    character(LEN=*) message
    ! Default stop with message without print or stop statements. 
    ! May need to be modified for, e.g. MPI codes 
    write(stdout,*)message
    end subroutine  
  
    subroutine cla_fatal(message)
    character(LEN=*) message
    ! Default stop with message without print or stop statements. 
    ! May need to be modified for, e.g. MPI codes 
    write(stderr,*)message
    stop 6
    end subroutine

    subroutine cla_read_str(cla_input_str)
      implicit none
      character(len=*) :: cla_input_str
      integer :: i,n,k
      character(CLALEN) :: outs
      cla_cla = trim(adjustl(cla_input_str))
      cla_cla_len = len_trim(cla_cla)
      ! Replace whitespace with single space:
      outs = cla_cla
      n = 0  ; k=cla_cla_len                  ! k=index last non-blank (may be null)
      DO i = 1,k-1                            ! dont process last char yet
         n = n+1 ; outs(n:n) = cla_cla(i:i)
         IF (cla_cla(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
      END DO
      n = n+1  ; outs(n:n)  = cla_cla(k:k)    ! last non-blank char output (may be null)
      IF (n < k) outs(n+1:) = ' '             ! pad trailing blanks
      cla_cla = outs
      cla_cla_len = len_trim(cla_cla)
    end subroutine cla_read_str
    
    subroutine cla_init_str(cla_input_str)
      implicit none
      character(len=*) :: cla_input_str
      integer :: i,n,k
      character(CLALEN) :: outs
      cla_cla = trim(adjustl(cla_input_str))
      cla_cla_len = len_trim(cla_cla)
      ! Replace whitespace with single space:
      outs = cla_cla
      n = 0  ; k=cla_cla_len                  ! k=index last non-blank (may be null)
      DO i = 1,k-1                            ! dont process last char yet
         n = n+1 ; outs(n:n) = cla_cla(i:i)
         IF (cla_cla(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
      END DO
      n = n+1  ; outs(n:n)  = cla_cla(k:k)    ! last non-blank char output (may be null)
      IF (n < k) outs(n+1:) = ' '             ! pad trailing blanks
      cla_cla = outs
      cla_cla_len = len_trim(cla_cla)

      ! Allocate a zero size registry, just so that it gets
      ! associated.
      cla_num = 0
      allocate(cla_registry(0))
      allocate(cla_posarg_registry(0))      
      cla_kindstr(cla_int)     = 'integer'
      cla_kindstr(cla_float)   = 'float'
      cla_kindstr(cla_char)    = 'character'
      cla_kindstr(cla_xchar)   = 'xcharacter' !NOT IMPLEMENTED
      cla_kindstr(cla_logical) = 'logical'
      cla_kindstr(cla_flag)    = 'flag'
      cla_empty='THIS_IS_THE_EMPTY_STRING'
      cla_true_str(1)='true'
      cla_true_str(2)='on'
      cla_true_str(3)='1'
      cla_true_str(4)='t'
      cla_true_str(5)='T'
      cla_true_str(6)='.true.'
    end subroutine cla_init_str
    
    subroutine cla_init_default
      ! Allocate a zero size registry, just so that it gets
      ! associated.
      cla_num = 0
      allocate(cla_registry(0))
      allocate(cla_posarg_registry(0))      
      cla_kindstr(cla_int)     = 'integer'
      cla_kindstr(cla_float)   = 'float'
      cla_kindstr(cla_char)    = 'character'
      cla_kindstr(cla_xchar)   = 'xcharacter' !NOT IMPLEMENTED
      cla_kindstr(cla_logical) = 'logical'
      cla_kindstr(cla_flag)    = 'flag'
      cla_empty='THIS_IS_THE_EMPTY_STRING'
      cla_true_str(1)='true'
      cla_true_str(2)='on'
      cla_true_str(3)='1'
      cla_true_str(4)='t'
      cla_true_str(5)='T'
      cla_true_str(6)='.true.'
      ! Set cla_cla_len = 0 in default case when command line is actually to be parsed
      ! rather than in input string:
      cla_cla_len = 0
    end subroutine cla_init_default

    subroutine cla_posarg_register(key,description,kkind,default)
      character(len=*) :: key
      character(len=*) :: description
      integer(kind=int_kind) :: kkind
      character(len=*) :: default
      type(cla_posarg_t), dimension(:), pointer :: cla_posarg_registry_tmp
      integer(kind=int_kind) :: i

      ! This is a dumb way to increase the size of the
      ! registry of command line arguments, but there
      ! should not be so many arguments that either speed
      ! or memory is an issue.
      allocate(cla_posarg_registry_tmp(cla_posarg_num+1))
      do i=1,cla_posarg_num
         cla_posarg_registry_tmp(i)%key         = cla_posarg_registry(i)%key
         cla_posarg_registry_tmp(i)%description = cla_posarg_registry(i)%description
         cla_posarg_registry_tmp(i)%kind        = cla_posarg_registry(i)%kind
         cla_posarg_registry_tmp(i)%default     = cla_posarg_registry(i)%default
         if (index(trim(key),' ') /= 0) then
            call cla_fatal('Error: cla_posarg key contains a space character.')
         end if
         if (cla_str_eq(trim(cla_posarg_registry(i)%key),trim(key))) then
            call cla_fatal('cla_posarg key already been registered' // &
                           cla_posarg_registry(i)%key)
         end if
      end do
      cla_posarg_num = cla_posarg_num + 1
      deallocate(cla_posarg_registry)
      allocate(cla_posarg_registry(cla_posarg_num))
      do i=1,cla_posarg_num-1
         cla_posarg_registry(i)%key         = cla_posarg_registry_tmp(i)%key
         cla_posarg_registry(i)%description = cla_posarg_registry_tmp(i)%description
         cla_posarg_registry(i)%kind        = cla_posarg_registry_tmp(i)%kind
         cla_posarg_registry(i)%default     = cla_posarg_registry_tmp(i)%default
      end do
      i = cla_posarg_num
      cla_posarg_registry(i)%key         = key
      cla_posarg_registry(i)%description = description
      cla_posarg_registry(i)%description = description
      cla_posarg_registry(i)%kind        = kkind
      cla_posarg_registry(i)%default     = default
      deallocate(cla_posarg_registry_tmp)
    end subroutine
    
    subroutine cla_register(key,longkey,description,kkind,default)
      character(len=2) :: key
      character(len=*) :: longkey
      character(len=*) :: description
      integer(kind=int_kind) :: kkind
      character(len=*) :: default
      type(cla_t), dimension(:), pointer :: cla_registry_tmp
      integer(kind=int_kind) :: i

      if (len(key) > 2)then
        call cla_fatal("The short key should be a dash plus one character (-e)")
      end if
      
      ! This is a dumb way to increase the size of the
      ! registry of command line arguments, but there
      ! should not be so many arguments that either speed
      ! or memory is an issue.
      allocate(cla_registry_tmp(cla_num+1))
      do i=1,cla_num
         cla_registry_tmp(i)%key         = cla_registry(i)%key
         cla_registry_tmp(i)%longkey     = cla_registry(i)%longkey
         cla_registry_tmp(i)%description = cla_registry(i)%description
         cla_registry_tmp(i)%kind        = cla_registry(i)%kind
         cla_registry_tmp(i)%default     = cla_registry(i)%default
         if (index(trim(key),' ') /= 0) then
            call cla_fatal('Attempt to register cla key containing a space.')
         end if
         if (index(trim(longkey),' ') /= 0) then
            call cla_fatal('Attempt to register long key containing space')
         end if         
         if (cla_str_eq(trim(cla_registry(i)%key),trim(key))) then
            call cla_fatal('Attempt to register cla key already registered'// &
                           trim(key))
         end if
      end do
      cla_num = cla_num + 1
      deallocate(cla_registry)
      allocate(cla_registry(cla_num))
      do i=1,cla_num-1
         cla_registry(i)%key         = cla_registry_tmp(i)%key
         cla_registry(i)%longkey     = cla_registry_tmp(i)%longkey
         cla_registry(i)%description = cla_registry_tmp(i)%description
         cla_registry(i)%kind        = cla_registry_tmp(i)%kind
         cla_registry(i)%default     = cla_registry_tmp(i)%default
      end do
      i = cla_num
      cla_registry(i)%key         = key
      cla_registry(i)%longkey     = longkey
      cla_registry(i)%description = description
      cla_registry(i)%kind        = kkind
      cla_registry(i)%default     = default
      deallocate(cla_registry_tmp)
    end subroutine

    subroutine cla_show
      integer(kind=int_kind) :: i
      character(len=STRLEN)  :: value
      character(len=STRLEN)  :: i_str
      call cla_message('General usage:')
      call cla_message('  command -[key] [value] --[longkey] [value] -[flag] [positional arguments]')
      call cla_message('  The key/value pairs must be matched if they appear.')
      call cla_message('  Key/value pairs and flags may be in any order.')
      call cla_message(' ')
      call cla_message('The following command line arguments and switches are expected:')
      do i=1,cla_num
         write(i_str,'(i5)')i
         call cla_message('---------- i: '// trim(i_str))
         call cla_message('         key: '// trim(cla_registry(i)%key))
         call cla_message('     longkey: '// trim(cla_registry(i)%longkey))
         call cla_message(' description: '// trim(cla_registry(i)%description))
         call cla_message('        kind: '// &
                          trim(cla_kindstr(cla_registry(i)%kind)))
         call cla_message('     default: '// trim(cla_registry(i)%default))
      end do
      call cla_message(' ')
      call cla_message('The following positional (non-keyword) arguments are expected:')
      do i=1,cla_posarg_num
         write(i_str,'(i5)')i
         call cla_message('---------- i: '// trim(i_str))
         call cla_message('         key: '// trim(cla_posarg_registry(i)%key))
         call cla_message(' description: '// trim(cla_posarg_registry(i)%description))
         call cla_message('        kind: '// &
                          trim(cla_kindstr(cla_posarg_registry(i)%kind)))
         call cla_message('     default: '// trim(cla_posarg_registry(i)%default))
         if ( cla_key_present(trim(cla_registry(i)%key)) ) then
            call cla_get_char(trim(cla_registry(i)%key),value)
            call cla_message('    present?: T')
            if (cla_registry(i)%kind == cla_flag) then
            else
               call cla_message('       value: '// trim(value))
            endif
         else
            call cla_message('    present?: F')
         endif
      end do
      
      call cla_message(' ')
      call cla_message('Also, -?, -h, -H, -help, --help, and --usage are recognized.')
      call cla_message(' ')
    end subroutine cla_show

    subroutine cla_help(cmd_name)
      character(len=*) :: cmd_name
      integer(kind=int_kind) :: i
      character(len=256) :: cmd_usage
      cmd_usage = ""
      do i=1,cla_num
        if (cla_registry(i)%kind == cla_flag) then
          cmd_usage = trim(cmd_usage) // " [" // trim(cla_registry(i)%key) // "]"
        else
          cmd_usage = trim(cmd_usage) // " [" // trim(cla_registry(i)%key) // "=" // &
          trim(cla_registry(i)%default) // "]"
        endif
      enddo
      do i=1,cla_posarg_num
        cmd_usage = trim(cmd_usage) // " " // trim(cla_posarg_registry(i)%key)
      enddo
      write(stdout,*)'General usage:'
      write(stdout,*)'  ',cmd_name, trim(cmd_usage)
      write(stdout,*)' '
      write(stdout,*)'Options and flags {default values}:'
      if (cla_num == 0) write(stdout,*)"None"
      do i=1,cla_num
         if (cla_registry(i)%kind == cla_flag) then
            write(stdout,'(1x,a,1x,a24,":",4x,a)')trim(cla_registry(i)%key), &
                                       trim(cla_registry(i)%longkey), &
                                       trim(cla_registry(i)%description)
         else
            write(stdout,'(1x,a,1x,a24,":",4x,a,2x,"{",a,"}")')trim(cla_registry(i)%key), &
                                 trim(cla_registry(i)%longkey), &   
                                 trim(cla_registry(i)%description), & 
                                 trim(cla_registry(i)%default) 
         endif
      end do
      write(stdout,*)' '

      write(stdout,*)'Positional arguments:'
      if (cla_posarg_num == 0) write(stdout,*)"None"
      do i=1,cla_posarg_num
        write(stdout,'(1x,a,":",1x,a,4x,a)')trim(cla_posarg_registry(i)%key), &
                                trim(cla_posarg_registry(i)%description)
      end do
      
      write(stdout,*)' '
      write(stdout,*)'Also, -?, -h, -H, -help, --help, and --usage are recognized.'
      write(stdout,*)' '
    end subroutine cla_help

    integer function cla_eq(str1,str2)
      implicit none
      character(*) :: str1, str2
      cla_eq = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
    end function cla_eq

    logical function cla_key_arg_match(key,longkey,arg)
      implicit none
      ! do a match that includes two alternate keys and possibility of = in arg
      integer :: iequal
      character(*) :: key,longkey,arg
      cla_key_arg_match = .false.
      cla_key_arg_match = cla_str_eq(trim(key),trim(arg)) .or. &
                      cla_str_eq(trim(longkey),trim(arg))
      if (cla_key_arg_match) return
      iequal = index(arg,"=")
      if (iequal > 1) &
         cla_key_arg_match = cla_str_eq(trim(key),arg(1:(iequal-1))) .or. &
                         cla_str_eq(trim(longkey),arg(1:(iequal-1)))      
    end function cla_key_arg_match   
    

    logical function cla_str_eq(str1,str2)
      implicit none
      character(*) :: str1, str2
      integer :: str_test
      str_test = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
      cla_str_eq = .false.
      if (str_test /= 0) cla_str_eq = .true.
    end function cla_str_eq
    
    
    subroutine cla_validate(cmd_name)
      implicit none
      character(len=*)      :: cmd_name
      character(len=STRLEN) :: arg
      character(len=STRLEN)  :: value, key
      integer(kind=int_kind) :: ncla, k, kk, iequal
      ncla = cla_command_argument_count()
      if (ncla == 0) return
      
      ! First check for -?, -h, -H, -help, or --help flags.
      call cla_get_command_argument(1,arg)
      key = trim(arg)
      if (cla_str_eq(trim(key),'-h')      .or. &
          cla_str_eq(trim(key),'-?')      .or. &
          cla_str_eq(trim(key),'/?')      .or. &
          cla_str_eq(trim(key),'-H')      .or. &
          cla_str_eq(trim(key),'-help')   .or. &
          cla_str_eq(trim(key),'--help')  .or. &
          cla_str_eq(trim(key),'--usage')      &
          ) then
         call cla_help(cmd_name)
         stop
      endif
    end subroutine cla_validate
    
    logical function cla_key_present(key)
      implicit none
      character(len=STRLEN) :: arg
      character(len=*)  :: key
      character(len=STRLEN) :: longkey
      character(len=2) :: shortkey
      character(len=STRLEN)  :: value
      
      integer(kind=int_kind) :: ncla, k, kk
!      integer :: cla_command_argument_count
!      external cla_command_argument_count

      !     Loop over the command line arguments to assign to
      !     value.
      !     Note that no error is reported if the key was NOT
      !     registered, but it is present on the command line.
 
      cla_key_present = .false.
      

!      print *,'Calling cla_key_present with key = ',trim(key)
      value = trim(cla_empty)
      do kk=1,cla_num
         ! must test for exact match, not just substring
         if (cla_key_arg_match(cla_registry(kk)%key, &
                           cla_registry(kk)%longkey, &
                           key))then
            value = trim(cla_registry(kk)%default)
            longkey = cla_registry(kk)%longkey
            shortkey = cla_registry(kk)%key
            exit
         end if
      end do
      
      if (index(trim(value),trim(cla_empty)) /= 0) then
         call cla_show
         call cla_fatal('Unknown command line argument: '//trim(key))
      endif
      
      ncla = cla_command_argument_count()
      if (ncla == 0) return
      
      do k=1,ncla
         call cla_get_command_argument(k,arg)
         ! test for exact match
         if (cla_key_arg_match(shortkey,longkey,arg))then
            cla_key_present = .true.
            return
         endif
      enddo
      
    end function cla_key_present

    subroutine cla_get_char(key,value)
      implicit none
      character(len=STRLEN)  :: arg
      character(len=*)       :: key
      character(len=2)       :: shortkey
      character(len=STRLEN)  :: longkey
      character(len=STRLEN)  :: value, pvalue
      character(len=STRLEN)  :: kkey
      integer(kind=int_kind) :: ncla, k, kkind, iequal
      integer :: kk, kmatch, ordinal
      logical :: just_matched, prev_matched, is_match, carryover
!      integer :: cla_command_argument_count
!      external cla_command_argument_count

      !     Loop over the command line arguments to assign to
      !     value.
      !     Note that no error is reported if the key was NOT
      !     registered, but it is present on the command line.
      if (index(key,"-") /= 1)then
         ! assume positional argument, confirm by matching name
         ordinal = -1
         do k=1,cla_posarg_num
            ! must test for exact match, not just substring
            if (cla_str_eq(trim(cla_posarg_registry(k)%key),trim(key))) then
               pvalue = trim(cla_posarg_registry(k)%default)
               ordinal = k
            end if
         end do

         ! It seems that value is not yet defined here:
         !         if (index(trim(value),trim(cla_empty)) /= 0) then
         if (ordinal == -1) then
            print *,'Error: You tried to retrieve an unknown command line argument: ',trim(key)
            call cla_show
            stop 5
         endif
         
         if (ordinal > 0) then
            ncla = cla_command_argument_count()
            if (ncla == 0) then
               value=pvalue
               return
            end if
            kmatch = 0
            prev_matched = .False.
            
            do k=1,ncla
               call cla_get_command_argument(k,arg)
               ! test for exact match among key args
               just_matched = .False.    
               do kk = 1, cla_num
                  kkey = cla_registry(kk)%key
                  kkind = cla_registry(kk)%kind
                  is_match = cla_key_arg_match(kkey,cla_registry(kk)%longkey,arg)
                  carryover = kkind/=cla_flag
                  iequal = index(arg,"=")
                  if (is_match .and. iequal > 1)then
                     carryover = .False.
                  end if
                  just_matched = just_matched .or. is_match
                  if (just_matched) exit  ! preserve kkey and kkind
               end do
               if (just_matched .or. prev_matched )then
                  ! current arg part of keyword arg constructal
                  ! so this is not positional
                  prev_matched = just_matched .and. carryover
                  carryover = .False.
                  cycle
               end if
               kmatch = kmatch + 1
               ! increment # of positionals
               if(kmatch == ordinal) then
                  pvalue = trim(arg)
                  value=pvalue
                  return
               end if
               
            end do
         end if
         value = pvalue
         return
      end if
      
      ! keyword
      do k=1,cla_num
         ! must test for exact match, not just substring
         if (cla_key_arg_match(cla_registry(k)%key, &
              cla_registry(k)%longkey, &
              key)) then
            shortkey = cla_registry(k)%key
            longkey = cla_registry(k)%longkey
            
            value = trim(cla_registry(k)%default)
            kkind = cla_registry(k)%kind
         end if
      end do
      
      if (index(trim(value),trim(cla_empty)) /= 0) then
         print *,'Error: You tried to retrieve an unknown command line argument: ',trim(key)
         call cla_show
         stop 5
      endif
      
      ncla = cla_command_argument_count()
      if (ncla == 0) return
      
      do k=1,ncla
         call cla_get_command_argument(k,arg)
         ! test for exact match
         if (cla_key_arg_match(shortkey,longkey,trim(arg))) then
            if (kkind == cla_flag) then
               value = 't'
               return
            else
               iequal = index(arg,"=")
               if (iequal < 1)then
                  call cla_get_command_argument(k+1,arg)       
                  value = trim(arg)
                  return
               else
                  value=arg(iequal+1:len_trim(arg))
                  return
               end if
            endif
         end if
      enddo
      
    end subroutine cla_get_char
    
    
    subroutine cla_get_float_r4(key,float_value)
      implicit none
      character(len=*)       :: key
      character(len=STRLEN)  :: value
      real(kind=4)           :: float_value
      
      
      call cla_get_char(key,value)
      if (index(trim(value),trim(cla_empty)) == 0) read(value,*)float_value
    end subroutine cla_get_float_r4
    
  subroutine cla_get_float_r8(key,float_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    real(kind=8)           :: float_value

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) &
       read(value,*,ERR=100)float_value
    return
 100  call cla_fatal("Input value not correct type: "//value)
  end subroutine cla_get_float_r8


  subroutine cla_get_int_i4(key,int_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    integer(kind=4)        :: int_value
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*,ERR=100)int_value
    return
 100  call cla_fatal("Input value not correct type: "//value)    
    
  end subroutine cla_get_int_i4

  subroutine cla_get_int_i8(key,int_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    integer(kind=8)        :: int_value
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*,ERR=100)int_value
    return
 100  call cla_fatal("Input value not correct type: "//value)        
  end subroutine cla_get_int_i8

  subroutine cla_get_logical(key,logical_value)
    implicit none
    character(len=*)  :: key
    character(len=STRLEN)  :: value
    logical :: logical_value
    integer(kind=int_kind) :: k
    
    logical_value = .false.

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) then
       do k=1,6
          if (index(trim(value),trim(cla_true_str(k))) /= 0) then
             logical_value = .true.
          endif
       end do
    end if
  end subroutine cla_get_logical

  subroutine cla_get_flag(key,logical_value)
    implicit none
    character(len=*)  :: key
    character(len=STRLEN)  :: value
    logical :: logical_value
    integer(kind=int_kind) :: k
    
    logical_value = .false.

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) then
       do k=1,6
          if (index(trim(value),trim(cla_true_str(k))) /= 0) then
             logical_value = .true.
          endif
       end do
    end if
  end subroutine cla_get_flag
  

end module cla
