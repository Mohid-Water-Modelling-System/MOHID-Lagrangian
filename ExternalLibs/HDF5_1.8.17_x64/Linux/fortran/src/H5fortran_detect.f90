!****h* ROBODoc/H5fortran_detect.f90
!
! NAME
!  H5fortran_detect
! 
! PURPOSE
!  This stand alone program is used at build time to generate the header file
!  H5fort_type_defines.h. The source code itself was automatically generated by
!  the program H5test_kind_STORAGE_SIZE.f90
!
! NOTES
!  This source code makes use of the Fortran intrinsic function STORAGE_SIZE because
!  the availability of the intrinsic function was determined to be available at
!  configure time
!
! COPYRIGHT
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the files COPYING and Copyright.html.  COPYING can be found at the root   *
!   of the source code distribution tree; Copyright.html can be found at the  *
!   root level of an installed copy of the electronic HDF5 document set and   *
!   is linked from the top-level documents page.  It can also be found at     *
!   http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
!   access to either file, you may request a copy from help@hdfgroup.org.     *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! AUTHOR
!  H5test_kind_C_SIZEOF.f90
!
!*****

 MODULE H5test_kind_STORAGE_SIZE_mod
 USE ISO_C_BINDING
 IMPLICIT NONE
 CONTAINS
 SUBROUTINE i00()
    IMPLICIT NONE
    INTEGER :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 4
    WRITE(*,*) "#define H5_FORTRAN_HAS_NATIVE_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE r00()
    IMPLICIT NONE
    REAL :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 4
    WRITE(*,*) "#define H5_FORTRAN_HAS_REAL_NATIVE_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE d00()
    IMPLICIT NONE
    DOUBLE PRECISION :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 8
    WRITE(*,*) "#define H5_FORTRAN_HAS_DOUBLE_NATIVE_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE i01()
    IMPLICIT NONE
   INTEGER(KIND=1) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 1
    WRITE(*,*) "#define H5_FORTRAN_HAS_INTEGER_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE i02()
    IMPLICIT NONE
   INTEGER(KIND=2) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 2
    WRITE(*,*) "#define H5_FORTRAN_HAS_INTEGER_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE i04()
    IMPLICIT NONE
   INTEGER(KIND=4) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 4
    WRITE(*,*) "#define H5_FORTRAN_HAS_INTEGER_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE i08()
    IMPLICIT NONE
   INTEGER(KIND=8) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 8
    WRITE(*,*) "#define H5_FORTRAN_HAS_INTEGER_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE i16()
    IMPLICIT NONE
   INTEGER(KIND=16) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 16
    WRITE(*,*) "#define H5_FORTRAN_HAS_INTEGER_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE r04()
    IMPLICIT NONE
    REAL(KIND= 4) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 4
    WRITE(*,*) "#define H5_FORTRAN_HAS_REAL_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE r08()
    IMPLICIT NONE
    REAL(KIND= 8) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 8
    WRITE(*,*) "#define H5_FORTRAN_HAS_REAL_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE r10()
    IMPLICIT NONE
    REAL(KIND= 10) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 10
    WRITE(*,*) "#define H5_FORTRAN_HAS_REAL_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 SUBROUTINE r16()
    IMPLICIT NONE
    REAL(KIND= 16) :: a
    INTEGER(C_SIZE_T) :: a_size
    CHARACTER(LEN=2) :: ichr2, jchr2
    a_size = STORAGE_SIZE(a, c_size_t)/STORAGE_SIZE(c_char_'a',c_size_t)
    WRITE(ichr2,'(I2)') a_size
    WRITE(jchr2,'(I2)') 16
    WRITE(*,*) "#define H5_FORTRAN_HAS_REAL_"//TRIM(ADJUSTL(ichr2))//"_KIND "//ADJUSTL(jchr2)
    RETURN
 END SUBROUTINE
 END MODULE H5test_kind_STORAGE_SIZE_mod
 
 PROGRAM H5test_kind_STORAGE_SIZE
 USE H5test_kind_STORAGE_SIZE_mod
 WRITE(*,*) " /*generating header file*/ "
 CALL i00()
 CALL r00()
 CALL d00()
 CALL i01()
 CALL i02()
 CALL i04()
 CALL i08()
 CALL i16()
 CALL r04()
 CALL r08()
 CALL r10()
 CALL r16()
 END PROGRAM H5test_kind_STORAGE_SIZE