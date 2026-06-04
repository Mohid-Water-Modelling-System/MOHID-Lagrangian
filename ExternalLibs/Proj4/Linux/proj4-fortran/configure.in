dnl Process this file with autoconf to produce a configure script.
AC_INIT([libfproj4], [1.0],[Magnus.Hagdorn@ed.ac.uk])
AM_INIT_AUTOMAKE(libfproj4, 1.0)
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_FC
dnl do some F95 checking
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[write(*,*) 'Hello world'])],,AC_MSG_FAILURE([Unexpected f95 failure]))
AC_LANG_POP(Fortran)
dnl look for maths lib
AC_CHECK_LIB(m,exp,,,)
dnl look for proj4
AC_ARG_WITH(proj4,AS_HELP_STRING([--with-proj4],[location of where libproj4 is installed]),[proj4=$withval])
if test "$proj4" != ""; then
  AC_MSG_RESULT(proj4 location: $proj4)
  CPPFLAGS="$CPPFLAGS -I$proj4/include"
  LDFLAGS="$LDFLAGS -L$proj4/lib"
fi

AC_CHECK_LIB(proj,pj_init,,AC_MSG_FAILURE([Can't find proj4]))

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
