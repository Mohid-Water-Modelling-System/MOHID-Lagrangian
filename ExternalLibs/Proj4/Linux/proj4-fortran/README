proj4-fortran
=============
This directory contains f77 and f90 wrappers for proj4, a cartograohic
projections library. The C library is wrapped using cfortran.h (which is
hopefully reasonably platform independent). To use these bindings you need

* proj4, which you can get from http://www.remotesensing.org/proj/
* GNU make, a C compiler and a f90 compiler

COMPILE
The following environment variables should be set:
CC:          the C compiler
F90:         the F90 compiler
PROJ_PREFIX: path to where the proj4 library is installed
CFLAGS:      C compiler flags
             you need to specify which F90 compiler you use,
             e.g. for NAG it is -DNAGf90Fortran
             see cfortran.h
FFLAGS:      F90 compiler flags

The library is compiled using make. You can also compile the test program
make test-proj

There is no API documentation yet (ever?). Have a look at the test program
test-proj.f90 and the proj4 documentation.