lib_LIBRARIES = libfproj4.a
libfproj4_a_SOURCES = fort-proj.c proj4.f90 proj4.inc cfortran.h
include_HEADERS = proj4.inc 
nodist_include_HEADERS = proj4.mod
check_PROGRAMS = test_proj
test_proj_SOURCES = test-proj.f90 
test_proj_LDADD = libfproj4.a
