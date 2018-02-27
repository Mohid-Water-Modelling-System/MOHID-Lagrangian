#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(PENF ${LIB_SRC})

SET_TARGET_PROPERTIES(PENF PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
