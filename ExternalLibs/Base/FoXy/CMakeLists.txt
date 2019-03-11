#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(Foxy ${LIB_SRC})

add_dependencies(Foxy ${PENF_PATH})
add_dependencies(Foxy ${StringiFor_PATH})

SET_TARGET_PROPERTIES(Foxy PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
