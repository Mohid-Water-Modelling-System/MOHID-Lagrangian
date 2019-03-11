#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(VecFor ${LIB_SRC})

add_dependencies(VecFor ${PENF_PATH})

SET_TARGET_PROPERTIES(VecFor PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
