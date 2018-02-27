#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(FoXy ${LIB_SRC})

add_dependencies(FoXy ${PENF_PATH})
add_dependencies(FoXy ${StringiFor_PATH})

SET_TARGET_PROPERTIES(FoXy PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
