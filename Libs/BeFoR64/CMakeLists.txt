#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(BeFoR64 ${LIB_SRC})

add_dependencies(BeFoR64 ${PENF_PATH})

SET_TARGET_PROPERTIES(BeFoR64 PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
