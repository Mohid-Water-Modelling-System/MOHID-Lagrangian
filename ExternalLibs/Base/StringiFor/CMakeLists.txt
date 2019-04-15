#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(StringiFor ${LIB_SRC})

add_dependencies(StringiFor ${PENF_PATH})
add_dependencies(StringiFor ${BeFoR64_PATH})

SET_TARGET_PROPERTIES(StringiFor PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
