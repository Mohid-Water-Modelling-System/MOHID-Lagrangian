#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(VTKFortran ${LIB_SRC})

add_dependencies(VTKFortran ${PENF_PATH})
add_dependencies(VTKFortran ${BeFoR64_PATH})
add_dependencies(VTKFortran ${FoXy_PATH})
add_dependencies(VTKFortran ${StringiFor_PATH})

SET_TARGET_PROPERTIES(VTKFortran PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
