#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group("LagrTracerLib" FILES ${LIB_SRC})

#################################################################
# Library target
#################################################################
# Also adds sources to the Solution Explorer
ADD_LIBRARY(${LIB} ${LIB_SRC})

# Properties->General->Output Directory
SET_TARGET_PROPERTIES(${LIB} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
