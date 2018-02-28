#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group("PENF" FILES ${LIB_SRC})

# Properties->C/C++->General->Additional Include Directories
include_directories (.)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(PENF ${LIB_SRC})

# Creates folder "libraries" and adds target project (math.vcproj)
set_property(TARGET PENF PROPERTY FOLDER "libraries")

# Properties->General->Output Directory
set_target_properties(PENF PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)