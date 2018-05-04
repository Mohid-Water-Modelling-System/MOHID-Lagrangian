#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group(${LIB_NAME} FILES ${LIB_SRC})

# Check platforms
if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
elseif(CMAKE_HOST_APPLE)
    set(LINUX 1)
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
endif()

# Properties->C/C++->General->Additional Include Directories
include_directories (.)
if(WINDOWS)
  include_directories (${PROJECT_SOURCE_DIR}/Libs)
  include_directories (${PROJECT_SOURCE_DIR}/Libs/build_win/modules)
  include_directories (${PROJECT_SOURCE_DIR}/Libs/build_win/lib)
elseif(LINUX)
  include_directories (${PROJECT_SOURCE_DIR}/Libs)
  include_directories (${PROJECT_SOURCE_DIR}/Libs/build_linux/modules)
  include_directories (${PROJECT_SOURCE_DIR}/Libs/build_linux/lib)
endif()

# Also adds sources to the Solution Explorer
ADD_LIBRARY(${LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "libraries")

# Properties->General->Output Directory
#set_target_properties(${LIB_NAME} PROPERTIES
#                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
