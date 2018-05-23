message("At Library CMAKE file")

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

message(STATUS "Library source files " ${LIB_SRC})

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
if(WINDOWS)
  set(EXTERNAL_LIBS_PATH ${PROJECT_SOURCE_DIR}/Libs/build_win)
  set(EXTERNAL_LIBS_MODULES_PATH ${PROJECT_SOURCE_DIR}/Libs/build_win/modules)
  set(EXTERNAL_LIBS_LIB_PATH ${PROJECT_SOURCE_DIR}/Libs/build_win/lib)
elseif(LINUX)
  set(EXTERNAL_LIBS_PATH ${PROJECT_SOURCE_DIR}/Libs/build_linux)
  set(EXTERNAL_LIBS_MODULES_PATH ${PROJECT_SOURCE_DIR}/Libs/build_linux/modules)
  set(EXTERNAL_LIBS_LIB_PATH ${PROJECT_SOURCE_DIR}/Libs/build_linux/lib)
endif()
include_directories (.)
include_directories(${EXTERNAL_LIBS_PATH})
include_directories(${EXTERNAL_LIBS_MODULES_PATH})
include_directories(${EXTERNAL_LIBS_LIB_PATH})

message(STATUS "Added External Library path " ${EXTERNAL_LIBS_PATH})
message(STATUS "Added External Library module path " ${EXTERNAL_LIBS_MODULES_PATH})
message(STATUS "Added External Library compile path " ${EXTERNAL_LIBS_LIB_PATH})

# Also adds sources to the Solution Explorer
ADD_LIBRARY(${LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "libraries")

# Properties->General->Output Directory
#set_target_properties(${LIB_NAME} PROPERTIES
#                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
