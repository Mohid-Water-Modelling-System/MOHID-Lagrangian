message("At Library CMAKE file")

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

FILE(GLOB_RECURSE C_LIB_SRC *.C *.cpp *.h)
SET(C_LIB_SRC ${C_LIB_SRC} PARENT_SCOPE)

#message(STATUS "Library source files are " ${LIB_SRC})

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
  set(SIMPLE_LIBS_PATH ${PROJECT_SOURCE_DIR}/Libs/Libs_build_win)  
  set(FOX_LIBS_PATH ${PROJECT_SOURCE_DIR}/Libs/fox/build_win)  
elseif(LINUX)
  set(SIMPLE_LIBS_PATH ${PROJECT_SOURCE_DIR}/Libs/Libs_build_linux)
  set(FOX_LIBS_PATH ${PROJECT_SOURCE_DIR}/Libs/fox/build_linux)  
endif()
set(SIMPLE_LIBS_MODULES_PATH ${SIMPLE_LIBS_PATH}/modules)
set(SIMPLE_LIBS_LIB_PATH ${SIMPLE_LIBS_PATH}/lib)
set(FOX_LIBS_MODULES_PATH ${FOX_LIBS_PATH}/modules)
set(FOX_LIBS_LIB_PATH ${FOX_LIBS_PATH}/lib)

include_directories(.)
include_directories(${SIMPLE_LIBS_PATH})
include_directories(${SIMPLE_LIBS_MODULES_PATH})
include_directories(${SIMPLE_LIBS_LIB_PATH})
include_directories(${FOX_LIBS_PATH})
include_directories(${FOX_LIBS_MODULES_PATH})
include_directories(${FOX_LIBS_LIB_PATH})

message(STATUS "Added simple Libraries path " ${SIMPLE_LIBS_PATH})
#message(STATUS "Added SIMPLE Library module path " ${SIMPLE_LIBS_MODULES_PATH})
#message(STATUS "Added SIMPLE Library compile path " ${SIMPLE_LIBS_LIB_PATH})
message(STATUS "Added FoX Library path " ${FOX_LIBS_PATH})

# Also adds sources to the Solution Explorer
add_library(${LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "libraries")
set_property(TARGET ${LIB_NAME} PROPERTY LINKER_LANGUAGE Fortran)


