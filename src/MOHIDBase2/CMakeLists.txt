message("At MOHIDBase2 Library CMAKE file")

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.f *.F *.inc)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

set(LOCAL_LIB_NAME MOHIDBase2)

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group(${LOCAL_LIB_NAME} FILES ${LIB_SRC})

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
	set(hdf5_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64/Lib)
	set(hdf5_MOD_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64/Include/fortran)
	set(netcdf_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Netcdf_4.4.1.1/VS2015/x64)
elseif(LINUX)
  
endif()

include_directories(${hdf5_LIBS_PATH})
include_directories(${hdf5_MOD_PATH})
include_directories(${netcdf_PATH})

message(STATUS "Added HDF5 Libraries path " ${hdf5_MOD_PATH})
message(STATUS "Added NetCDF Libraries path " ${netcdf_PATH})

# Also adds sources to the Solution Explorer
add_library(${LOCAL_LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${LOCAL_LIB_NAME} PROPERTY FOLDER "libraries")
set_property(TARGET ${LOCAL_LIB_NAME} PROPERTY LINKER_LANGUAGE Fortran)


