message("At MOHIDBase2 Library CMAKE file")

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.f *.F *.inc)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group(${MOHIDBase2_LIB_NAME} FILES ${LIB_SRC})

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
    set(hdf5_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64)
	set(hdf5_LIBS_PATH ${hdf5_PATH}/Lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/Include/fortran)
	set(netcdf_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Netcdf_4.4.1.1/VS2015/x64)
elseif(LINUX)
    set(hdf5_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64/Linux/hdf5)
    set(hdf5_LIBS_PATH ${hdf5_PATH}/Lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/include)
	set(netcdf_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Netcdf_4.4.1.1/Linux/build)
    set(netcdfc_inc_PATH ${netcdf_PATH}/include)
    set(netcdfc_lib_PATH ${netcdf_PATH}/lib)
    set(netcdff_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/netcdf-fortran-4.4.4/build)
    set(netcdff_mod_PATH ${netcdff_PATH}/include)
    set(netcdff_lib_PATH ${netcdff_PATH}/lib)
endif()

include_directories(${hdf5_LIBS_PATH})
include_directories(${hdf5_MOD_PATH})
include_directories(${netcdf_PATH})
if(LINUX)
    include_directories(${netcdfc_inc_PATH})
    include_directories(${netcdfc_lib_PATH})
    include_directories(${netcdff_mod_PATH})
    include_directories(${netcdff_lib_PATH})
endif()

message(STATUS "Added HDF5 Libraries path " ${hdf5_PATH})
message(STATUS "Added NetCDF Libraries path " ${netcdf_PATH})
if(LINUX)
    message(STATUS "Added NetCDF-Fortran Libraries path " ${netcdff_PATH})
    #message(STATUS "Added NetCDF include path           " ${netcdfc_inc_PATH})
    #message(STATUS "Added NetCDF lib path               " ${netcdfc_lib_PATH})
    #message(STATUS "Added NetCDF-Fortran include path   " ${netcdff_mod_PATH})
    #message(STATUS "Added NetCDF-Fortran lib path       " ${netcdff_lib_PATH})
endif()

# Also adds sources to the Solution Explorer
add_library(${MOHIDBase2_LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${MOHIDBase2_LIB_NAME} PROPERTY FOLDER "libraries")
set_property(TARGET ${MOHIDBase2_LIB_NAME} PROPERTY LINKER_LANGUAGE Fortran)