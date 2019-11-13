message("At Library CMAKE file")

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

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
set(base_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Base/build)
set(datetime_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/datetime-fortran/build)
set(cvs_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/fortran-csv-module/build)
set(FOX_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/fox/build)
set(MOHIDBase_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/MOHID/build)
if(WINDOWS)
	set(hdf5_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64)
	set(hdf5_DLL_PATH ${hdf5_PATH}/bin)
	set(hdf5_LIBS_PATH ${hdf5_PATH}/Lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/Include/fortran)
	set(proj4_inc_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Proj4/Include)
	set(netcdf_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Netcdf_4.4.1.1/VS2015/x64)
elseif(LINUX)
    set(zlib_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/zlib-1.2.11/build)
    set(zlib_LIBS_PATH ${zlib_PATH}/lib)
    set(zlib_MOD_PATH ${zlib_PATH}/include)
    #set(curl_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/curl-7.26.0/build)
    #set(curl_LIBS_PATH ${curl_PATH}/lib)
    #set(curl_MOD_PATH ${curl_PATH}/include/curl)
	set(hdf5_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64/Linux/hdf5)
	set(hdf5_LIBS_PATH ${hdf5_PATH}/lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/include)
	set(proj4_inc_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Proj4/Linux/proj4-fortran/build/include)
	set(netcdf_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Netcdf_4.4.1.1/Linux/build)
    set(netcdfc_inc_PATH ${netcdf_PATH}/include)
    set(netcdfc_lib_PATH ${netcdf_PATH}/lib)
    set(netcdff_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/netcdf-fortran-4.4.4/build)
    set(netcdff_mod_PATH ${netcdff_PATH}/include)
    set(netcdff_lib_PATH ${netcdff_PATH}/lib)
endif()

include_directories(${zlib_LIBS_PATH})
include_directories(${zlib_MOD_PATH})
#include_directories(${curl_LIBS_PATH})
#include_directories(${curl_MOD_PATH})
include_directories(${hdf5_LIBS_PATH})
include_directories(${hdf5_MOD_PATH})
include_directories(${proj4_inc_PATH})
include_directories(${netcdf_PATH})
if(WINDOWS)
	include_directories(${hdf5_DLL_PATH})
elseif(LINUX)
    include_directories(${netcdfc_inc_PATH})
    include_directories(${netcdfc_lib_PATH})
    include_directories(${netcdff_mod_PATH})
    include_directories(${netcdff_lib_PATH})
endif()

include_directories(${base_LIBS_PATH}/modules)
include_directories(${base_LIBS_PATH}/lib)
include_directories(${datetime_LIBS_PATH}/modules)
include_directories(${datetime_LIBS_PATH}/lib)
include_directories(${cvs_LIBS_PATH}/modules)
include_directories(${cvs_LIBS_PATH}/lib)
include_directories(${FOX_LIBS_PATH}/modules)
include_directories(${FOX_LIBS_PATH}/lib)
include_directories(${MOHIDBase_PATH}/modules)
include_directories(${MOHIDBase_PATH}/lib)

message(STATUS "Added Base Libraries path " ${base_LIBS_PATH})
message(STATUS "Added DateTime Libraries path " ${datetime_LIBS_PATH})
message(STATUS "Added cvs Libraries path " ${cvs_LIBS_PATH})
message(STATUS "Added FoX Library path " ${FOX_LIBS_PATH})
message(STATUS "Added MOHID Library path " ${MOHIDBase_PATH})

# Also adds sources to the Solution Explorer
add_library(${LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "libraries")
set_property(TARGET ${LIB_NAME} PROPERTY LINKER_LANGUAGE Fortran)

