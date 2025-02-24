message("At MOHIDBase1 Library CMAKE file")

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

# Collect sources into the variable LIB_SRC
FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.f *.F *.inc)

FOREACH(item ${LIB_SRC})
  IF(${item} MATCHES "mpif.F90")
    LIST(REMOVE_ITEM LIB_SRC ${item})
  ENDIF(${item} MATCHES "mpif.F90")
ENDFOREACH(item)
FOREACH(item ${LIB_SRC})
  IF(${item} MATCHES "ModulePhreeqC.F90")
    LIST(REMOVE_ITEM LIB_SRC ${item})
  ENDIF(${item} MATCHES "ModulePhreeqC.F90")
ENDFOREACH(item)
FOREACH(item ${LIB_SRC})
  IF(${item} MATCHES "ModulePhreeqCRM.F90")
    LIST(REMOVE_ITEM LIB_SRC ${item})
  ENDIF(${item} MATCHES "ModulePhreeqCRM.F90")
ENDFOREACH(item)
FOREACH(item ${LIB_SRC})
  IF(${item} MATCHES "IPhreeqc_interface.F90")
    LIST(REMOVE_ITEM LIB_SRC ${item})
  ENDIF(${item} MATCHES "IPhreeqc_interface.F90")
ENDFOREACH(item)
FOREACH(item ${LIB_SRC})
  IF(${item} MATCHES "ModuleCuda.F90")
    LIST(REMOVE_ITEM LIB_SRC ${item})
  ENDIF(${item} MATCHES "ModuleCuda.F90")
ENDFOREACH(item)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

# Also adds sources to the Solution Explorer
add_library(MOHIDBase1 ${LIB_SRC})

# Check platforms
if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
elseif(CMAKE_HOST_APPLE)
    set(LINUX 1)
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
endif()

# Properties->C/C++->General->Additional Include Directories
get_filename_component(ExternalLibsPath ../../ ABSOLUTE)
if(WINDOWS)
	set(hdf5_PATH ${ExternalLibsPath}/HDF5_1.8.17_x64)
	set(hdf5_DLL_PATH ${hdf5_PATH}/bin)
	set(hdf5_LIBS_PATH ${hdf5_PATH}/Lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/Include/fortran)
	set(proj4_inc_PATH ${ExternalLibsPath}/Proj4/Include)
	set(netcdf_PATH ${ExternalLibsPath}/Netcdf_4.4.1.1/VS2015/x64)
elseif(LINUX)
	set(zlib_PATH ${ExternalLibsPath}/zlib-1.2.11/build)
	set(zlib_LIBS_PATH ${zlib_PATH}/lib)
	set(zlib_MOD_PATH ${zlib_PATH}/include)
	#set(curl_PATH ${ExternalLibsPath}/curl-7.26.0/build)
	#set(curl_LIBS_PATH ${curl_PATH}/lib)
	#set(curl_MOD_PATH ${curl_PATH}/include/curl)
	set(hdf5_PATH ${ExternalLibsPath}/HDF5_1.8.17_x64/Linux/hdf5)
	set(hdf5_LIBS_PATH ${hdf5_PATH}/lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/include)
	set(proj4_inc_PATH ${ExternalLibsPath}/Proj4/Linux/proj4-fortran/build/include)
	set(netcdf_PATH ${ExternalLibsPath}/Netcdf_4.4.1.1/Linux/build)
	set(netcdfc_inc_PATH ${netcdf_PATH}/include)
	set(netcdfc_lib_PATH ${netcdf_PATH}/lib)
	set(netcdff_PATH ${ExternalLibsPath}/netcdf-fortran-4.4.4/build)
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
#if(WINDOWS)
	#include_directories(${hdf5_DLL_PATH})
#endif()

# Properties->Linker->Input->Additional Dependencies
if(WINDOWS)  
  file(GLOB HDF5_LIBS "${hdf5_PATH}/Lib/*.lib")
  file(GLOB HDF5_DLL_LIBS "${hdf5_DLL_PATH}/*.dll")
  file(GLOB NETCDF_LIBS "${netcdf_PATH}/*.lib")
  file(GLOB NETCDF_DLL_LIBS "${netcdf_PATH}/*.dll")
elseif(LINUX)  
  file(GLOB HDF5_LIBS "${hdf5_PATH}/lib/*.a")
  file(GLOB ZLIB_LIBS "${zlib_PATH}/lib/*.a")
  file(GLOB NETCDFC_LIBS "${netcdfc_lib_PATH}/*.a")
  file(GLOB NETCDFF_LIBS "${netcdff_lib_PATH}/*.a")
endif()

if(WINDOWS)
  target_link_libraries (MOHIDBase1 ${NETCDF_LIBS})
elseif(OSX)  
  target_link_libraries (MOHIDBase1 ${ZLIB_LIBS})
  #target_link_libraries (MOHIDBase1 ${CURL_LIBS})
  target_link_libraries (MOHIDBase1 ${HDF5_LIBS})
  target_link_libraries (MOHIDBase1 ${NETCDFC_LIBS})
  target_link_libraries (MOHIDBase1 ${NETCDFF_LIBS})
elseif(LINUX)
  target_link_libraries (MOHIDBase1  -Wl,--start-group ${ZLIB_LIBS} -Wl,--end-group)
  #target_link_libraries (MOHIDBase1  -Wl,--start-group ${CURL_LIBS} -Wl,--end-group)
  target_link_libraries (MOHIDBase1  -Wl,--start-group ${HDF5_LIBS} -Wl,--end-group)
  target_link_libraries (MOHIDBase1  -Wl,--start-group ${NETCDFC_LIBS} -Wl,--end-group)
  target_link_libraries (MOHIDBase1 dl )
  target_link_libraries (MOHIDBase1 -Wl,--start-group ${ZLIB_LIBS} ${HDF5_LIBS} ${NETCDFC_LIBS} ${NETCDFF_LIBS} -Wl,--end-group)
endif()