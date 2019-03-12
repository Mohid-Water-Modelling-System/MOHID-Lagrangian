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


# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group(${MOHIDBase1_LIB_NAME} FILES ${LIB_SRC})

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
	set(proj4_inc_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Proj4/Include)
elseif(LINUX)
  set(hdf5_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/HDF5_1.8.17_x64/Linux/hdf5)
  set(hdf5_LIBS_PATH ${hdf5_PATH}/Lib)
	set(hdf5_MOD_PATH ${hdf5_PATH}/include)
	set(proj4_inc_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Proj4/Linux/proj4-fortran-old/build/include)
endif()

include_directories(${hdf5_LIBS_PATH})
include_directories(${hdf5_MOD_PATH})
include_directories(${proj4_inc_PATH})

message(STATUS "Added HDF5 Libraries path " ${hdf5_PATH})
message(STATUS "Added Proj4 Libraries path " ${proj4_inc_PATH})

# Also adds sources to the Solution Explorer
add_library(${MOHIDBase1_LIB_NAME} ${LIB_SRC})

# Creates folder "libraries" and adds target project (MOHIDLagrangianLib.vcproj)
set_property(TARGET ${MOHIDBase1_LIB_NAME} PROPERTY FOLDER "libraries")
set_property(TARGET ${MOHIDBase1_LIB_NAME} PROPERTY LINKER_LANGUAGE Fortran)


