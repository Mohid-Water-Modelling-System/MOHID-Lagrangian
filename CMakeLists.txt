cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(MOHIDLagrangian Fortran) 
enable_language(Fortran)

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release")
endif()

# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

SET(BIN_NAME MOHIDLagrangian)
SET(LIB_NAME MOHIDLagrangianLib)
set(MOHIDBase1_LIB_NAME MOHIDBase1)
set(MOHIDBase2_LIB_NAME MOHIDBase2)

message("Project name " ${BIN_NAME})
message("Library name " ${LIB_NAME})

#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

#################################################################
# BUILD PATHS
#################################################################

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

message(STATUS "Building paths")
message(STATUS "Executable output path " ${EXECUTABLE_OUTPUT_PATH})
message(STATUS "Library output path " ${LIBRARY_OUTPUT_PATH})
message(STATUS "Modules output path " ${CMAKE_Fortran_MODULE_DIRECTORY})


#################################################################
# FFLAGS depend on the compiler and the build type
#################################################################

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

ADD_DEFINITIONS(-D${CMAKE_Fortran_COMPILER_ID})

message(STATUS "COMPILER INFO: ${CMAKE_Fortran_COMPILER_ID} - ${Fortran_COMPILER_NAME}")

IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran
  set(FORTRAN_FLAGS "-fopenmp -fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -fbounds-check -fprofile-arcs -ftest-coverage -Wimplicit-interface ${EXTRA_FLAGS} " CACHE STRING "" FORCE)

ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort
  set(FORTRAN_FLAGS "-openmp -fpp -heap-arrays -nostandard-realloc-lhs ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-traceback -debug all -check all -ftrapuv -warn nointerfaces ${EXTRA_FLAGS} " CACHE STRING "" FORCE)

ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL" OR Fortran_COMPILER_NAME MATCHES "xlf*")
  # xlf (untested)
  set(FORTRAN_FLAGS "-q64 -qrealsize=8 -qsuffix=f=f90:cpp=f90 ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -qfullpath -qkeepparm ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
ELSE ()
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
ENDIF ()

message (STATUS "Fortran compiler full path: " ${CMAKE_Fortran_COMPILER})
message (STATUS "Fortran compiler flags: " ${CMAKE_Fortran_FLAGS})
message (STATUS "Fortran extra flags: " ${CMAKE_Fortran_FLAGS_RELEASE})

#################################################################
# ADD SOURCE SUBDIRS
#################################################################

SET(LIBRARY_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/MOHIDLagrangianLib)
SET(APP_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/MOHIDLagrangianApp)
set(MOHIDBase1_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/MOHIDBase1)
set(MOHIDBase2_SOURCE_PATH ${CMAKE_SOURCE_DIR}/src/MOHIDBase2)

message(STATUS "Library source path " ${LIBRARY_SOURCE_PATH})
message(STATUS "Application source path " ${APP_SOURCE_PATH})
message(STATUS "MOHIDBase1 source path " ${MOHIDBase1_SOURCE_PATH})
message(STATUS "MOHIDBase2 source path " ${MOHIDBase2_SOURCE_PATH})

ADD_SUBDIRECTORY(${LIBRARY_SOURCE_PATH})
ADD_SUBDIRECTORY(${APP_SOURCE_PATH})
ADD_SUBDIRECTORY(${MOHIDBase1_SOURCE_PATH})
ADD_SUBDIRECTORY(${MOHIDBase2_SOURCE_PATH})


#add dependencies between libs

add_dependencies(${MOHIDBase2_LIB_NAME} ${MOHIDBase1_LIB_NAME})
add_dependencies(${LIB_NAME} ${MOHIDBase1_LIB_NAME})
add_dependencies(${LIB_NAME} ${MOHIDBase2_LIB_NAME})