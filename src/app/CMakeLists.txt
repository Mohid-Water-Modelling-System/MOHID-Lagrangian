message("At Application CMAKE file")

# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer

add_executable(${BIN_NAME} MOHIDLagrangian.f90)
add_executable(nc_reader read_nc.f90)

if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
elseif(CMAKE_HOST_APPLE)
    set(LINUX 1)
    set(OSX 1)
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
endif()

# Properties->Fortran->General->Additional Include Directories
#include_directories (${LIB_PATH})
include_directories (${LIBRARY_OUTPUT_PATH})

message(STATUS "Added " ${LIB_NAME} " output path " ${LIBRARY_OUTPUT_PATH})

if(WINDOWS)
  set(EXTERNAL_LIBS_MODULES ${PROJECT_SOURCE_DIR}/Libs/Libs_build_win/modules)
  set(EXTERNAL_LIBS_LIB ${PROJECT_SOURCE_DIR}/Libs/Libs_build_win/lib)
elseif(LINUX)
  set(EXTERNAL_LIBS_MODULES ${PROJECT_SOURCE_DIR}/Libs/Libs_build_linux/modules)
  set(EXTERNAL_LIBS_LIB ${PROJECT_SOURCE_DIR}/Libs/Libs_build_linux/lib)
endif()

include_directories (${EXTERNAL_LIBS_MODULES})
include_directories (${EXTERNAL_LIBS_LIB})

message(STATUS "Added External Library module path " ${EXTERNAL_LIBS_MODULES})
message(STATUS "Added External Library .lib path " ${EXTERNAL_LIBS_LIB})

# Properties->Linker->Input->Additional Dependencies

if(WINDOWS)
  file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/RELEASE/*.lib")
  file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/*.lib")
  file(GLOB EXTERNAL_LIBS "${PROJECT_SOURCE_DIR}/Libs/Libs_build_win/lib/RELEASE/*.lib")
elseif(LINUX)
  file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/RELEASE/*.a")
  file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/*.a")
  file(GLOB EXTERNAL_LIBS "${PROJECT_SOURCE_DIR}/Libs/Libs_build_linux/lib/*.a")
endif()

message(STATUS "Targeting library " ${LIB_NAME})
message(STATUS "Targeting " ${LIBS})
message(STATUS "Targeting external libraries " ${EXTERNAL_LIBS})

target_link_libraries (${BIN_NAME}  ${LIB_NAME})
target_link_libraries (nc_reader  ${LIB_NAME})
if(WINDOWS)
  target_link_libraries (${BIN_NAME} ${LIBS})
  target_link_libraries (${BIN_NAME} ${EXTERNAL_LIBS})
  target_link_libraries (nc_reader ${LIBS})
  target_link_libraries (nc_reader ${EXTERNAL_LIBS})
elseif(OSX)
  target_link_libraries (${BIN_NAME} ${LIBS})
  target_link_libraries (${BIN_NAME} ${EXTERNAL_LIBS})
  target_link_libraries (nc_reader ${LIBS})
  target_link_libraries (nc_reader ${EXTERNAL_LIBS})
elseif(LINUX)
  target_link_libraries (${BIN_NAME}  -Wl,--start-group ${LIBS} -Wl,--end-group)
  target_link_libraries (${BIN_NAME}  -Wl,--start-group ${EXTERNAL_LIBS} -Wl,--end-group)
  target_link_libraries (nc_reader -Wl,--start-group ${LIBS} -Wl,--end-group)
  target_link_libraries (nc_reader -Wl,--start-group ${EXTERNAL_LIBS} -Wl,--end-group)
endif()

# Creates a folder "executables" and adds target
# project (app.vcproj) under it
set_property(TARGET ${BIN_NAME} PROPERTY FOLDER "executables")
set_property(TARGET ${BIN_NAME} PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET nc_reader PROPERTY FOLDER "executables")
set_property(TARGET nc_reader PROPERTY LINKER_LANGUAGE Fortran)

# Properties->General->Output Directory
set_target_properties(${BIN_NAME} PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
					  
set_target_properties(nc_reader PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)