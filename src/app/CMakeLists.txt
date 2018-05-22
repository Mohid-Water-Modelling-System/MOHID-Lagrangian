# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (MOHIDLagrangian MOHIDLagrangian.f90)

if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
elseif(CMAKE_HOST_APPLE)
    set(LINUX 1)
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
endif()

# Properties->Fortran->General->Additional Include Directories
#include_directories (${LIB_PATH})
include_directories (${LIBRARY_OUTPUT_PATH})

if(WINDOWS)
  set(EXTERNAL_LIBS_MODULES ${CMAKE_SOURCE_DIR}/Libs/build_win/modules)
  set(EXTERNAL_LIBS_LIB ${CMAKE_SOURCE_DIR}/Libs/build_win/lib)
elseif(LINUX)
  set(EXTERNAL_LIBS_MODULES ${CMAKE_SOURCE_DIR}/Libs/build_linux/modules)
  set(EXTERNAL_LIBS_LIB ${CMAKE_SOURCE_DIR}/Libs/build_linux/lib)
endif()

include_directories (${EXTERNAL_LIBS_MODULES})
include_directories (${EXTERNAL_LIBS_LIB})

# Properties->Linker->Input->Additional Dependencies
file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/RELEASE/*.lib")
file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/DEBUG/*.lib")
file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/*.lib")
if(WINDOWS)
  file(GLOB EXTERNAL_LIBS "${CMAKE_SOURCE_DIR}/Libs/build_win/lib/RELEASE/*.lib")
elseif(LINUX)
  file(GLOB EXTERNAL_LIBS "${CMAKE_SOURCE_DIR}/Libs/build_linux/lib/*.lib")
endif()

target_link_libraries (MOHIDLagrangian  ${LIB_NAME})
target_link_libraries (MOHIDLagrangian  ${LIBS})
target_link_libraries (MOHIDLagrangian  ${EXTERNAL_LIBS})

# Creates a folder "executables" and adds target
# project (app.vcproj) under it
set_property(TARGET MOHIDLagrangian PROPERTY FOLDER "executables")

# Properties->General->Output Directory
set_target_properties(MOHIDLagrangian PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Adds logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS MOHIDLagrangian
         RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/_install)
