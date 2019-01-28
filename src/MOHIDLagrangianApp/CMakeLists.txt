message("At Application CMAKE file")

# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer

add_executable(${BIN_NAME} MOHIDLagrangian.f90)

if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
elseif(CMAKE_HOST_APPLE)
    set(LINUX 1)
    set(OSX 1)
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
endif()

# Properties->Fortran->General->Additional Include Directories
include_directories (${LIBRARY_OUTPUT_PATH})
message(STATUS "Added " ${LIB_NAME} " output path " ${LIBRARY_OUTPUT_PATH})

set(base_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/Base/build)
set(datetime_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/datetime-fortran/build)
set(cvs_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/fortran-csv-module/build)
set(FOX_LIBS_PATH ${PROJECT_SOURCE_DIR}/ExternalLibs/fox/build)

include_directories(${base_LIBS_PATH}/modules)
include_directories(${base_LIBS_PATH}/lib)
include_directories(${datetime_LIBS_PATH}/modules)
include_directories(${datetime_LIBS_PATH}/lib)
include_directories(${cvs_LIBS_PATH}/modules)
include_directories(${cvs_LIBS_PATH}/lib)
include_directories(${FOX_LIBS_PATH}/modules)
include_directories(${FOX_LIBS_PATH}/lib)

message(STATUS "Added Base Libraries path " ${base_LIBS_PATH})
message(STATUS "Added DateTime Libraries path " ${datetime_LIBS_PATH})
message(STATUS "Added cvs Libraries path " ${cvs_LIBS_PATH})
message(STATUS "Added FoX Library path " ${FOX_LIBS_PATH})

# Properties->Linker->Input->Additional Dependencies
if(WINDOWS)
  file(GLOB base_LIBS "${base_LIBS_PATH}/lib/RELEASE/*.lib")
  file(GLOB datetime_LIBS "${datetime_LIBS_PATH}/lib/RELEASE/*.lib")
  file(GLOB csv_LIBS "${cvs_LIBS_PATH}/lib/RELEASE/*.lib")
  file(GLOB FOX_LIBS "${FOX_LIBS_PATH}/lib/RELEASE/*.lib")
elseif(LINUX)
  file(GLOB base_LIBS "${base_LIBS_PATH}/lib/*.a")
  file(GLOB datetime_LIBS "${datetime_LIBS_PATH}/lib/*.a")
  file(GLOB csv_LIBS "${cvs_LIBS_PATH}/lib/*.a")
  file(GLOB FOX_LIBS "${FOX_LIBS_PATH}/lib/*.a")
endif()

message(STATUS "Targeting library " ${LIB_NAME})

target_link_libraries (${BIN_NAME}  ${LIB_NAME})
if(WINDOWS)
  target_link_libraries (${BIN_NAME} ${base_LIBS})
  target_link_libraries (${BIN_NAME} ${datetime_LIBS})
  target_link_libraries (${BIN_NAME} ${csv_LIBS})
  target_link_libraries (${BIN_NAME} ${FOX_LIBS})
elseif(OSX)
  target_link_libraries (${BIN_NAME} ${base_LIBS})
  target_link_libraries (${BIN_NAME} ${datetime_LIBS})
  target_link_libraries (${BIN_NAME} ${csv_LIBS})
  target_link_libraries (${BIN_NAME} ${FOX_LIBS})
  #target_link_libraries (nc_reader ${SIMPLE_LIBS})
elseif(LINUX)
  target_link_libraries (${BIN_NAME}  -Wl,--start-group ${base_LIBS} -Wl,--end-group)
  target_link_libraries (${BIN_NAME}  -Wl,--start-group ${datetime_LIBS} -Wl,--end-group)
  target_link_libraries (${BIN_NAME}  -Wl,--start-group ${csv_LIBS} -Wl,--end-group)
  target_link_libraries (${BIN_NAME}  -Wl,--start-group ${FOX_LIBS} -Wl,--end-group)
endif()

# Creates a folder "executables" and adds target project (app.vcproj) under it
set_property(TARGET ${BIN_NAME} PROPERTY FOLDER "executables")
set_property(TARGET ${BIN_NAME} PROPERTY LINKER_LANGUAGE Fortran)

# Properties->General->Output Directory
set_target_properties(${BIN_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)