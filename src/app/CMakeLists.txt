# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (MOHIDLagrangian main.f90)

# Properties->Fortran->General->Additional Include Directories
include_directories (${LIB_PATH})
include_directories (${LIBRARY_OUTPUT_PATH})
set(EXTERNAL_LIBS_MODULES ${CMAKE_SOURCE_DIR}/Libs/build/modules)
include_directories (${EXTERNAL_LIBS_MODULES})

# Properties->Linker->Input->Additional Dependencies
file(GLOB LIBS "${LIBRARY_OUTPUT_PATH}/RELEASE/*.lib") #This needs to be changed to enable library debug
file(GLOB EXTERNAL_LIBS "${CMAKE_SOURCE_DIR}/Libs/build/lib/RELEASE/*.lib") #sub-ideal linking solution, want to avoid external libs at application level
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