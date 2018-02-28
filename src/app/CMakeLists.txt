# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (MOHIDLagrangian main.f90)

# Properties->Fortran->General->Additional Include Directories
include_directories (${PROJECT_SOURCE_DIR}/src/lib)
include_directories (${PROJECT_SOURCE_DIR}/Libs/build/lib/RELEASE)

# Properties->Linker->Input->Additional Dependencies
target_link_libraries (MOHIDLagrangian  MOHIDLagrangianLib)
target_link_libraries (MOHIDLagrangian  PENF)

# Creates a folder "executables" and adds target 
# project (app.vcproj) under it
set_property(TARGET MOHIDLagrangian PROPERTY FOLDER "executables")

# Properties->General->Output Directory
set_target_properties(MOHIDLagrangian PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
					  
# Adds logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS MOHIDLagrangian
         RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/_install)