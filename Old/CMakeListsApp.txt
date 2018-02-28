# Properties->Fortran->General->Additional Include Directories
include_directories (${PROJECT_SOURCE_DIR}/src/lib)

# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (MOHID_LagrTracer main.f90)

link_directories(${CMAKE_SOURCE_DIR}/build/lib/RELEASE)
# Properties->Linker->Input->Additional Dependencies
target_link_libraries (MOHID_LagrTracer  LagrTracer)

link_directories(${CMAKE_SOURCE_DIR}/Libs/Build/lib/RELEASE)

# Creates a folder "executables" and adds target 
# project (app.vcproj) under it
set_property(TARGET MOHID_LagrTracer PROPERTY FOLDER "executables")

# Properties->General->Output Directory
set_target_properties(MOHID_LagrTracer PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)