# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.13.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.13.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build

# Include any dependencies generated for this target.
include BeFoR64/CMakeFiles/BeFoR64.dir/depend.make

# Include the progress variables for this target.
include BeFoR64/CMakeFiles/BeFoR64.dir/progress.make

# Include the compile flags for this target's objects.
include BeFoR64/CMakeFiles/BeFoR64.dir/flags.make

BeFoR64/CMakeFiles/BeFoR64.dir/befor64.F90.o: BeFoR64/CMakeFiles/BeFoR64.dir/flags.make
BeFoR64/CMakeFiles/BeFoR64.dir/befor64.F90.o: ../BeFoR64/befor64.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object BeFoR64/CMakeFiles/BeFoR64.dir/befor64.F90.o"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64/befor64.F90 -o CMakeFiles/BeFoR64.dir/befor64.F90.o

BeFoR64/CMakeFiles/BeFoR64.dir/befor64.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/BeFoR64.dir/befor64.F90.i"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64/befor64.F90 > CMakeFiles/BeFoR64.dir/befor64.F90.i

BeFoR64/CMakeFiles/BeFoR64.dir/befor64.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/BeFoR64.dir/befor64.F90.s"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64/befor64.F90 -o CMakeFiles/BeFoR64.dir/befor64.F90.s

BeFoR64/CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.o: BeFoR64/CMakeFiles/BeFoR64.dir/flags.make
BeFoR64/CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.o: ../BeFoR64/befor64_pack_data_m.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object BeFoR64/CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.o"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64/befor64_pack_data_m.F90 -o CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.o

BeFoR64/CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.i"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64/befor64_pack_data_m.F90 > CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.i

BeFoR64/CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.s"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && /usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64/befor64_pack_data_m.F90 -o CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.s

# Object files for target BeFoR64
BeFoR64_OBJECTS = \
"CMakeFiles/BeFoR64.dir/befor64.F90.o" \
"CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.o"

# External object files for target BeFoR64
BeFoR64_EXTERNAL_OBJECTS =

lib/libBeFoR64.a: BeFoR64/CMakeFiles/BeFoR64.dir/befor64.F90.o
lib/libBeFoR64.a: BeFoR64/CMakeFiles/BeFoR64.dir/befor64_pack_data_m.F90.o
lib/libBeFoR64.a: BeFoR64/CMakeFiles/BeFoR64.dir/build.make
lib/libBeFoR64.a: BeFoR64/CMakeFiles/BeFoR64.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran static library ../lib/libBeFoR64.a"
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && $(CMAKE_COMMAND) -P CMakeFiles/BeFoR64.dir/cmake_clean_target.cmake
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BeFoR64.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
BeFoR64/CMakeFiles/BeFoR64.dir/build: lib/libBeFoR64.a

.PHONY : BeFoR64/CMakeFiles/BeFoR64.dir/build

BeFoR64/CMakeFiles/BeFoR64.dir/clean:
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 && $(CMAKE_COMMAND) -P CMakeFiles/BeFoR64.dir/cmake_clean.cmake
.PHONY : BeFoR64/CMakeFiles/BeFoR64.dir/clean

BeFoR64/CMakeFiles/BeFoR64.dir/depend:
	cd /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/BeFoR64 /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64 /Users/rbc-laptop/Documents/GitHub/MOHID-Lagrangian/ExternalLibs/Base/build/BeFoR64/CMakeFiles/BeFoR64.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : BeFoR64/CMakeFiles/BeFoR64.dir/depend
