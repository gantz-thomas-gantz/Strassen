# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/cmake/1425/bin/cmake

# The command to remove a file.
RM = /snap/cmake/1425/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/janmarxen/MU4IN901/MOD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/janmarxen/MU4IN901/MOD/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/codegen:
.PHONY : CMakeFiles/main.dir/codegen

CMakeFiles/main.dir/src/main.c.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/main.c.o: /home/janmarxen/MU4IN901/MOD/src/main.c
CMakeFiles/main.dir/src/main.c.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/janmarxen/MU4IN901/MOD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/main.dir/src/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/main.dir/src/main.c.o -MF CMakeFiles/main.dir/src/main.c.o.d -o CMakeFiles/main.dir/src/main.c.o -c /home/janmarxen/MU4IN901/MOD/src/main.c

CMakeFiles/main.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/main.dir/src/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/janmarxen/MU4IN901/MOD/src/main.c > CMakeFiles/main.dir/src/main.c.i

CMakeFiles/main.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/main.dir/src/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/janmarxen/MU4IN901/MOD/src/main.c -o CMakeFiles/main.dir/src/main.c.s

CMakeFiles/main.dir/src/IO.c.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/IO.c.o: /home/janmarxen/MU4IN901/MOD/src/IO.c
CMakeFiles/main.dir/src/IO.c.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/janmarxen/MU4IN901/MOD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/main.dir/src/IO.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/main.dir/src/IO.c.o -MF CMakeFiles/main.dir/src/IO.c.o.d -o CMakeFiles/main.dir/src/IO.c.o -c /home/janmarxen/MU4IN901/MOD/src/IO.c

CMakeFiles/main.dir/src/IO.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/main.dir/src/IO.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/janmarxen/MU4IN901/MOD/src/IO.c > CMakeFiles/main.dir/src/IO.c.i

CMakeFiles/main.dir/src/IO.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/main.dir/src/IO.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/janmarxen/MU4IN901/MOD/src/IO.c -o CMakeFiles/main.dir/src/IO.c.s

CMakeFiles/main.dir/src/naive_matmat.c.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/naive_matmat.c.o: /home/janmarxen/MU4IN901/MOD/src/naive_matmat.c
CMakeFiles/main.dir/src/naive_matmat.c.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/janmarxen/MU4IN901/MOD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/main.dir/src/naive_matmat.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/main.dir/src/naive_matmat.c.o -MF CMakeFiles/main.dir/src/naive_matmat.c.o.d -o CMakeFiles/main.dir/src/naive_matmat.c.o -c /home/janmarxen/MU4IN901/MOD/src/naive_matmat.c

CMakeFiles/main.dir/src/naive_matmat.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/main.dir/src/naive_matmat.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/janmarxen/MU4IN901/MOD/src/naive_matmat.c > CMakeFiles/main.dir/src/naive_matmat.c.i

CMakeFiles/main.dir/src/naive_matmat.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/main.dir/src/naive_matmat.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/janmarxen/MU4IN901/MOD/src/naive_matmat.c -o CMakeFiles/main.dir/src/naive_matmat.c.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/src/main.c.o" \
"CMakeFiles/main.dir/src/IO.c.o" \
"CMakeFiles/main.dir/src/naive_matmat.c.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/src/main.c.o
main: CMakeFiles/main.dir/src/IO.c.o
main: CMakeFiles/main.dir/src/naive_matmat.c.o
main: CMakeFiles/main.dir/build.make
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/janmarxen/MU4IN901/MOD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main
.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/janmarxen/MU4IN901/MOD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/janmarxen/MU4IN901/MOD /home/janmarxen/MU4IN901/MOD /home/janmarxen/MU4IN901/MOD/build /home/janmarxen/MU4IN901/MOD/build /home/janmarxen/MU4IN901/MOD/build/CMakeFiles/main.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/main.dir/depend

