# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cczhong/Codes/MANA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cczhong/Codes/MANA

# Include any dependencies generated for this target.
include CMakeFiles/bin/grasp2-GenGraph.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bin/grasp2-GenGraph.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bin/grasp2-GenGraph.dir/flags.make

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.o: GenGraph/main_GenGraph.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.o -c /home/cczhong/Codes/MANA/GenGraph/main_GenGraph.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/main_GenGraph.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/main_GenGraph.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.o: GenGraph/overlap.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.o -c /home/cczhong/Codes/MANA/GenGraph/overlap.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/overlap.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/overlap.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.o: GenGraph/sfa_build.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.o -c /home/cczhong/Codes/MANA/GenGraph/sfa_build.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/sfa_build.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/sfa_build.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.o: GenGraph/gsa.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.o -c /home/cczhong/Codes/MANA/GenGraph/gsa.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/gsa.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/gsa.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.o: GenGraph/sfa.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.o -c /home/cczhong/Codes/MANA/GenGraph/sfa.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/sfa.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/sfa.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.o: GenGraph/interval_array.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.o -c /home/cczhong/Codes/MANA/GenGraph/interval_array.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/interval_array.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/interval_array.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.o: GenGraph/database_index.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.o -c /home/cczhong/Codes/MANA/GenGraph/database_index.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/database_index.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/database_index.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.o: GenGraph/kmer_filtering.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.o -c /home/cczhong/Codes/MANA/GenGraph/kmer_filtering.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/kmer_filtering.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/kmer_filtering.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.o: GenGraph/sequence_search.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.o -c /home/cczhong/Codes/MANA/GenGraph/sequence_search.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/sequence_search.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/sequence_search.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.o: GenGraph/align_batch.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.o -c /home/cczhong/Codes/MANA/GenGraph/align_batch.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/align_batch.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/align_batch.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.o: GenGraph/string_graph.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.o -c /home/cczhong/Codes/MANA/GenGraph/string_graph.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/string_graph.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/string_graph.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.o: GenGraph/bwt.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.o -c /home/cczhong/Codes/MANA/GenGraph/bwt.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/bwt.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/bwt.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.o: GenGraph/bio_alphabet.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.o -c /home/cczhong/Codes/MANA/GenGraph/bio_alphabet.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/bio_alphabet.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/bio_alphabet.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.o: GenGraph/loader.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.o -c /home/cczhong/Codes/MANA/GenGraph/loader.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/loader.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/loader.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.o: GenGraph/concatenator.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.o -c /home/cczhong/Codes/MANA/GenGraph/concatenator.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/concatenator.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/concatenator.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.o: GenGraph/contig_refinement.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.o -c /home/cczhong/Codes/MANA/GenGraph/contig_refinement.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/contig_refinement.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/contig_refinement.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.o: GenGraph/bwt_search.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.o -c /home/cczhong/Codes/MANA/GenGraph/bwt_search.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/bwt_search.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/bwt_search.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.o: GenGraph/kmer_unitcoder.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.o -c /home/cczhong/Codes/MANA/GenGraph/kmer_unitcoder.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/kmer_unitcoder.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/kmer_unitcoder.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.o: GenGraph/minimizer_sort.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_19) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.o -c /home/cczhong/Codes/MANA/GenGraph/minimizer_sort.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/minimizer_sort.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/minimizer_sort.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.o: GenGraph/divsufsort.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_20) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.o -c /home/cczhong/Codes/MANA/GenGraph/divsufsort.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/divsufsort.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/divsufsort.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.o: GenGraph/sssort.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_21) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.o -c /home/cczhong/Codes/MANA/GenGraph/sssort.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/sssort.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/sssort.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.o: GenGraph/trsort.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_22) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.o -c /home/cczhong/Codes/MANA/GenGraph/trsort.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/trsort.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/trsort.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.s

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.o: CMakeFiles/bin/grasp2-GenGraph.dir/flags.make
CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.o: GenGraph/utils.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_23) "Building CXX object CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.o"
	/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.o -c /home/cczhong/Codes/MANA/GenGraph/utils.cc

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/cczhong/Codes/MANA/GenGraph/utils.cc > CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.i

CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/cczhong/Codes/MANA/GenGraph/utils.cc -o CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.s

# Object files for target bin/grasp2-GenGraph
bin/grasp2__GenGraph_OBJECTS = \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.o" \
"CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.o"

# External object files for target bin/grasp2-GenGraph
bin/grasp2__GenGraph_EXTERNAL_OBJECTS =

bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/main_GenGraph.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/overlap.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa_build.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/gsa.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sfa.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/interval_array.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/database_index.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_filtering.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sequence_search.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/align_batch.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/string_graph.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bio_alphabet.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/loader.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/concatenator.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/contig_refinement.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/bwt_search.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/kmer_unitcoder.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/minimizer_sort.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/divsufsort.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/sssort.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/trsort.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/GenGraph/utils.cc.o
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/build.make
bin/grasp2-GenGraph: /usr/lib64/libboost_iostreams-mt.so
bin/grasp2-GenGraph: /usr/lib64/libboost_filesystem-mt.so
bin/grasp2-GenGraph: /usr/lib64/libboost_system-mt.so
bin/grasp2-GenGraph: /usr/lib64/libboost_program_options-mt.so
bin/grasp2-GenGraph: /usr/lib64/libboost_regex-mt.so
bin/grasp2-GenGraph: CMakeFiles/bin/grasp2-GenGraph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/cczhong/Codes/MANA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_24) "Linking CXX executable bin/grasp2-GenGraph"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bin/grasp2-GenGraph.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bin/grasp2-GenGraph.dir/build: bin/grasp2-GenGraph

.PHONY : CMakeFiles/bin/grasp2-GenGraph.dir/build

CMakeFiles/bin/grasp2-GenGraph.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bin/grasp2-GenGraph.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bin/grasp2-GenGraph.dir/clean

CMakeFiles/bin/grasp2-GenGraph.dir/depend:
	cd /home/cczhong/Codes/MANA && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cczhong/Codes/MANA /home/cczhong/Codes/MANA /home/cczhong/Codes/MANA /home/cczhong/Codes/MANA /home/cczhong/Codes/MANA/CMakeFiles/bin/grasp2-GenGraph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bin/grasp2-GenGraph.dir/depend

