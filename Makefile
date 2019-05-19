# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/cczhong/Codes/MANA/CMakeFiles /home/cczhong/Codes/MANA/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/cczhong/Codes/MANA/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named bin/grasp2-assemble

# Build rule for target.
bin/grasp2-assemble: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 bin/grasp2-assemble
.PHONY : bin/grasp2-assemble

# fast build rule for target.
bin/grasp2-assemble/fast:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/build
.PHONY : bin/grasp2-assemble/fast

#=============================================================================
# Target rules for targets named bin/grasp2-build

# Build rule for target.
bin/grasp2-build: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 bin/grasp2-build
.PHONY : bin/grasp2-build

# fast build rule for target.
bin/grasp2-build/fast:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/build
.PHONY : bin/grasp2-build/fast

GenGraph/align_batch.o: GenGraph/align_batch.cc.o

.PHONY : GenGraph/align_batch.o

# target to build an object file
GenGraph/align_batch.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/align_batch.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/align_batch.cc.o
.PHONY : GenGraph/align_batch.cc.o

GenGraph/align_batch.i: GenGraph/align_batch.cc.i

.PHONY : GenGraph/align_batch.i

# target to preprocess a source file
GenGraph/align_batch.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/align_batch.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/align_batch.cc.i
.PHONY : GenGraph/align_batch.cc.i

GenGraph/align_batch.s: GenGraph/align_batch.cc.s

.PHONY : GenGraph/align_batch.s

# target to generate assembly for a file
GenGraph/align_batch.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/align_batch.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/align_batch.cc.s
.PHONY : GenGraph/align_batch.cc.s

GenGraph/bio_alphabet.o: GenGraph/bio_alphabet.cc.o

.PHONY : GenGraph/bio_alphabet.o

# target to build an object file
GenGraph/bio_alphabet.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bio_alphabet.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bio_alphabet.cc.o
.PHONY : GenGraph/bio_alphabet.cc.o

GenGraph/bio_alphabet.i: GenGraph/bio_alphabet.cc.i

.PHONY : GenGraph/bio_alphabet.i

# target to preprocess a source file
GenGraph/bio_alphabet.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bio_alphabet.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bio_alphabet.cc.i
.PHONY : GenGraph/bio_alphabet.cc.i

GenGraph/bio_alphabet.s: GenGraph/bio_alphabet.cc.s

.PHONY : GenGraph/bio_alphabet.s

# target to generate assembly for a file
GenGraph/bio_alphabet.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bio_alphabet.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bio_alphabet.cc.s
.PHONY : GenGraph/bio_alphabet.cc.s

GenGraph/bwt.o: GenGraph/bwt.cc.o

.PHONY : GenGraph/bwt.o

# target to build an object file
GenGraph/bwt.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bwt.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bwt.cc.o
.PHONY : GenGraph/bwt.cc.o

GenGraph/bwt.i: GenGraph/bwt.cc.i

.PHONY : GenGraph/bwt.i

# target to preprocess a source file
GenGraph/bwt.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bwt.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bwt.cc.i
.PHONY : GenGraph/bwt.cc.i

GenGraph/bwt.s: GenGraph/bwt.cc.s

.PHONY : GenGraph/bwt.s

# target to generate assembly for a file
GenGraph/bwt.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bwt.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bwt.cc.s
.PHONY : GenGraph/bwt.cc.s

GenGraph/bwt_search.o: GenGraph/bwt_search.cc.o

.PHONY : GenGraph/bwt_search.o

# target to build an object file
GenGraph/bwt_search.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bwt_search.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bwt_search.cc.o
.PHONY : GenGraph/bwt_search.cc.o

GenGraph/bwt_search.i: GenGraph/bwt_search.cc.i

.PHONY : GenGraph/bwt_search.i

# target to preprocess a source file
GenGraph/bwt_search.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bwt_search.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bwt_search.cc.i
.PHONY : GenGraph/bwt_search.cc.i

GenGraph/bwt_search.s: GenGraph/bwt_search.cc.s

.PHONY : GenGraph/bwt_search.s

# target to generate assembly for a file
GenGraph/bwt_search.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/bwt_search.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/bwt_search.cc.s
.PHONY : GenGraph/bwt_search.cc.s

GenGraph/concatenator.o: GenGraph/concatenator.cc.o

.PHONY : GenGraph/concatenator.o

# target to build an object file
GenGraph/concatenator.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/concatenator.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/concatenator.cc.o
.PHONY : GenGraph/concatenator.cc.o

GenGraph/concatenator.i: GenGraph/concatenator.cc.i

.PHONY : GenGraph/concatenator.i

# target to preprocess a source file
GenGraph/concatenator.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/concatenator.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/concatenator.cc.i
.PHONY : GenGraph/concatenator.cc.i

GenGraph/concatenator.s: GenGraph/concatenator.cc.s

.PHONY : GenGraph/concatenator.s

# target to generate assembly for a file
GenGraph/concatenator.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/concatenator.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/concatenator.cc.s
.PHONY : GenGraph/concatenator.cc.s

GenGraph/contig_refinement.o: GenGraph/contig_refinement.cc.o

.PHONY : GenGraph/contig_refinement.o

# target to build an object file
GenGraph/contig_refinement.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/contig_refinement.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/contig_refinement.cc.o
.PHONY : GenGraph/contig_refinement.cc.o

GenGraph/contig_refinement.i: GenGraph/contig_refinement.cc.i

.PHONY : GenGraph/contig_refinement.i

# target to preprocess a source file
GenGraph/contig_refinement.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/contig_refinement.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/contig_refinement.cc.i
.PHONY : GenGraph/contig_refinement.cc.i

GenGraph/contig_refinement.s: GenGraph/contig_refinement.cc.s

.PHONY : GenGraph/contig_refinement.s

# target to generate assembly for a file
GenGraph/contig_refinement.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/contig_refinement.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/contig_refinement.cc.s
.PHONY : GenGraph/contig_refinement.cc.s

GenGraph/database_index.o: GenGraph/database_index.cc.o

.PHONY : GenGraph/database_index.o

# target to build an object file
GenGraph/database_index.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/database_index.cc.o
.PHONY : GenGraph/database_index.cc.o

GenGraph/database_index.i: GenGraph/database_index.cc.i

.PHONY : GenGraph/database_index.i

# target to preprocess a source file
GenGraph/database_index.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/database_index.cc.i
.PHONY : GenGraph/database_index.cc.i

GenGraph/database_index.s: GenGraph/database_index.cc.s

.PHONY : GenGraph/database_index.s

# target to generate assembly for a file
GenGraph/database_index.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/database_index.cc.s
.PHONY : GenGraph/database_index.cc.s

GenGraph/divsufsort.o: GenGraph/divsufsort.cc.o

.PHONY : GenGraph/divsufsort.o

# target to build an object file
GenGraph/divsufsort.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/divsufsort.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/divsufsort.cc.o
.PHONY : GenGraph/divsufsort.cc.o

GenGraph/divsufsort.i: GenGraph/divsufsort.cc.i

.PHONY : GenGraph/divsufsort.i

# target to preprocess a source file
GenGraph/divsufsort.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/divsufsort.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/divsufsort.cc.i
.PHONY : GenGraph/divsufsort.cc.i

GenGraph/divsufsort.s: GenGraph/divsufsort.cc.s

.PHONY : GenGraph/divsufsort.s

# target to generate assembly for a file
GenGraph/divsufsort.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/divsufsort.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/divsufsort.cc.s
.PHONY : GenGraph/divsufsort.cc.s

GenGraph/gsa.o: GenGraph/gsa.cc.o

.PHONY : GenGraph/gsa.o

# target to build an object file
GenGraph/gsa.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/gsa.cc.o
.PHONY : GenGraph/gsa.cc.o

GenGraph/gsa.i: GenGraph/gsa.cc.i

.PHONY : GenGraph/gsa.i

# target to preprocess a source file
GenGraph/gsa.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/gsa.cc.i
.PHONY : GenGraph/gsa.cc.i

GenGraph/gsa.s: GenGraph/gsa.cc.s

.PHONY : GenGraph/gsa.s

# target to generate assembly for a file
GenGraph/gsa.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/gsa.cc.s
.PHONY : GenGraph/gsa.cc.s

GenGraph/interval_array.o: GenGraph/interval_array.cc.o

.PHONY : GenGraph/interval_array.o

# target to build an object file
GenGraph/interval_array.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/interval_array.cc.o
.PHONY : GenGraph/interval_array.cc.o

GenGraph/interval_array.i: GenGraph/interval_array.cc.i

.PHONY : GenGraph/interval_array.i

# target to preprocess a source file
GenGraph/interval_array.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/interval_array.cc.i
.PHONY : GenGraph/interval_array.cc.i

GenGraph/interval_array.s: GenGraph/interval_array.cc.s

.PHONY : GenGraph/interval_array.s

# target to generate assembly for a file
GenGraph/interval_array.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/interval_array.cc.s
.PHONY : GenGraph/interval_array.cc.s

GenGraph/kmer_filtering.o: GenGraph/kmer_filtering.cc.o

.PHONY : GenGraph/kmer_filtering.o

# target to build an object file
GenGraph/kmer_filtering.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/kmer_filtering.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/kmer_filtering.cc.o
.PHONY : GenGraph/kmer_filtering.cc.o

GenGraph/kmer_filtering.i: GenGraph/kmer_filtering.cc.i

.PHONY : GenGraph/kmer_filtering.i

# target to preprocess a source file
GenGraph/kmer_filtering.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/kmer_filtering.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/kmer_filtering.cc.i
.PHONY : GenGraph/kmer_filtering.cc.i

GenGraph/kmer_filtering.s: GenGraph/kmer_filtering.cc.s

.PHONY : GenGraph/kmer_filtering.s

# target to generate assembly for a file
GenGraph/kmer_filtering.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/kmer_filtering.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/kmer_filtering.cc.s
.PHONY : GenGraph/kmer_filtering.cc.s

GenGraph/kmer_unitcoder.o: GenGraph/kmer_unitcoder.cc.o

.PHONY : GenGraph/kmer_unitcoder.o

# target to build an object file
GenGraph/kmer_unitcoder.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/kmer_unitcoder.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/kmer_unitcoder.cc.o
.PHONY : GenGraph/kmer_unitcoder.cc.o

GenGraph/kmer_unitcoder.i: GenGraph/kmer_unitcoder.cc.i

.PHONY : GenGraph/kmer_unitcoder.i

# target to preprocess a source file
GenGraph/kmer_unitcoder.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/kmer_unitcoder.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/kmer_unitcoder.cc.i
.PHONY : GenGraph/kmer_unitcoder.cc.i

GenGraph/kmer_unitcoder.s: GenGraph/kmer_unitcoder.cc.s

.PHONY : GenGraph/kmer_unitcoder.s

# target to generate assembly for a file
GenGraph/kmer_unitcoder.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/kmer_unitcoder.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/kmer_unitcoder.cc.s
.PHONY : GenGraph/kmer_unitcoder.cc.s

GenGraph/loader.o: GenGraph/loader.cc.o

.PHONY : GenGraph/loader.o

# target to build an object file
GenGraph/loader.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/loader.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/loader.cc.o
.PHONY : GenGraph/loader.cc.o

GenGraph/loader.i: GenGraph/loader.cc.i

.PHONY : GenGraph/loader.i

# target to preprocess a source file
GenGraph/loader.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/loader.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/loader.cc.i
.PHONY : GenGraph/loader.cc.i

GenGraph/loader.s: GenGraph/loader.cc.s

.PHONY : GenGraph/loader.s

# target to generate assembly for a file
GenGraph/loader.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/loader.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/loader.cc.s
.PHONY : GenGraph/loader.cc.s

GenGraph/main_assemble.o: GenGraph/main_assemble.cc.o

.PHONY : GenGraph/main_assemble.o

# target to build an object file
GenGraph/main_assemble.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/main_assemble.cc.o
.PHONY : GenGraph/main_assemble.cc.o

GenGraph/main_assemble.i: GenGraph/main_assemble.cc.i

.PHONY : GenGraph/main_assemble.i

# target to preprocess a source file
GenGraph/main_assemble.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/main_assemble.cc.i
.PHONY : GenGraph/main_assemble.cc.i

GenGraph/main_assemble.s: GenGraph/main_assemble.cc.s

.PHONY : GenGraph/main_assemble.s

# target to generate assembly for a file
GenGraph/main_assemble.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/main_assemble.cc.s
.PHONY : GenGraph/main_assemble.cc.s

GenGraph/main_build.o: GenGraph/main_build.cc.o

.PHONY : GenGraph/main_build.o

# target to build an object file
GenGraph/main_build.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/main_build.cc.o
.PHONY : GenGraph/main_build.cc.o

GenGraph/main_build.i: GenGraph/main_build.cc.i

.PHONY : GenGraph/main_build.i

# target to preprocess a source file
GenGraph/main_build.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/main_build.cc.i
.PHONY : GenGraph/main_build.cc.i

GenGraph/main_build.s: GenGraph/main_build.cc.s

.PHONY : GenGraph/main_build.s

# target to generate assembly for a file
GenGraph/main_build.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/main_build.cc.s
.PHONY : GenGraph/main_build.cc.s

GenGraph/minimizer_sort.o: GenGraph/minimizer_sort.cc.o

.PHONY : GenGraph/minimizer_sort.o

# target to build an object file
GenGraph/minimizer_sort.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/minimizer_sort.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/minimizer_sort.cc.o
.PHONY : GenGraph/minimizer_sort.cc.o

GenGraph/minimizer_sort.i: GenGraph/minimizer_sort.cc.i

.PHONY : GenGraph/minimizer_sort.i

# target to preprocess a source file
GenGraph/minimizer_sort.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/minimizer_sort.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/minimizer_sort.cc.i
.PHONY : GenGraph/minimizer_sort.cc.i

GenGraph/minimizer_sort.s: GenGraph/minimizer_sort.cc.s

.PHONY : GenGraph/minimizer_sort.s

# target to generate assembly for a file
GenGraph/minimizer_sort.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/minimizer_sort.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/minimizer_sort.cc.s
.PHONY : GenGraph/minimizer_sort.cc.s

GenGraph/sequence_search.o: GenGraph/sequence_search.cc.o

.PHONY : GenGraph/sequence_search.o

# target to build an object file
GenGraph/sequence_search.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/sequence_search.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sequence_search.cc.o
.PHONY : GenGraph/sequence_search.cc.o

GenGraph/sequence_search.i: GenGraph/sequence_search.cc.i

.PHONY : GenGraph/sequence_search.i

# target to preprocess a source file
GenGraph/sequence_search.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/sequence_search.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sequence_search.cc.i
.PHONY : GenGraph/sequence_search.cc.i

GenGraph/sequence_search.s: GenGraph/sequence_search.cc.s

.PHONY : GenGraph/sequence_search.s

# target to generate assembly for a file
GenGraph/sequence_search.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/sequence_search.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sequence_search.cc.s
.PHONY : GenGraph/sequence_search.cc.s

GenGraph/sfa.o: GenGraph/sfa.cc.o

.PHONY : GenGraph/sfa.o

# target to build an object file
GenGraph/sfa.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sfa.cc.o
.PHONY : GenGraph/sfa.cc.o

GenGraph/sfa.i: GenGraph/sfa.cc.i

.PHONY : GenGraph/sfa.i

# target to preprocess a source file
GenGraph/sfa.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sfa.cc.i
.PHONY : GenGraph/sfa.cc.i

GenGraph/sfa.s: GenGraph/sfa.cc.s

.PHONY : GenGraph/sfa.s

# target to generate assembly for a file
GenGraph/sfa.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sfa.cc.s
.PHONY : GenGraph/sfa.cc.s

GenGraph/sfa_build.o: GenGraph/sfa_build.cc.o

.PHONY : GenGraph/sfa_build.o

# target to build an object file
GenGraph/sfa_build.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sfa_build.cc.o
.PHONY : GenGraph/sfa_build.cc.o

GenGraph/sfa_build.i: GenGraph/sfa_build.cc.i

.PHONY : GenGraph/sfa_build.i

# target to preprocess a source file
GenGraph/sfa_build.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sfa_build.cc.i
.PHONY : GenGraph/sfa_build.cc.i

GenGraph/sfa_build.s: GenGraph/sfa_build.cc.s

.PHONY : GenGraph/sfa_build.s

# target to generate assembly for a file
GenGraph/sfa_build.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sfa_build.cc.s
.PHONY : GenGraph/sfa_build.cc.s

GenGraph/sssort.o: GenGraph/sssort.cc.o

.PHONY : GenGraph/sssort.o

# target to build an object file
GenGraph/sssort.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/sssort.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sssort.cc.o
.PHONY : GenGraph/sssort.cc.o

GenGraph/sssort.i: GenGraph/sssort.cc.i

.PHONY : GenGraph/sssort.i

# target to preprocess a source file
GenGraph/sssort.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/sssort.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sssort.cc.i
.PHONY : GenGraph/sssort.cc.i

GenGraph/sssort.s: GenGraph/sssort.cc.s

.PHONY : GenGraph/sssort.s

# target to generate assembly for a file
GenGraph/sssort.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/sssort.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/sssort.cc.s
.PHONY : GenGraph/sssort.cc.s

GenGraph/string_graph.o: GenGraph/string_graph.cc.o

.PHONY : GenGraph/string_graph.o

# target to build an object file
GenGraph/string_graph.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/string_graph.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/string_graph.cc.o
.PHONY : GenGraph/string_graph.cc.o

GenGraph/string_graph.i: GenGraph/string_graph.cc.i

.PHONY : GenGraph/string_graph.i

# target to preprocess a source file
GenGraph/string_graph.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/string_graph.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/string_graph.cc.i
.PHONY : GenGraph/string_graph.cc.i

GenGraph/string_graph.s: GenGraph/string_graph.cc.s

.PHONY : GenGraph/string_graph.s

# target to generate assembly for a file
GenGraph/string_graph.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/string_graph.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/string_graph.cc.s
.PHONY : GenGraph/string_graph.cc.s

GenGraph/trsort.o: GenGraph/trsort.cc.o

.PHONY : GenGraph/trsort.o

# target to build an object file
GenGraph/trsort.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/trsort.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/trsort.cc.o
.PHONY : GenGraph/trsort.cc.o

GenGraph/trsort.i: GenGraph/trsort.cc.i

.PHONY : GenGraph/trsort.i

# target to preprocess a source file
GenGraph/trsort.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/trsort.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/trsort.cc.i
.PHONY : GenGraph/trsort.cc.i

GenGraph/trsort.s: GenGraph/trsort.cc.s

.PHONY : GenGraph/trsort.s

# target to generate assembly for a file
GenGraph/trsort.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/trsort.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/trsort.cc.s
.PHONY : GenGraph/trsort.cc.s

GenGraph/utils.o: GenGraph/utils.cc.o

.PHONY : GenGraph/utils.o

# target to build an object file
GenGraph/utils.cc.o:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/utils.cc.o
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/utils.cc.o
.PHONY : GenGraph/utils.cc.o

GenGraph/utils.i: GenGraph/utils.cc.i

.PHONY : GenGraph/utils.i

# target to preprocess a source file
GenGraph/utils.cc.i:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/utils.cc.i
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/utils.cc.i
.PHONY : GenGraph/utils.cc.i

GenGraph/utils.s: GenGraph/utils.cc.s

.PHONY : GenGraph/utils.s

# target to generate assembly for a file
GenGraph/utils.cc.s:
	$(MAKE) -f CMakeFiles/bin/grasp2-assemble.dir/build.make CMakeFiles/bin/grasp2-assemble.dir/GenGraph/utils.cc.s
	$(MAKE) -f CMakeFiles/bin/grasp2-build.dir/build.make CMakeFiles/bin/grasp2-build.dir/GenGraph/utils.cc.s
.PHONY : GenGraph/utils.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... bin/grasp2-assemble"
	@echo "... bin/grasp2-build"
	@echo "... GenGraph/align_batch.o"
	@echo "... GenGraph/align_batch.i"
	@echo "... GenGraph/align_batch.s"
	@echo "... GenGraph/bio_alphabet.o"
	@echo "... GenGraph/bio_alphabet.i"
	@echo "... GenGraph/bio_alphabet.s"
	@echo "... GenGraph/bwt.o"
	@echo "... GenGraph/bwt.i"
	@echo "... GenGraph/bwt.s"
	@echo "... GenGraph/bwt_search.o"
	@echo "... GenGraph/bwt_search.i"
	@echo "... GenGraph/bwt_search.s"
	@echo "... GenGraph/concatenator.o"
	@echo "... GenGraph/concatenator.i"
	@echo "... GenGraph/concatenator.s"
	@echo "... GenGraph/contig_refinement.o"
	@echo "... GenGraph/contig_refinement.i"
	@echo "... GenGraph/contig_refinement.s"
	@echo "... GenGraph/database_index.o"
	@echo "... GenGraph/database_index.i"
	@echo "... GenGraph/database_index.s"
	@echo "... GenGraph/divsufsort.o"
	@echo "... GenGraph/divsufsort.i"
	@echo "... GenGraph/divsufsort.s"
	@echo "... GenGraph/gsa.o"
	@echo "... GenGraph/gsa.i"
	@echo "... GenGraph/gsa.s"
	@echo "... GenGraph/interval_array.o"
	@echo "... GenGraph/interval_array.i"
	@echo "... GenGraph/interval_array.s"
	@echo "... GenGraph/kmer_filtering.o"
	@echo "... GenGraph/kmer_filtering.i"
	@echo "... GenGraph/kmer_filtering.s"
	@echo "... GenGraph/kmer_unitcoder.o"
	@echo "... GenGraph/kmer_unitcoder.i"
	@echo "... GenGraph/kmer_unitcoder.s"
	@echo "... GenGraph/loader.o"
	@echo "... GenGraph/loader.i"
	@echo "... GenGraph/loader.s"
	@echo "... GenGraph/main_assemble.o"
	@echo "... GenGraph/main_assemble.i"
	@echo "... GenGraph/main_assemble.s"
	@echo "... GenGraph/main_build.o"
	@echo "... GenGraph/main_build.i"
	@echo "... GenGraph/main_build.s"
	@echo "... GenGraph/minimizer_sort.o"
	@echo "... GenGraph/minimizer_sort.i"
	@echo "... GenGraph/minimizer_sort.s"
	@echo "... GenGraph/sequence_search.o"
	@echo "... GenGraph/sequence_search.i"
	@echo "... GenGraph/sequence_search.s"
	@echo "... GenGraph/sfa.o"
	@echo "... GenGraph/sfa.i"
	@echo "... GenGraph/sfa.s"
	@echo "... GenGraph/sfa_build.o"
	@echo "... GenGraph/sfa_build.i"
	@echo "... GenGraph/sfa_build.s"
	@echo "... GenGraph/sssort.o"
	@echo "... GenGraph/sssort.i"
	@echo "... GenGraph/sssort.s"
	@echo "... GenGraph/string_graph.o"
	@echo "... GenGraph/string_graph.i"
	@echo "... GenGraph/string_graph.s"
	@echo "... GenGraph/trsort.o"
	@echo "... GenGraph/trsort.i"
	@echo "... GenGraph/trsort.s"
	@echo "... GenGraph/utils.o"
	@echo "... GenGraph/utils.i"
	@echo "... GenGraph/utils.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
