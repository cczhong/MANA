#include "gsa.h"
#include "sfa_build.h"
#include "clump.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <sstream>
#include <limits>

#ifndef _OVERLAP_H_
#define _OVERLAP_H_

class Overlap   {
 public:
  Overlap(void);
  ~Overlap(void);

  void DetectOverlaps(
    SFABuild& seqs,                 // the original set of sequences
    const std::string& dir,         // the folder that contains the index files
    const std::string& file_stem,   // the file stem
    const int& min_overlap          // the minimum overlap length
  );
  
 private:
  
  void CountNumIndex(
    const std::string& dir,             // the folder that contains the index files
    const std::string& file_stem,       // the file stem to load 
    int &num_index                      // (output) the number of index pieces
  );

  void OpenIndexFiles(
    const std::string& dir, 
    const std::string& file_stem, 
    const int num_index,                    // the number of indexes to load
    std::ifstream* GSAfh,  // (output) the array that contains the file handles to the generalized suffix array indexes
    std::ifstream* LCPfh   // (output) the array that contains the file handles to the LCPs
  );

  bool isFileExists (const std::string& file_name) {
    std::ifstream f(file_name.c_str());
    return f.good();
  }

  // find the perfect suffix-prefix overlap from the set of clumps sharing the same k-mer
  void ResolvePerfectOverlap(
    std::vector<CLUMPTYPE>& m_clump,
    SFABuild& seqs
  );

};

#endif