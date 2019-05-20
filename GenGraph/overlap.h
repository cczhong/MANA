#include "gsa.h"

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

#ifndef _OVERLAP_
#define _OVERLAP_

class Overlap   {
 public:
  Overlap(void);
  ~Overlap(void);

  void DetectOverlaps(
    const SFACHARTYPE** seq,    // the original set of sequences
    const std::string& dir,     // the folder that contains the index files
    const std::string& file_stem,   // the file stem
    const int num_errors        // the number of errors allowed during overlapping
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
    const int num_index,                // the number of indexes to load
    std::fstream* GSAfh,                // (output) the array that contains the file handles to the generalized suffix array indexes
    std::fstream* LCPfh                 // (output) the array that contains the file handles to the LCPs
  );

  bool isFileExists (const std::string& file_name) {
    std::ifstream f(file_name.c_str());
    return f.good();
  }
};

#endif