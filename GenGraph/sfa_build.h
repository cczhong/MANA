#include "sequence.h"
#include "gsa.h"

#include <boost/filesystem.hpp>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <string>
#include <list>

#ifndef _SFA_BUILD_
#define _SFA_BUILD_

struct PositionType {
  RIDTYPE rid;
  LCPTYPE pos;
};

class SFABuild {
 public:
  SFABuild(void);
  SFABuild(std::string &file_name);
  SFABuild(SFABuild &seq_obj);
  ~SFABuild(void);
  void DestructSFA(void);
  void DestructSequences(void);
  void LoadSequences(std::string& file_name);
  void CopySequences(int num_sequences, char **source, char **target);
  void InPlaceReverse(void);
  void BuildSFADefault(void);
  void BuildKeyArrayDefault(void);
  void DumpSFA(std::string& dir, std::string& file_stem);
  void LoadSFA(std::string& dir, std::string& file_stem);
  void CountDBSize(void);
  double GetDBSizeInMegaBase(void);
  void PrintAllSeqs(void);
  std::pair<SFAIDXTYPE, SFAIDXTYPE> SearchSFA(std::string& search_seed);
  void GetMaxExtInfoWithinRange(
      std::pair<SFAIDXTYPE, SFAIDXTYPE>& range, 
      std::list<PositionType>& pos_list
  );
  std::string GetSequence(int index);
  std::string GetHeader(int index);
  std::string GetSuffixSFA(int index);
  std::string GetFileStem(const std::string& path);
  friend class DatabaseIndex;
  friend class ReachableReads;
  friend class ReadAlignment;
  friend class GreedyAssembly;
  friend class AssembleExtend;
  friend class ReMap;
  friend class Unitiger;
 protected:
  // data begin
  int num_seqs_;
  double db_size_MB_;
  char** header_;
  char** sequence_;
  GSA* suffix_array_;
  std::vector<SFAIDXTYPE> key_array_;
  bool is_header_loaded_, is_sequence_loaded_, is_sfa_built_, is_k_array_built_;
  bool is_size_counted_;
  // methods begin
  void BuildSuffixArray(char** sequences, GSA* suffix_array);
  void BuildKeyArray(GSA* suffix_array, std::vector<SFAIDXTYPE>& key_array);
};

#endif
