#ifndef _SFA_BUILD_
#define _SFA_BUILD_

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

#define MAX_BLOCK 32

class SFABuild {
 public:
  SFABuild(void);
  SFABuild(std::string& seq_file);
  SFABuild(std::string& seq_file, std::string& block_file);
  SFABuild(SFABuild &seq_obj);
  ~SFABuild(void);
  void DestructSFA(void);
  void DestructSequences(void);
  void LoadSequences(std::string& seq_file);
  void LoadBlockSize(std::string& block_file);
  void CopySequences(int num_sequences, char **source, char **target);
  void InPlaceReverse(void);
  void BuildSFADefault(void);
  void BuildSFAMulti(const SFAIDXTYPE& max_size, std::string &dir, std::string &file_stem);
  void BuildKeyArrayDefault(void);
  void DumpSFA(std::string& dir, std::string& file_stem);
  void LoadSFA(std::string& dir, std::string& file_stem);
  void CountDBSize(void);
  double GetDBSizeInMegaBase(void);
  void PrintAllSeqs(void);
  std::pair<SFAIDXTYPE, SFAIDXTYPE> SearchSFA(std::string& search_seed);
  void GetMaxExtInfoWithinRange(
      std::pair<SFAIDXTYPE, SFAIDXTYPE>& range, 
      std::list<GSATYPE>& pos_list
  );
  std::string GetSequence(const RIDTYPE& index);
  std::string GetSequence(const int& block_ID, const RIDTYPE& index);

  std::string GetHeader(int index);
  std::string GetSuffixSFA(int index);
  
  std::string GetSuffixSeq(
    const GSATYPE& s,   // the location of the sequence needs to be copied
    const int& l        // the length of the sequence needs to be copied
  );

  std::string GetSuffixSeq(
    const int& block_ID,  // the ID of the block if multiple SFA was generated
    const GSATYPE& s,     // the location of the sequence needs to be copied
    const int& l          // the length of the sequence needs to be copied
  );

  bool CheckMultiParam(const SFAIDXTYPE& max_size);

  int GetNumBlocks(void); // access the number of index blocks

  int GetSeqLen(const RIDTYPE& r) {
    assert(is_sequence_loaded_);
    assert(r < num_seqs_);
    return seq_len_[r];
  }

  // returns the length of the suffix
  int GetSufLen(const int& block_ID, const GSATYPE& s); 

  RIDTYPE GetFullRID(const int& block_ID, const RIDTYPE& r);
  
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
  int* seq_len_;
  GSA* suffix_array_;
  std::vector<SFAIDXTYPE> key_array_;
  bool is_header_loaded_, is_sequence_loaded_, is_sfa_built_, is_k_array_built_;
  bool is_size_counted_;
  bool is_multi_;
  int num_blocks_;
  std::vector<RIDTYPE> block_size_; 
  // methods begin
  void BuildSuffixArray(char** sequences, GSA* suffix_array);
  void BuildKeyArray(GSA* suffix_array, std::vector<SFAIDXTYPE>& key_array);
};

#endif
