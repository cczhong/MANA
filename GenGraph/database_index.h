#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <tuple>
#include <string>
#include <list>

#include "sfa_build.h"
#include "reduced_alphabet.h"
#include "scoring_prot.h"

#ifndef _DATABASE_INDEX_H_
#define _DATABASE_INDEX_H_

struct OVERLAPTYPE  {
  RIDTYPE doc;  // the target document ID
  POSTYPE len;  // the length of the overlap
};

struct READPAIRTYPE {
  RIDTYPE doc_fw;
  RIDTYPE doc_re;
  bool init_fw;
  bool init_re;
  int q_pos_fw, r_pos_fw; // the position of the match in the query and the read, respectively
  int q_pos_re, r_pos_re; // the position of the match in the query and the read, respectively
  int overlap;
  int score;  // matching score for the seeds in the query and target
  double aln_e_value;
};

class DatabaseIndex {
 public:
  DatabaseIndex(int in_alph_id, int in_seed_len, int in_overlap_len);
  DatabaseIndex();
  ~DatabaseIndex();
  void BuildSeedmerMap(SFABuild& seq_obj); // map from seed-mer to sequence position
  // the key is the reduced alphabet string, the value is
  // a list of RID (read ID) - POS (position) pairs of the 
  // corresponding strings in the original alphabet
  
  /*
  void CreateReducedMap(
      std::unordered_map<std::string, std::list<std::string> >& reduc_alph_map
  );
  void DumpReducedMap(
      std::unordered_map<std::string, std::list<std::string> >& reduc_alph_map,
      std::string& out_file
  );
  void LoadReducedMap(
      std::string& in_file, 
      std::unordered_map<std::string, std::unordered_map<std::string, bool> >& reduc_alph_map
  );
  */
  
  // the key is a pair of RID-POS that specify a given seed k-mer
  // the value is a list of RIDs which can be used to extend the seed
  void CreateSeedExt(
      int min_seed_coverage, SFABuild& seq_obj, SFABuild& rev_seq_obj, 
      std::unordered_map<std::string, std::list<READPAIRTYPE> >& ext_seed
  );
  void DumpSeedExt(
      std::unordered_map<std::string, std::list<READPAIRTYPE> >& ext_seed,
      std::string& out_file
  );
  void LoadSeedExt(
      std::string& in_file,
      std::unordered_map<std::string, std::list<READPAIRTYPE> >& ext_seed
  );
  // the key is a RID, and the value is a list of RIDs that have
  // significant overlap with the key read
  void CreateReadExt(
      int min_ext_coverage,
      SFABuild& seq_obj, 
      SFABuild& rev_seq_obj,
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& fw_ext_read,
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& re_ext_read
  );
  void DumpReadExt(
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& fw_ext_read,
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& re_ext_read,
      std::string& out_file
  );
  void LoadReadExt(
      std::string& in_file,
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& fw_ext_read,
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& re_ext_read
  );
  /*
  void CreateHighScoreMer(
      std::unordered_map<std::string, std::list<std::string> >& high_score_match
  );
  void DumpHighScoreMer(
      std::unordered_map<std::string, std::list<std::string> >& high_score_match,
      std::string& out_file
  );
  void LoadHighScoreMer(
      std::string& in_file,
      std::unordered_map<std::string, std::list<std::string> >& high_score_match
  );
  */

  int GetSeedLen(void);
  int GetOverlapLen(void);
  int GetAlphID(void);
  void GetSeedExt(
      SFABuild &seq_obj, int min_seed_coverage,
      std::unordered_map<std::string, std::list<GSATYPE> > &ext
  );
  void GetSeedExtRev(
      SFABuild &rev_seq_obj, int min_seed_coverage,
      std::unordered_map<std::string, std::list<GSATYPE> > &rev_ext
  );
  void MatchSeedPair(
      SFABuild &seq_obj,
      std::unordered_map<std::string, std::list<GSATYPE> > &ext,
      std::unordered_map<std::string, std::list<GSATYPE> > &rev_ext,
      std::unordered_map<std::string, std::list<READPAIRTYPE> >& ext_pair
  );
  void CreateReadExtWorker(
      int min_ext_coverage, SFABuild& seq_obj, 
      std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> >& ext_read
  );
 protected:
  bool is_smm_built_; // check if the seed_mer_map is built
  bool is_alph_set_;
  bool is_seed_len_set_;
  bool is_overlap_len_set_;
  int alph_id_;
  int seed_len_;
  int overlap_len_;
  std::unordered_map<std::string, GSATYPE> seed_mer_map_;
  bool IsExtRedundant(
      std::vector<std::list<GSATYPE> >& fw_ext, 
      std::vector<std::list<GSATYPE> >& re_ext, 
      std::unordered_map<RIDTYPE, std::set<int> >& read_map_table, 
      std::list<GSATYPE>& fw_phase, std::list<GSATYPE>& re_phase, 
      int& map_ID
  );
  bool IsPositionListMatch(
      std::list<GSATYPE>& l1, std::list<GSATYPE>& l2
  );
  bool IsSeqCompatible(
      SFABuild& seq_obj, int seed_len,
      RIDTYPE fw_rid, POSTYPE fw_pos,
      RIDTYPE re_rid, POSTYPE re_pos
  );
  void MatchSeedPairSingle(
      SFABuild& seq_obj,
      std::list<GSATYPE>& fw_sp, std::list<GSATYPE>& re_sp,
      std::list<READPAIRTYPE>& read_pair
  );
};

#endif
