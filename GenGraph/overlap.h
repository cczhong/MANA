#include "gsa.h"
#include "sfa_build.h"
#include "clump.h"
#include "util_func.h"

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
#include <unordered_map>

#ifndef _OVERLAP_H_
#define _OVERLAP_H_

#define LARGE_U16 32000

// a group of reads sharing a common k-mer 
struct KREADGROUPTYPE  {
    std::vector<RIDTYPE> group_ID;
    std::vector<POSTYPE> kmer_begin;
    std::vector<POSTYPE> read_length;
};

struct OVERLAPTYPE  {
    RIDTYPE source;
    RIDTYPE target;
    LCPTYPE len;
};

struct PARTIALOLPTYPE   {
    RIDTYPE target;         // the target read ID (the source should be indicated by the index of the array)
    LCPTYPE len;            // the length of the merged sequene length
    bool is_transitive;     // indicate whether the edge is transitive (and therefore should be removed)
};

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

    void DetectUnitigs(
        const std::string& dir,                 // the index folder
        const std::string& file_stem,           // the file stem
        SFABuild& seqs 
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
        const CLUMPTYPE& clump,
        SFABuild& seqs,
        std::ofstream& OLPfh
    );

    // eiliminate the transitive edges and detect contained reads based on prefix overlap
    void InsertSource(
        const int& index_ID,
        const GSATYPE& suf,
        const LCPTYPE& lcp,
        const LCPTYPE& len,
        SFABuild& seqs,
        std::vector<GSATYPE>& suf_stack,
        std::vector<LCPTYPE>& lcp_stack,
        std::vector<LCPTYPE>& len_stack,
        int& s_index
    );

    // load information from the overlap index file
    void LoadOverlapInfo(
        const std::string& file,
        SFABuild& seqs,
        std::vector<std::vector<PARTIALOLPTYPE> >& olp_info
    );  

    void RemoveTransitiveEdges(
        std::vector<std::vector<PARTIALOLPTYPE> >& olp_info
    );

    bool is_overlap_init_;
    std::vector< std::unordered_map<RIDTYPE, OVERLAPTYPE> > overlap_fw_;  // source-indexed overlap information
    std::vector< std::unordered_map<RIDTYPE, OVERLAPTYPE> > overlap_re_;  // target-indexed overlap information

};

#endif