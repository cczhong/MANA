#include "gsa.h"
#include "sfa_build.h"
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
#include <functional>
#include <queue>

#ifndef _CLUMP_H_
#define _CLUMP_H_

struct CLUMPTYPE    {
    std::string key_prefix;             // the common prefix of the suffix block
    std::vector<GSATYPE> suf_block;     // the suffix block
    std::vector<LCPTYPE> lcp_block;     // the lcp block
    int index_ID;                       // the index where this clump was retrieved

    void operator= (const CLUMPTYPE& a)
    {
        this->key_prefix = a.key_prefix;
        this->suf_block = a.suf_block;
        this->lcp_block = a.lcp_block;
        this->index_ID = a.index_ID;
        return;
    }

};

class ClumpCmp  {
public:
    ClumpCmp()  {}
    bool operator() (const CLUMPTYPE &a, const CLUMPTYPE &b) const  
    {
        assert(a.key_prefix.length() == b.key_prefix.length());
        return (a.key_prefix > b.key_prefix);
    }
};


class Clump {

public:
    Clump(void) {   is_empty_ = true; is_init_ = false;   }
    ~Clump(void)    
    {   
        return;
    }

    // insert a clump into the priority queue
    void InsertClump(const CLUMPTYPE &a);

    // pop a clump from the priority queue
    CLUMPTYPE PopClump(void);

    // check if the priority queue is empty
    bool IsPQueueEmpty(void);

    // check if two clumps have the same key prefix
    bool IsClumpIdentical(const CLUMPTYPE& a, const CLUMPTYPE& b);

    // append clumps to "current_" (i.e. merge their suffix blocks)
    void AppendClump(const CLUMPTYPE& a);

    // initilize the clump priority queue from a vector of file handles
    void InitClumpPQueue(
        const int& num_index,               // the number of index files
        std::ifstream* GSAfh,               // the array that contains the file handles to the generalized suffix array indexes
        std::ifstream* LCPfh,               // the array that contains the file handles to the LCPs
        SFABuild &seqs,                     // the sequences
        const int& min_lcp                  // the minimum length of the LCP (or overlap) 
    );

    // retrieves the next clump from a given index 
    bool GetSuffixClump(
        std::ifstream* GSAfh,                // the array that contains the file handles to the generalized suffix array indexes
        std::ifstream* LCPfh,                // the array that contains the file handles to the LCPs
        const int& index_ID,                 // the ID of the index file
        SFABuild &seqs,                      // the sequences
        const int& min_lcp,                  // the minimum length of the LCP (or overlap)
        CLUMPTYPE &clump,                    // (output) the detected clump
        GSATYPE &buffer                      // (output) contains the first suffix of the next clump (if applicable)
    );  // returns true if clump recorded, returns false if reaches the end of file

    // get the full clump from multiple index blocks
    bool NextClump(
        std::ifstream* GSAfh,                // the array that contains the file handles to the generalized suffix array indexes
        std::ifstream* LCPfh,                // the array that contains the file handles to the LCPs
        SFABuild &seqs,                      // the sequences
        const int& min_lcp                   // the minimum length of the LCP (overlap)                 
    ); 

    // merge multiple clumps from different indexes
    void MergeMultiClumps(
        SFABuild& seqs                      // the sequences
    );

    // merge a pair of clumps
    CLUMPTYPE MergeTwoClumps(
        const CLUMPTYPE& a,
        const CLUMPTYPE& b,
        SFABuild& seqs
    );

    friend class Overlap;
    
private:
    std::priority_queue<CLUMPTYPE, std::vector<CLUMPTYPE>, ClumpCmp> clump_queue_;
    std::vector<CLUMPTYPE> multi_current_;  // the lexicographically smallest clump (could contain clumps from multiple index files)
    CLUMPTYPE current_;                     // the merged current block of suffixes
    bool is_empty_;                         // label whether "current_" contains information
    bool is_init_;                          // label whether the "is_EOF_" array is initialized
    std::vector<bool> is_EOF_;              // label whether the corresponding index has reached the end of file
    std::vector<CLUMPTYPE> blk_clumps_;     // the clumps retrieved for each index
    std::vector<GSATYPE> blk_buffers_;      // the buffer for each index
};

#endif