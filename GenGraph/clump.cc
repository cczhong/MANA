#include "clump.h"

using namespace std;

void Clump::InsertClump(
    const CLUMPTYPE &a
)   {
    clump_queue_.push(a);
    return;
}

CLUMPTYPE Clump::PopClump(void)
{
    assert(!clump_queue_.empty());
    CLUMPTYPE t = clump_queue_.top();
    clump_queue_.pop();
    return t;    
}

bool Clump::IsPQueueEmpty(void)
{
    return clump_queue_.empty();
}

bool Clump::IsClumpIdentical(const CLUMPTYPE& a, const CLUMPTYPE& b)    
{
    return (a.key_prefix == b.key_prefix);
}

void Clump::AppendClump(const CLUMPTYPE& a)   
{
    current_.push_back(a);
    is_empty_ = false;
    return;
}

void Clump::InitClumpPQueue(
    const int& num_index,                       // the number of index files
    std::ifstream* GSAfh,                       // the array that contains the file handles to the generalized suffix array indexes
    std::ifstream* LCPfh,                       // the array that contains the file handles to the LCPs
    SFABuild& seqs,                             // the sequences
    const int& min_lcp                          // the minimum length of the LCP (or overlap) 
)   {
    assert(!is_init_);
    assert(num_index > 0);
    // initilize the blocks
    is_EOF_.resize(num_index, false);
    blk_clumps_.resize(num_index);
    blk_buffers_.resize(num_index);
    for(int i = 0; i < num_index; ++ i) {
        blk_buffers_[i].doc = blk_buffers_[i].pos = 0;
    }
    // get the first clump for each index file
    for(int i = 0; i < num_index; ++ i)   {
        bool g = GetSuffixClump(
            GSAfh, LCPfh, i, seqs, min_lcp, blk_clumps_[i], blk_buffers_[i]
        );  // g indicates whether a clump is retrieved (hence not EOF)
        is_EOF_[i] = !g;
        // insert the clump into priority queue if clump detection is successful
        if(g)    {  InsertClump(blk_clumps_[i]);    }
    }
    return;
}

bool Clump::GetSuffixClump(
    std::ifstream* GSAfh,                // the array that contains the file handles to the generalized suffix array indexes
    std::ifstream* LCPfh,                // the array that contains the file handles to the LCPs
    const int& index_ID,                 // the ID of the index file
    SFABuild &seqs,                      // the sequences
    const int& min_lcp,                  // the minimum length of the LCP (or overlap)
    CLUMPTYPE &clump,                    // (output) the detected clump
    GSATYPE &buffer                      // (output) contains the first suffix of the next clump (if applicable)
)   {
    // check if reaching the end of file
    if(GSAfh[index_ID].eof() || LCPfh[index_ID].eof())    {
        return false;
    }
    // initialize the clump vector
    clump.index_ID = index_ID;
    clump.suf_block.resize(1000);
    clump.suf_block[0] = buffer;  // boundary case OK because the first LCP is always 0
    clump.key_prefix = seqs.GetSuffixSeq(index_ID, buffer, min_lcp);
    // load the clumps
    GSATYPE h;  // a temporary info holder
    LCPTYPE l;  // a temporary info holder
    int c_index = 0;
    bool is_c_open = false;

    // DEBUG
    //cout << "Printing block:    " << index_ID << endl;
    while(true) {
        GSAfh[index_ID].read((char*) &h.doc, sizeof(RIDTYPE));
        GSAfh[index_ID].read((char*) &h.pos, sizeof(POSTYPE));
        LCPfh[index_ID].read((char*) &l, sizeof(LCPTYPE));
        
        // DEBUG
        //string tmp = seqs.GetSuffixSeq(index_ID, h, min_lcp);
        //cout << "SF: (" << tmp << ")\tLCP:    " << l << endl;
        //delete [] tmp;

        // check LCP info to define clump
        if(l < min_lcp)    {
            if(is_c_open)    {
                // record the current suffix info into buffer
                buffer = h;
                clump.suf_block.resize(c_index);
                clump.key_prefix = seqs.GetSuffixSeq(index_ID, clump.suf_block[0], min_lcp);
                // terminate the current clump
                return true;
            }   else    {
                // update the head of the clump
                clump.suf_block[0] = h;
            }        
        }   else    {
            // indicates a continuation of the clump
            is_c_open = true;
            if(clump.suf_block.size() <= c_index - 1) {
                // double the clump size
                clump.suf_block.resize(2 * clump.suf_block.size());
            }
            clump.suf_block[++ c_index] = h;

        }
        if(GSAfh[index_ID].eof() || LCPfh[index_ID].eof())    {   break;  }
        
    }

    // DEBUG
    cout << "clump detected:    " << is_c_open << endl;

    if(is_c_open)   {   
        clump.suf_block.resize(c_index);  
        return true;
    }   else    {
        return false;
    }
}

bool Clump::NextClump(
    std::ifstream* GSAfh,                // the array that contains the file handles to the generalized suffix array indexes
    std::ifstream* LCPfh,                // the array that contains the file handles to the LCPs
    SFABuild &seqs,                      // the sequences
    const int& min_lcp                   // the minimum length of the LCP (overlap)
)   {
    if(clump_queue_.empty())    {  return false;  }
    // clear the "current_" data structure
    if(!is_empty_)  {   current_.clear();   }
    // DEBUG
    //cout << "Good here 0" << endl;
    // get the lexicographically smallest index among all indexes
    do
    {
        // DEBUG
        //cout << "seq:   " << clump_queue_.top().key_prefix << endl;
        //cout << "index  " << clump_queue_.top().index_ID << endl;

        AppendClump(clump_queue_.top());
        int idx = clump_queue_.top().index_ID;
        clump_queue_.pop();

        // DEBUG
        //cout << "Good here 1" << endl;

        // insert the new clump into queue
        if(!is_EOF_[idx])   {
            bool g = GetSuffixClump(
                GSAfh, LCPfh, idx, seqs, min_lcp, blk_clumps_[idx], blk_buffers_[idx]
            );  // g indicates whether a clump is retrieved (hence not EOF)

            // DEBUG
            //cout << "Good here 2" << endl;

            is_EOF_[idx] = !g;
            // insert the clump into priority queue if clump detection is successful
            if(g)    {  InsertClump(blk_clumps_[idx]);  }
            
            // DEBUG
            //cout << "Good here 3" << endl;
        }
    } while (!clump_queue_.empty() && current_[0].key_prefix == clump_queue_.top().key_prefix);
    
    // DEBUG
    //cout << "!!! NEW CLUMP" << endl;
    //for(int i = 0; i < current_.size(); ++ i)   {
    //    cout << "   seq:   " << current_[i].key_prefix << " ID: " << current_[i].index_ID << endl;
        //for(int j = 0; j < current_[i].suf_block.size(); ++ j)   {
        //    cout << "   " << current_[i].suf_block[j].doc << "   " << current_[i].suf_block[j].pos << endl;
        //}
    //}
    return true;
}