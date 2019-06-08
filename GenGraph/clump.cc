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
    multi_current_.push_back(a);
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
        blk_buffers_[i].pos = seqs.GetSufLen(i, blk_buffers_[i]);   // set the position to the last character
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
    // DEBUG
    //cout << "=========GetSuffix called:  " << index_ID << endl;
    // check if reaching the end of file
    if(GSAfh[index_ID].eof() || LCPfh[index_ID].eof())    {
        return false;
    }
    
    // load the clumps
    GSATYPE h;  // a temporary info holder
    LCPTYPE l;  // a temporary info holder
    int c_index = 0;
    // if the buffer suffix is sufficiently long, it contains a valid key prefix and is an open clump
    bool is_c_open = seqs.GetSufLen(index_ID, buffer) >= min_lcp ? true : false;
    if(is_c_open)   {
        // initialize the clump vector if buffer suffix is a part of the clump
        clump.index_ID = index_ID;
        clump.suf_block.resize(1000);
        clump.suf_block[0] = buffer;  // boundary case OK because the first LCP is always 0
        clump.lcp_block.resize(1000);
        clump.lcp_block[0] = 0;
        clump.key_prefix = seqs.GetSuffixSeq(index_ID, buffer, min_lcp);
        ++ c_index;
    }

    // DEBUG
    //cout << "Printing block:    " << index_ID << endl;
    //cout << "Initial block status:  " << is_c_open << endl;
    //cout << "buffer suffix length:  " << seqs.GetSufLen(index_ID, buffer) << endl;

    bool print_tag = false;
    while(true) {

        if(GSAfh[index_ID].peek() == EOF || LCPfh[index_ID].peek() == EOF)    {   break;  }

        GSAfh[index_ID].read((char*) &h.doc, sizeof(RIDTYPE));
        GSAfh[index_ID].read((char*) &h.pos, sizeof(POSTYPE));
        LCPfh[index_ID].read((char*) &l, sizeof(LCPTYPE));

        // DEBUG
        //RIDTYPE ck_id = seqs.GetFullRID(index_ID, h.doc);
        //if((ck_id == 100193 && h.pos == 0) || (ck_id == 59945 && h.pos == 1))    {
        //    print_tag = true;
        //    cout << "==========" << index_ID << "   " << h.doc << endl;
        //    for(int x = 0; x <= c_index; ++ x)   {
        //        string tmp = seqs.GetSuffixSeq(index_ID, clump.suf_block[x]);
        //        cout << "SF: (" << tmp << ")\tLCP:    " << l << endl;
        //    }
        //}
        //if(print_tag)    {
        //    string tmp = seqs.GetSuffixSeq(index_ID, h);
        //    cout << "SF: (" << tmp << ")\tLCP:    " << l << endl;
        //}

        // check whether the suffix is sufficiently long
        if(seqs.GetSufLen(index_ID, h) < min_lcp)    {
            if(is_c_open)    {
                // DEBUG
                //if(print_tag)    cout << "terminate the clump" << endl;  

                // record the current suffix info into buffer
                buffer = h;
                clump.suf_block.resize(c_index);
                clump.lcp_block.resize(c_index);
                return true;
            }   
            // otherwise do nothing (the suffix is not useful)      
        }   else    {
            // the suffix is sufficiently long
            if(is_c_open && l >= min_lcp)    {
                // is the current clump is open and the LCP is at least the minimum required overlap length
                // continue the clump
                // DEBUG
                //if(print_tag)   cout << "continue the clump" << endl;  
                if(clump.suf_block.size() <= c_index - 1) {
                    // double the clump size is necessary
                    clump.suf_block.resize(2 * clump.suf_block.size());
                    clump.lcp_block.resize(2 * clump.lcp_block.size());
                    assert(clump.suf_block.size() == clump.lcp_block.size());
                }
                clump.suf_block[c_index] = h;
                clump.lcp_block[c_index] = l;
                ++ c_index;
            }   else if (is_c_open && l < min_lcp)  {
                // DEBUG
                //if(print_tag)   cout << "terminate the clump" << endl;  
                // termination of the current clump
                buffer = h;
                clump.suf_block.resize(c_index);
                clump.lcp_block.resize(c_index);
                return true;
            }   else if (!is_c_open)  {
                // open the clump
                // DEBUG
                //if(print_tag)   cout << "open the clump" << endl;  
                is_c_open = true;
                clump.index_ID = index_ID;
                clump.suf_block.resize(1000);
                clump.suf_block[0] = h;
                clump.lcp_block.resize(1000);
                clump.lcp_block[0] = 0;     // this will be the first LCP and will never be used
                clump.key_prefix = seqs.GetSuffixSeq(index_ID, h, min_lcp);
                ++ c_index;
            }
        }
    }
    assert(!is_c_open || c_index >= 1);
    // DEBUG
    //if(print_tag)   cout << "clump status:    " << is_c_open << endl;
    //if(print_tag)   exit(0);

    if(is_c_open)   {   
        clump.suf_block.resize(c_index);  
        clump.lcp_block.resize(c_index);
        return true;
    }   else    {
        return false;
    }
}


/* OBSOLETE version
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
    clump.lcp_block.resize(1000);
    clump.lcp_block[0] = 0;
    clump.key_prefix = seqs.GetSuffixSeq(index_ID, buffer, min_lcp);
    // load the clumps
    GSATYPE h;  // a temporary info holder
    LCPTYPE l;  // a temporary info holder
    int c_index = 0;
    bool is_c_open = false;

    // DEBUG
    //cout << "Printing block:    " << index_ID << endl;
    while(true) {

        if(GSAfh[index_ID].peek() == EOF || LCPfh[index_ID].peek() == EOF)    {   break;  }

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
                clump.suf_block.resize(c_index + 1);
                clump.lcp_block.resize(c_index + 1);
                clump.key_prefix = seqs.GetSuffixSeq(index_ID, clump.suf_block[0], min_lcp);
                // terminate the current clump
                return true;
            }   else    {
                // update the head of the clump
                clump.suf_block[0] = h;
                clump.lcp_block[0] = l;
            }        
        }   else    {
            // indicates a continuation of the clump
            is_c_open = true;
            if(clump.suf_block.size() <= c_index - 1) {
                // double the clump size
                clump.suf_block.resize(2 * clump.suf_block.size());
                clump.lcp_block.resize(2 * clump.lcp_block.size());
                assert(clump.suf_block.size() == clump.lcp_block.size());
            }
            ++ c_index;
            clump.suf_block[c_index] = h;
            clump.lcp_block[c_index] = l;

        }
    }
    assert(!is_c_open || c_index >= 1);
    // DEBUG
    //cout << "clump detected:    " << is_c_open << endl;

    if(is_c_open)   {   
        clump.suf_block.resize(c_index + 1);  
        clump.lcp_block.resize(c_index + 1);
        return true;
    }   else    {
        return false;
    }
}
*/


bool Clump::NextClump(
    std::ifstream* GSAfh,                // the array that contains the file handles to the generalized suffix array indexes
    std::ifstream* LCPfh,                // the array that contains the file handles to the LCPs
    SFABuild &seqs,                      // the sequences
    const int& min_lcp                   // the minimum length of the LCP (overlap)
)   {
    if(clump_queue_.empty())    {  return false;  }
    // clear the "multi_current_" data structure
    if(!is_empty_)  {   multi_current_.clear();   }
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
        
        // DEBUG
        //cout << "The top of the queue:  "  << clump_queue_.top().key_prefix <<  "   " << clump_queue_.top().index_ID << endl;
        

        clump_queue_.pop();
        // insert the new clump into queue
        if(!is_EOF_[idx])   {
            bool g = GetSuffixClump(
                GSAfh, LCPfh, idx, seqs, min_lcp, blk_clumps_[idx], blk_buffers_[idx]
            );  // g indicates whether a clump is retrieved (hence not EOF)

            // DEBUG
            //cout << "New clump: " << endl;
            //for(int x = 0; x < blk_clumps_[idx].suf_block.size(); ++ x) {
                //cout << seqs.GetSuffixSeq(idx, blk_clumps_[idx].suf_block[x]) << endl;
            //    if(seqs.GetSuffixSeq(idx, blk_clumps_[idx].suf_block[x]).length() < min_lcp)    {
            //        cout << "Error!!!" << endl;
            //        cout << blk_clumps_[idx].suf_block[x].pos << endl;
            //        cout << seqs.GetSequence(idx, blk_clumps_[idx].suf_block[x].doc) << endl;
            //        cout << seqs.GetSufLen(idx, blk_clumps_[idx].suf_block[x]) << endl;
            //        cout << seqs.GetSuffixSeq(idx, blk_clumps_[idx].suf_block[x]) << endl;
            //        exit(0);
            //    }
            //}

            is_EOF_[idx] = !g;
            // insert the clump into priority queue if clump detection is successful
            if(g)    {  InsertClump(blk_clumps_[idx]);  }
            
            // DEBUG
            //cout << "Insert into queue: " << blk_clumps_[idx].key_prefix << "   " << idx << endl;
        }
    } while (!clump_queue_.empty() && multi_current_[0].key_prefix == clump_queue_.top().key_prefix);
    
    // DEBUG
    //if(multi_current_[0].key_prefix == "YHRYLAEFAS")    {
    //    cout << "!!! NEW CLUMP" << endl;
    //    for(int i = 0; i < multi_current_.size(); ++ i)   {
    //        cout << "   seq:   " << multi_current_[i].key_prefix << " block ID: " << multi_current_[i].index_ID << endl;
    //        for(int j = 0; j < multi_current_[i].suf_block.size(); ++ j)   {
    //            cout << "   " << multi_current_[i].suf_block[j].doc << "   " << multi_current_[i].suf_block[j].pos << endl;
    //            cout << seqs.GetSuffixSeq(multi_current_[i].index_ID, multi_current_[i].suf_block[j]) << endl;
    //        }
    //    }
    //}

    if(multi_current_.size() > 1)   {
        MergeMultiClumps(seqs);
    }   else if(multi_current_.size() == 1)    {
        current_ = multi_current_[0];
    }
    
    // DEBUG
    //cout << "!!! Merged CLUMP" << endl;
    //for(int j = 0; j < current_.suf_block.size(); ++ j)   {
    //    cout << "   " << current_.suf_block[j].doc << "   " << current_.suf_block[j].pos << endl;
    //    cout << seqs.GetSuffixSeq(current_.index_ID, current_.suf_block[j]) << endl;
    //}
    
    
    return true;
}

void Clump::MergeMultiClumps(
    SFABuild &seqs                      // the sequences
)   {
    // DEBUG
    //cout << "MergeMultiClumps:  " << multi_current_.size() << endl;
    for(int i = 0; i < ceil(log2(multi_current_.size())); ++ i)   {
        int k = (int) pow(2, i);
        for(int j = 0; j < multi_current_.size(); j += 2 * k)   {
            if(j + k < multi_current_.size())    {
                // merge the two clumps;
                // DEBUG
                //cout << "Merge pair clump:  " << j << ":" << multi_current_[j].index_ID << " " << j + k << ":" << multi_current_[j + k].index_ID << endl;
                multi_current_[j] = MergeTwoClumps(multi_current_[j], multi_current_[j + k], seqs);
            }
        }
    }
    // update the merged clump
    current_ = multi_current_[0];

    // DEBUG verify the merged clump is correct
    //for(int i = 0; i < current_.suf_block.size() - 1; ++ i)   {
    //    string str1 = seqs.GetSuffixSeq(current_.suf_block[i]);
    //    string str2 = seqs.GetSuffixSeq(current_.suf_block[i + 1]);
    //    if(str1.compare(str2) > 0 || str1.substr(0, current_.lcp_block[i + 1]) != str2.substr(0, current_.lcp_block[i + 1]))  {
    //        for(int i = 0; i < current_.suf_block.size(); ++ i)   {
    //            cout << seqs.GetSuffixSeq(current_.suf_block[i]) << "\t" << current_.lcp_block[i] << endl;
    //        }
    //        exit(0);
    //    }
    //}
    return;
}

CLUMPTYPE Clump::MergeTwoClumps(
    const CLUMPTYPE& a,
    const CLUMPTYPE& b,
    SFABuild& seqs
)   {
    assert(a.suf_block.size() == a.lcp_block.size());
    assert(b.suf_block.size() == b.lcp_block.size());
    assert(seqs.is_sequence_loaded_);
    assert(seqs.is_multi_);
    assert(a.index_ID < seqs.block_size_.size());
    assert(b.index_ID < seqs.block_size_.size());
    // boundary case
    if(b.suf_block.size() == 0) return a;
    if(a.suf_block.size() == 0) return b;
    UtilFunc util;
    CLUMPTYPE r;
    r.key_prefix = a.key_prefix;
    r.index_ID = 0;
    r.suf_block.resize(a.suf_block.size() + b.suf_block.size());
    r.lcp_block.resize(a.lcp_block.size() + b.lcp_block.size());
    int ridx = 0;
    int ia = 0;                         // the pointer for a given suffix in a
    int ib = 0;                         // the pointer for a given suffix in b
    int pab = a.key_prefix.length();    // the LCP between a[ia] and b[ib]
    LCPTYPE lcp_buffer = 0;             // the LCP from the last comparison
    bool is_prev_a = false;             // indicates whether the previous suffix comes from a
    while(ia < a.suf_block.size() && ib < b.suf_block.size()) {
        RIDTYPE doc_a = seqs.block_size_[a.index_ID] + a.suf_block[ia].doc;
        RIDTYPE doc_b = seqs.block_size_[b.index_ID] + b.suf_block[ib].doc;

        // DEBUG
        //cout << "comparison:    " << endl;
        //cout << (char*) (seqs.sequence_[doc_a] + a.suf_block[ia].pos) << endl;
        //cout << (char*) (seqs.sequence_[doc_b] + b.suf_block[ib].pos) << endl;

        pair<int, int> c = util.CmpStrWithLCP(
            (char*) (seqs.sequence_[doc_a] + a.suf_block[ia].pos), 
            (char*) (seqs.sequence_[doc_b] + b.suf_block[ib].pos), 
            pab
        );
        
        // DEBUG
        //cout << "results:   " << c.first << "   " << c.second << endl;

        if(c.first < 0 || (c.first == 0 && doc_a < doc_b))  {
            // push the suffix in a into the new clump
            r.suf_block[ridx].doc = doc_a;
            r.suf_block[ridx].pos = a.suf_block[ia].pos;
            r.lcp_block[ridx] = is_prev_a ? a.lcp_block[ia] : lcp_buffer;
            ++ ia;
            is_prev_a = true;
            lcp_buffer = c.second;                                  // record the LCP
            pab = pab > a.lcp_block[ia] ? a.lcp_block[ia] : pab;    // pick the smaller one
        }   else    {
            // push the suffix in b into the new clump
            r.suf_block[ridx].doc = doc_b;
            r.suf_block[ridx].pos = b.suf_block[ib].pos;
            r.lcp_block[ridx] = !is_prev_a ? b.lcp_block[ib] : lcp_buffer;
            ++ ib;
            is_prev_a = false;
            lcp_buffer = c.second;                                  // record the LCP
            pab = pab > b.lcp_block[ib] ? b.lcp_block[ib] : pab;    // pick the smaller one
        }
        ++ ridx;  
    }
    // copy the rest information if either clump is not exhausted
    if(ia < a.suf_block.size()) {
        int tidx = ridx;
        for(ia; ia < a.suf_block.size(); ++ ia)   {
            RIDTYPE doc_a = seqs.block_size_[a.index_ID] + a.suf_block[ia].doc;
            r.suf_block[ridx].doc = doc_a;
            r.suf_block[ridx].pos = a.suf_block[ia].pos;
            r.lcp_block[ridx] = a.lcp_block[ia];
            ++ ridx;
        }
        r.lcp_block[0] = 0;                 // the first LCP value is meaningless, set it to 0
        r.lcp_block[tidx] = lcp_buffer;     // the first LCP for the first suffix directly copied, coresponds to the last comparison
    }   else if(ib < b.suf_block.size())   {
        int tidx = ridx;
        for(ib; ib < b.suf_block.size(); ++ ib)   {
            RIDTYPE doc_b = seqs.block_size_[b.index_ID] + b.suf_block[ib].doc;
            r.suf_block[ridx].doc = doc_b;
            r.suf_block[ridx].pos = b.suf_block[ib].pos;
            r.lcp_block[ridx] = b.lcp_block[ib];
            ++ ridx;
        }
        r.lcp_block[0] = 0;                 // the first LCP value is meaningless, set it to 0
        r.lcp_block[tidx] = lcp_buffer;     // the first LCP for the first suffix directly copied, coresponds to the last comparison
    }
    
    assert(r.suf_block.size() == r.lcp_block.size());
    
    // DEBUG
    //if(a.key_prefix == "YHRYLAEFAS")    {
    //    for(int i = 0; i < r.suf_block.size() - 1; ++ i)   {
    //        string str = seqs.GetSuffixSeq(r.suf_block[i]);
    //        cout << ":::    " << str << endl;
    //    }
    //}

    // DEBUG verify the merged clump is correct
    //for(int i = 0; i < r.suf_block.size() - 1; ++ i)   {
    //    string str1 = seqs.GetSuffixSeq(r.suf_block[i]);
    //    string str2 = seqs.GetSuffixSeq(r.suf_block[i + 1]);
    //    if(str1.compare(str2) > 0 || str1.substr(0, r.lcp_block[i + 1]) != str2.substr(0, r.lcp_block[i + 1]))  {
    //        cout << "Clump A:  " << endl;
    //        for(int i = 0; i < a.suf_block.size(); ++ i)   {
    //            cout << seqs.GetSuffixSeq(a.index_ID, a.suf_block[i]) << "\t" << a.lcp_block[i] << endl;
    //        }
    //        cout << "Clump B:  " << endl;
    //        for(int i = 0; i < b.suf_block.size(); ++ i)   {
    //            cout << seqs.GetSuffixSeq(b.index_ID, b.suf_block[i]) << "\t" << b.lcp_block[i] << endl;
    //        }
    //        cout << "Merged Clump:  " << endl;
    //        for(int i = 0; i < r.suf_block.size(); ++ i)   {
    //            cout << seqs.GetSuffixSeq(r.suf_block[i]) << "\t" << r.lcp_block[i] << endl;
    //        }
    //        exit(0);
    //    }
    //}

    return r;
}