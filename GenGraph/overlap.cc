#include "overlap.h"

using namespace std;

Overlap::Overlap(void)  {
    is_overlap_init_ = false;
    return;
}

Overlap::~Overlap(void) {
    return;
}

void Overlap::DetectOverlaps(
    SFABuild& seqs,                 // the original set of sequences
    const std::string& dir,         // the folder that contains the index files
    const std::string& file_stem,   // the file stem
    const int& min_overlap          // the minimum overlap length
)   {
    int num_blocks = seqs.GetNumBlocks();
    seqs.InitContainedInfo();

    //CountNumIndex(dir, file_stem, num_index);
    std::ifstream* GSAfh = new std::ifstream[num_blocks];
    std::ifstream* LCPfh = new std::ifstream[num_blocks];

    // DEBUG
    //cout << "Good here 0" << endl;

    OpenIndexFiles(dir, file_stem, num_blocks, GSAfh, LCPfh);

    // DEBUG
    //cout << "Good here 1" << endl;


    // the overlap information needs to be dumped into disk
    string overlap_file = dir + "/" + file_stem + ".olp";
    std::ofstream OLPfh(overlap_file, ios_base::out | ios_base::binary);
    if (!OLPfh.good()) {
		std::cerr << "MANA::GenGraph::Overlap::DetectOverlaps: Cannot write into overlap file: " << overlap_file << "\n";
		exit (1);
	}
        
    Clump suf_clump;
    suf_clump.InitClumpPQueue(num_blocks, GSAfh, LCPfh, seqs, min_overlap);    // initilization

    bool t;
    do
    {
        t = suf_clump.NextClump(GSAfh, LCPfh, seqs, min_overlap);
        if(t)   {
            ResolvePerfectOverlap(suf_clump.current_, seqs, OLPfh);
        }
    } while (t);
    
    OLPfh.close();

    // DEBUG
    //seqs.PrintContainedInfo();

    //for(int i = 0; i < num_index; ++ i)   {
    //    bool success;
    //    do
    //    {
    //        success = GetSuffixClump(GSAfh[i], LCPfh[i], 10, clump, buffer);
            //cout << "is success:    " << success << endl;
    //    } while (success);
         
        
        //if(success) {
        //    for(int j = 0; j < clump.size(); ++ j)   {
        //        cout << "DOC:   " << clump[j].doc << "  POS:    " << clump[j].pos << endl;
        //    }
        //}
    //}
    for(int i = 0; i < num_blocks; ++ i)   {
        GSAfh[i].close();
        LCPfh[i].close();
    }
    delete [] GSAfh;
    delete [] LCPfh;
    return;
}

void Overlap::CountNumIndex(
    const std::string& dir,             // the folder that contains the index files
    const std::string& file_stem,       // the file stem to load 
    int &num_index                      // (output) the number of index pieces
)    {
    num_index = 0;
    while(true) {
        string GSAfile = dir + "/" + file_stem + "." + std::to_string(num_index) + ".gsa";
        string LCPfile = dir + "/" + file_stem + "." + std::to_string(num_index) + ".lcp";
        if(!isFileExists(GSAfile) || !isFileExists(LCPfile))    {
            break;
        }
        ++ num_index;        
    }
    return;
}

void Overlap::OpenIndexFiles(
    const std::string& dir,                 // the folder that contains the index files
    const std::string& file_stem,           // the file stem to load
    const int num_index,                    // the number of indexes to load
    std::ifstream* GSAfh,  // (output) the array that contains the file handles to the generalized suffix array indexes
    std::ifstream* LCPfh   // (output) the array that contains the file handles to the LCPs
)   {

    SFAIDXTYPE foo_size;

    for(int i = 0; i < num_index; ++ i)   {
        string gsa_file = dir + "/" + file_stem + "." + std::to_string(i) + ".gsa";
        //DEBUG
        //cout << gsa_file << endl;
        
        GSAfh[i].open(gsa_file, std::ios::in | std::ios::binary);
        if (!GSAfh[i].good()) {
			std::cerr << "MANA::GenGraph::Overlap::OverlapIndexFIles: Cannot open GSA index file: " << gsa_file << "\n";
			exit (1);
		}
        GSAfh[i].read((char*) &foo_size, sizeof(SFAIDXTYPE));   // get rid of the header information
        
        // DEBUG
        //cout << "info read: " << foo_size << endl;

        string lcp_file = dir + "/" + file_stem + "." + std::to_string(i) + ".lcp";
        LCPfh[i].open(lcp_file, std::ios::in | std::ios::binary);
        if (!LCPfh[i].good()) {
			std::cerr << "MANA::GenGraph::Overlap::OverlapIndexFIles: Cannot open LCP index file: " << lcp_file << "\n";
			exit (1);
		}

        //cout << "loaded:    " << i << endl;
    }
    
    return;
}


void Overlap::ResolvePerfectOverlap(const CLUMPTYPE& clump, SFABuild& seqs, std::ofstream& OLPfh)    {
    
    assert(OLPfh.good());
    assert(seqs.is_contained_init_);
    assert(clump.suf_block.size() == clump.lcp_block.size());
    
    if(clump.suf_block.size() <= 1) return;

    vector<GSATYPE> suf_stack(clump.suf_block.size());  // stack that contains full suffixes
    vector<LCPTYPE> lcp_stack(clump.lcp_block.size());  // stach that contains the LCP of the prefix of the source reads
    vector<LCPTYPE> len_stack(clump.lcp_block.size());  // stack that contains the length of the suffix
    int s_index = 0;                                    // the size of the current index

    LCPTYPE prev_prefix_len = 0;                        // the length of the previous prefix (target)
    GSATYPE prefix_candidate;                           // a possible prefix (target) candidate, need to check if contained 
    bool is_prefix_found = false;                       // indicate whether a prefix (target) has been detected

    // DEBUG
    //if(clump.key_prefix == "TDEEAKKLFP")    {
    //    cout << "===============CLUMP BEGIN===============" << clump.key_prefix << endl;
    //    for(int j = 0; j < clump.suf_block.size(); ++ j)   {
    //        RIDTYPE a = seqs.GetFullRID(clump.index_ID, clump.suf_block[j].doc);
    //        cout << seqs.GetSequence(a) << "    " << a << endl;
    //    }
        //exit(0);
    //}

    // check if first sequence is a prefix to consider the case it is contained
    if(clump.suf_block[0].pos == 0)    {
        is_prefix_found = true;
        prefix_candidate = clump.suf_block[0];
        prev_prefix_len = seqs.GetSufLen(clump.index_ID, clump.suf_block[0]);
    }   
    // EXPERIMENT
    //else    {
        // otherwise if the first suffix is a regular suffix, insert it into the stack
    //    InsertSource(clump.index_ID, 
    //        clump.suf_block[0], clump.lcp_block[0], seqs.GetSufLen(clump.index_ID, clump.suf_block[0]),
    //        seqs, suf_stack, lcp_stack, len_stack, s_index
    //    );
    //}
    

    int last_target = clump.suf_block.size() - 1;
    // EXPERIMENT: whether two passes with pre-termination would work better
    //for(int j = 1; j < clump.suf_block.size(); ++ j)   {
    //    // OK to terminate earlier if no more target read presents
    //    if(clump.suf_block[j].pos == 0) {   last_target = j;    }   
    //}
    
    for(int j = 1; j <= last_target; ++ j)   {

        RIDTYPE rid_current = seqs.GetFullRID(clump.index_ID, clump.suf_block[j].doc);
        // DEBUG
        //if(clump.key_prefix == "TDEEAKKLFP") {
        //    cout << "*** working on read:   " << rid_current << endl;   
        //}

        if(seqs.is_contained_[rid_current])    {
            // DEBUG
            //if(clump.key_prefix == "TDEEAKKLFP") cout << "the read is contained: " << rid_current << endl;
            
            // needs to update suffix stack as well
            while(s_index > 0 && clump.lcp_block[j] < len_stack[s_index - 1]) {
                // DEBUG
                //if(clump.key_prefix == "EVCPAGWTPG") cout << "Source suffix invalid " << s_index << endl;
                -- s_index;     // mimic stack popping, this source read becomes invalid
            }
            continue;
        }


        // check if the candidate prefix is contained
        // condition: if the LCP is greater than or equal to the length of the candidate prefix
        if(is_prefix_found && clump.lcp_block[j] >= prev_prefix_len)    {
            // DEBUG
            //if(clump.key_prefix == "TDEEAKKLFP") cout << "Prefix found && contained " << j << endl;

            RIDTYPE c = seqs.GetFullRID(clump.index_ID, prefix_candidate.doc);
            RIDTYPE d = seqs.GetFullRID(clump.index_ID, clump.suf_block[j].doc);
            seqs.is_contained_[c] = true;
            seqs.contained_by_[c] = d;
            is_prefix_found = false;       // the candidate prefix becomes invalid as it is contained
        }   else if(is_prefix_found && clump.lcp_block[j] < prev_prefix_len)    {
            // DEBUG
            //if(clump.key_prefix == "TDEEAKKLFP") cout << "Prefix found && not contained && overlap recorded " << j << endl;

            // detect complete prefix overlap
            // condition, the LCP drops to less than the candidate prefix length (indicating that the candidate is unique)
            for(int k = 0; k < s_index; ++ k)   {
                // only record if the target prefix contributes to unique sequence in extension
                if(len_stack[k] < prev_prefix_len) {
                    RIDTYPE source = seqs.GetFullRID(clump.index_ID, suf_stack[k].doc);
                    RIDTYPE target = seqs.GetFullRID(clump.index_ID, prefix_candidate.doc);
                    if(!seqs.is_contained_[source] && !seqs.is_contained_[target])    {
                        // calculate the length of the merged sequence
                        OLPfh.write((char*) &source, sizeof(RIDTYPE));
                        OLPfh.write((char*) &target, sizeof(RIDTYPE));
                        LCPTYPE mlen = seqs.GetSeqLen(source) + seqs.GetSeqLen(target) - len_stack[k];
                        OLPfh.write((char*) &mlen, sizeof(LCPTYPE));

                        // DEBUG
                        //if(clump.key_prefix == "TDEEAKKLFP") cout << "overlap recorded  " << source << "    " << target << endl;
                        
                        //OVERLAPTYPE o;
                        //o.source = source; o.target = target; o.len = len_stack[k];
                        //overlap_fw_[source].insert(std::make_pair(target, o));
                        //overlap_re_[target].insert(std::make_pair(source, o));
                    }
                    //cout << "Overlap detected:  " << source << "->" << target << endl;
                    //cout << seqs.GetSequence(source) << endl;
                    //cout << seqs.GetSequence(target) << endl;
                }
            }
            is_prefix_found = false;        // the candidate prefix becomes invalid as it has been used for overlap
        }

        // maintain stack validity (the source is no longer valid because of LCP decrease)
        while(s_index > 0 && clump.lcp_block[j] < len_stack[s_index - 1]) {
            // DEBUG
            //if(clump.key_prefix == "TDEEAKKLFP") cout << "Source suffix invalid " << s_index << endl;
            -- s_index;     // mimic stack popping, this source read becomes invalid
        }

        // DEBUG
        //if(clump.key_prefix == "TDEEAKKLFP")    {
        //    cout << "LCP info:  " << clump.lcp_block[j] << endl;
        //    for(int x = 0; x < s_index; ++ x)   {
        //        cout << "???    " << seqs.GetSequence(clump.index_ID, suf_stack[x].doc) << "    " << len_stack[x] << endl;
        //    }        
        //}
        // if the previous sequence is a candidate suffix (source), record it
        if(clump.lcp_block[j] >= seqs.GetSufLen(clump.index_ID, clump.suf_block[j - 1]) && clump.suf_block[j - 1].pos > 0)    {
            // DEBUG
            //if(clump.key_prefix == "TDEEAKKLFP") cout << "Source suffix recorded " << j << endl;

            InsertSource(clump.index_ID, 
                clump.suf_block[j - 1], clump.lcp_block[j - 1], seqs.GetSufLen(clump.index_ID, clump.suf_block[j - 1]),
                seqs, suf_stack, lcp_stack, len_stack, s_index
            );

            //suf_stack[s_index] = clump.suf_block[j - 1];
            //len_stack[s_index] = seqs.GetSufLen(clump.index_ID, clump.suf_block[j - 1]);
            //++ s_index;
        }

        // if the current sequence is a prefix, record it
        if(clump.suf_block[j].pos == 0)    {
            // DEBUG
            //if(clump.key_prefix == "TDEEAKKLFP") cout << "Target prefix recorded " << j << endl;
            if(clump.lcp_block[j] >= seqs.GetSufLen(clump.index_ID, clump.suf_block[j]) && clump.suf_block[j - 1].pos > 0)    {
                // the current read is contained by the previous read
                RIDTYPE c = seqs.GetFullRID(clump.index_ID, clump.suf_block[j].doc);
                RIDTYPE d = seqs.GetFullRID(clump.index_ID, clump.suf_block[j - 1].doc);
                seqs.is_contained_[c] = true;
                seqs.contained_by_[c] = d;
            }   else    {
                // record as a candidate target read
                is_prefix_found = true;
                prefix_candidate = clump.suf_block[j];
                prev_prefix_len = seqs.GetSufLen(clump.index_ID, clump.suf_block[j]);
            }
        }

    }
    // handle the last sequence if it is a prefix
    if(is_prefix_found)    {
        for(int k = 0; k < s_index; ++ k)   {
            // only record if the target prefix contributes to unique sequence in extension
            if(len_stack[k] < prev_prefix_len) {
                
                RIDTYPE source = seqs.GetFullRID(clump.index_ID, suf_stack[k].doc);
                RIDTYPE target = seqs.GetFullRID(clump.index_ID, prefix_candidate.doc);
                
                if(!seqs.is_contained_[source] && !seqs.is_contained_[target])    {
                    OLPfh.write((char*) &source, sizeof(RIDTYPE));
                    OLPfh.write((char*) &target, sizeof(RIDTYPE));
                    LCPTYPE mlen = seqs.GetSeqLen(source) + seqs.GetSeqLen(target) - len_stack[k];
                    OLPfh.write((char*) &mlen, sizeof(LCPTYPE));
                    
                    // DEBUG
                    //if(clump.key_prefix == "TDEEAKKLFP") cout << "overlap recorded  " << source << "    " << target << endl;

                    //OVERLAPTYPE o;
                    //o.source = source; o.target = target; o.len = len_stack[k];
                    //overlap_fw_[source].insert(std::make_pair(target, o));
                    //overlap_re_[target].insert(std::make_pair(source, o));
                }
                //cout << "Overlap detected:  " << source << "->" << target << endl;
                //cout << seqs.GetSequence(source) << endl;
                //cout << seqs.GetSequence(target) << endl;
            }
        }
    }

    // DEBUG
    // "AAGFLTRDSR"
    //if(clump.key_prefix == "TDEEAKKLFP")    {
    //    exit(0);
    //}
    return;
}

void Overlap::DetectUnitigs (
    const std::string& dir,                 // the index folder
    const std::string& file_stem,           // the file stem
    SFABuild& seqs 
)   {
    string olp_file = dir + "/" + file_stem + ".olp";
    std::vector<std::vector<PARTIALOLPTYPE> > olp_info;
    LoadOverlapInfo(olp_file, seqs, olp_info);
    return;
}

void Overlap::LoadOverlapInfo(
    const std::string& file,
    SFABuild& seqs,
    std::vector<std::vector<PARTIALOLPTYPE> >& olp_info
)   {
    assert(seqs.is_contained_init_);

    olp_info.resize(seqs.num_seqs_);
    std::ifstream in_fh(file, std::ios::in | std::ios::binary);
    if(!in_fh.good())  {
        cout << "MANA::GenGraph::Overlap::LoadOverlapInfo: Cannot read overlap index file " << file << endl;
        exit(1);
    }

    while(true) {
        if(in_fh.peek() == EOF) {   break;  }
        
        // loads in the information
        RIDTYPE source, target;
        LCPTYPE len;
        in_fh.read((char*) &source, sizeof(RIDTYPE));
        in_fh.read((char*) &target, sizeof(RIDTYPE));
        in_fh.read((char*) &len, sizeof(LCPTYPE));

        // record the information is it does not involve any contained read
        if(!seqs.is_contained_[source] && !seqs.is_contained_[target])  {
            PARTIALOLPTYPE po;
            po.target = target; po.len = len; po.is_transitive = false;
            olp_info[source].push_back(po);
        }

    }

    // DEBUG
    for(int i = 0; i < seqs.num_seqs_; ++ i)   {
        for(int j = 0; j < olp_info[i].size(); ++ j)   {
            cout << "======" << endl;
            cout << "ID:    " << seqs.GetHeader(i) << "\t" << seqs.GetHeader(olp_info[i][j].target) << endl;
            cout << "source:    " << seqs.GetSequence(i) << endl;
            cout << "target:    " << seqs.GetSequence(olp_info[i][j].target) << endl;
            cout << "len:   " << olp_info[i][j].len << endl;
        }
    }

    in_fh.close();
    return;
}  

void Overlap::InsertSource(
    const int& index_ID,
    const GSATYPE& suf,
    const LCPTYPE& lcp,                 // note that this is the lcp of the forward direction
    const LCPTYPE& len,
    SFABuild& seqs,
    std::vector<GSATYPE>& suf_stack,
    std::vector<LCPTYPE>& lcp_stack,    // note that this is the lcp of the reverse direction
    std::vector<LCPTYPE>& len_stack,
    int& s_index
)   {
    
    // DEBUG
    //cout << "=== begin of insert ===    " << suf.doc << "   " << lcp << "   " << len << "   " << s_index << endl;  

    assert(seqs.is_contained_init_);
    assert(suf_stack.size() == lcp_stack.size());
    assert(lcp_stack.size() == len_stack.size());

    if(s_index == 0)    {
        suf_stack[s_index] = suf;
        lcp_stack[s_index] = 0;
        len_stack[s_index] = len;
        ++ s_index;
        return;
    }

    // indicate whether we can delete this suffix
    std::vector<bool> del_mark(suf_stack.size(), false);
    string str = seqs.GetSequence(index_ID, suf.doc);
    str.resize(suf.pos);
    RIDTYPE c = seqs.GetFullRID(index_ID, suf.doc);

    // DEBUG
    //cout << "sequence:  " << str << "   " << c << endl;

    UtilFunc util;
    LCPTYPE pre_lcp = 0;        // the number of letter we can safely skip
    bool is_last_set = false;   // indicate whether the last lcp is computed
    LCPTYPE last_lcp;
    for(int i = s_index - 1; i >= 0; -- i)   {
        RIDTYPE d = seqs.GetFullRID(index_ID, suf_stack[i].doc);
        string t_str = seqs.GetSequence(index_ID, suf_stack[i].doc);
        t_str.resize(suf_stack[i].pos);

        if(seqs.is_contained_[d])    {
            del_mark[i] = true;  
            continue;   
        }

        // DEBUG
        //cout << "cmp sequence:  " << t_str << " " << d << endl;

        std::pair<int, int> r = util.CmpWithLCPRev(str.c_str(), t_str.c_str(), pre_lcp);

        // DEBUG
        //cout << "Compare results:   "   << r.first << " " << r.second << endl;
        //cout << "skip LCP:  " << pre_lcp << endl;
        //cout << "LCP:   " << lcp << endl;
        //cout << "length suffix: " << len_stack[i] << endl;


        if(r.first == -1 && r.second >= str.length() && len <= len_stack[i])    {
            // DEBUG
            //cout << "the current read is contained" << endl;
            //if(c == 10001)  {
            //    cout << "Here 1!!!" << endl;
            //    exit(0);
            //}
            // the current read is contained by the existing read
            seqs.is_contained_[c] = true;  
            seqs.contained_by_[c] = d;
            return;  // search no more because the current read is no longer valid
        }   else if(r.first >= 0 && r.second >= t_str.length())   {
            // DEBUG
            //cout << "the existing read is contained" << endl;
            //if(d == 10001)  {
            //    cout << "Here 2!!!" << endl;
            //    exit(0);
            //}
            // the existing read is contained
            seqs.is_contained_[d] = true;
            seqs.contained_by_[d] = c;
        }   else if(r.first == -1 && r.second >= str.length() && len > len_stack[i])   {
            // DEBUG
            //cout << "the existing read needs to be deleted" << endl;

            // the existing read is covered
            del_mark[i] = true;
        }
        if(!is_last_set) last_lcp = r.second;
        pre_lcp = r.second < lcp_stack[i] ? r.second : lcp_stack[i];
    }

    

    // delete all marked entries
    int n_index = 0;
    for(int i = 0; i < s_index; ++ i)   {
        if(!del_mark[i] && n_index < i)    {
            suf_stack[n_index] = suf_stack[i];
            lcp_stack[n_index] = lcp_stack[i];
            len_stack[n_index] = lcp_stack[i];
            ++ n_index;
        }   else if(!del_mark[i] && n_index == i)   {
            ++ n_index;
        }
        // if the entrie is marked for deletion then nothing needs to done (just let i increase)
    }
    s_index = n_index;

    // DEBUG
    //cout << "index: " << s_index << endl;
    //cout << "======out of loop" << endl;

    // finally insert the current entrie
    suf_stack[s_index] = suf;
    lcp_stack[s_index] = last_lcp;
    len_stack[s_index] = len;
    ++ s_index;

    return;
}