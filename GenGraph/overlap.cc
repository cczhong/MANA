#include "overlap.h"

using namespace std;

Overlap::Overlap(void)  {
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
    //CountNumIndex(dir, file_stem, num_index);
    std::ifstream* GSAfh = new std::ifstream[num_blocks];
    std::ifstream* LCPfh = new std::ifstream[num_blocks];

    // DEBUG
    //cout << "Good here 0" << endl;

    OpenIndexFiles(dir, file_stem, num_blocks, GSAfh, LCPfh);

    // DEBUG
    //cout << "Good here 1" << endl;

    Clump suf_clump;
    suf_clump.InitClumpPQueue(num_blocks, GSAfh, LCPfh, seqs, min_overlap);    // initilization

    bool t;
    do
    {
        t = suf_clump.NextClump(GSAfh, LCPfh, seqs, min_overlap);
        if(t)   {
            ResolvePerfectOverlap(suf_clump.current_, seqs);
        }
    } while (t);
    


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
        cout << gsa_file << endl;
        GSAfh[i].open(gsa_file, std::ios::in | std::ios::binary);
        cout << "file opened" << endl;
        if (!GSAfh[i].good()) {
			std::cerr << "MANA::GenGraph::Overlap::OverlapIndexFIles: Cannot open GSA index file: " << gsa_file << "\n";
			exit (1);
		}
        GSAfh[i].read((char*) &foo_size, sizeof(SFAIDXTYPE));   // get rid of the header information
        cout << "info read: " << foo_size << endl;

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


void Overlap::ResolvePerfectOverlap(std::vector<CLUMPTYPE>& m_clump, SFABuild &seqs)    {
    assert(m_clump.size() > 0);
    for(int i = 0; i < m_clump.size(); ++ i)   {
        assert(m_clump[i].suf_block.size() > 1);
        assert(m_clump[i].suf_block.size() == m_clump[i].lcp_block.size());
        vector<GSATYPE> suf_stack(m_clump[i].suf_block.size());  // stack that contains full suffixes
        vector<LCPTYPE> lcp_stack(m_clump[i].lcp_block.size());  // stach that contains the length of the lcp
        int s_index = 0;            // the size of the current index
        for(int j = 1; j < m_clump[i].suf_block.size(); ++ j)   {
            // maintain stack validity
            while(s_index > 0 && m_clump[i].lcp_block[j] < lcp_stack[s_index]) {
                -- s_index;     // mimic stack popping, this source read becomes invalid
            }
            // detect complete suffix overlap
            if(m_clump[i].lcp_block[j] == seqs.GetSufLen(m_clump[i].index_ID, m_clump[i].suf_block[j]))    {
                // TODO: more adavanced insert function to eliminate redundancy
                suf_stack[s_index] = m_clump[i].suf_block[j - 1];
                lcp_stack[s_index] = m_clump[i].lcp_block[j - 1];
                ++ s_index;
            }
            // detect complete prefix overlap
            if(m_clump[i].suf_block[j].pos == 0)    {
                // it should overlap with all suffixes in the stack
                for(int k = 0; k < s_index; ++ k)   {
                    if(lcp_stack[k] < seqs.GetSufLen(m_clump[i].index_ID, m_clump[i].suf_block[j])) {
                        RIDTYPE source = seqs.GetFullRID(m_clump[i].index_ID, suf_stack[k].doc);
                        RIDTYPE target = seqs.GetFullRID(m_clump[i].index_ID, m_clump[i].suf_block[j].doc);
                        cout << "Overlap detected:  " << endl;
                        cout << seqs.GetSequence(source) << endl;
                        cout << seqs.GetSequence(target) << endl;
                    }
                }
            }
        }
    }
    return;
}