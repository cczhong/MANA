#include "overlap.h"

using namespace std;

Overlap::Overlap(void)  {
    return;
}

Overlap::~Overlap(void) {
    return;
}

void Overlap::DetectOverlaps(
    const SFACHARTYPE** seq,    // the original set of sequences
    const std::string& dir,     // the folder that contains the index files
    const std::string& file_stem,   // the file stem
    const int num_errors        // the number of errors allowed during overlapping
)   {
    int num_index;
    CountNumIndex(dir, file_stem, num_index);
    std::fstream* GSAfh = new std::fstream[num_index];
    std::fstream* LCPfh = new std::fstream[num_index];
    OpenIndexFiles(dir, file_stem, num_index, GSAfh, LCPfh);

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
    const std::string& dir,             // the folder that contains the index files
    const std::string& file_stem,       // the file stem to load
    const int num_index,                // the number of indexes to load
    std::fstream* GSAfh,                // (output) the array that contains the file handles to the generalized suffix array indexes
    std::fstream* LCPfh                 // (output) the array that contains the file handles to the LCPs
)   {

    for(int i = 0; i < num_index; ++ i)   {
        string GSAfile = dir + "/" + file_stem + "." + std::to_string(i) + ".gsa";
        GSAfh[i].open(GSAfile, std::ios::in | std::ios::binary);
        if (!GSAfh[i]) {
			std::cerr << "MANA::GenGraph::Overlap::OverlapIndexFIles: Cannot open GSA index file: " << GSAfile << "\n";
			exit (1);
		}

        string LCPfile = dir + "/" + file_stem + "." + std::to_string(i) + ".lcp";
        LCPfh[i].open(LCPfile, std::ios::in | std::ios::binary);
        if (!LCPfh[i]) {
			std::cerr << "MANA::GenGraph::Overlap::OverlapIndexFIles: Cannot open LCP index file: " << LCPfile << "\n";
			exit (1);
		}

        cout << "loaded:    " << i << endl;
    }
    
    return;
}