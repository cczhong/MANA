#include "sfa_build.h"

using namespace std;

SFABuild::SFABuild(void)  {
  is_header_loaded_ = is_sequence_loaded_ = false;
  is_sfa_built_ = is_k_array_built_ = false; 
  is_size_counted_ = false;
  return;
}

SFABuild::SFABuild(std::string& file_name)  {
  is_header_loaded_ = is_sequence_loaded_ = is_sfa_built_ = false;
  is_size_counted_ = false;
  LoadSequences(file_name);
  return;
}

SFABuild::SFABuild(SFABuild &seq_obj)  {
  this->is_sfa_built_ = false;
  this->is_size_counted_ = false;
  this->num_seqs_ = seq_obj.num_seqs_;
  this->header_ = new char* [num_seqs_];
  this->sequence_ = new char* [num_seqs_];
  for(int i = 0; i < num_seqs_; ++ i) {
    if(seq_obj.is_header_loaded_)  {
      // copy the header
      int n = strlen(seq_obj.header_[i]) + 1;
      this->header_[i] = new char [n];
      memcpy(this->header_[i], seq_obj.header_[i], n);
    }
    if(seq_obj.is_sequence_loaded_)  {
      // copy the sequence
      int n = strlen(seq_obj.sequence_[i]) + 1;
      this->sequence_[i] = new char [n];
      memcpy(this->sequence_[i], seq_obj.sequence_[i], n);
    }
  }
  
  this->is_header_loaded_ = seq_obj.is_header_loaded_;
  this->is_sequence_loaded_ = seq_obj.is_sequence_loaded_;
  return;
}


SFABuild::~SFABuild(void) {
  if(is_sequence_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] sequence_[i];
    }
    delete [] sequence_;
  }
  if(is_header_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] header_[i];
    }
    delete [] header_;
  }
  if(is_sfa_built_)  {
    delete suffix_array_;
  }
  return;
}

// Accessing the contents in the object

void SFABuild::PrintAllSeqs(void)  {
  cout << "Sequences in object:" << endl;
  for(int i = 0; i < num_seqs_; ++ i) {
    cout << sequence_[i] << endl;
  }
  return;
}

void SFABuild::DestructSFA(void)  {
  // clear suffix array and key array to save memory
  suffix_array_->clear();
  suffix_array_->purgeSA();
  suffix_array_->purgeDoc();
  suffix_array_->purgeLCP();
  suffix_array_->purgeMLCP();
  key_array_.clear();
  is_sfa_built_ = false;
  is_k_array_built_ = false;
}

void SFABuild::DestructSequences(void) {
  if(is_sequence_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] sequence_[i];
    }
    is_sequence_loaded_ = false;
  }
  if(is_header_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] header_[i];
    }
    is_header_loaded_ = false;
  }
  num_seqs_ = 0;
  return;
}



std::string SFABuild::GetFileStem(const std::string& path)  {
  // going backward untill the '\/' character
  int i;
  for(i = path.length() - 1; i >= 0; -- i) {
    if(path[i] == '/')  {
      break;
    }
  }
  return path.substr(i + 1, path.length() - i - 1);
}

void SFABuild::DumpSFA(std::string& dir, std::string& file_stem) {
  string lcp_file = dir + "/" + file_stem + ".lcp";
  string mcp_file = dir + "/" + file_stem + ".mcp";
  string gsa_file = dir + "/" + file_stem + ".gsa";
  // ignoring SFA and DOC files
  suffix_array_->dump(
      "", "", lcp_file.c_str(), mcp_file.c_str(), gsa_file.c_str()
  );
  return;
}

void SFABuild::LoadSFA(std::string& dir, std::string& file_stem) {
  string lcp_file = dir + "/" + file_stem + ".lcp";
  string mcp_file = dir + "/" + file_stem + ".mcp";
  string gsa_file = dir + "/" + file_stem + ".gsa";
  // ignoring SFA and DOC files
  if(!boost::filesystem::exists(lcp_file) ||
     !boost::filesystem::exists(mcp_file) ||
     !boost::filesystem::exists(gsa_file))  {
    cerr << "Error: SFABuild::LoadSFA: Indexing files not found: have you run \'grasp-build\' first?" << endl;
    exit(0);
  } else  {
    // load the existing suffix array
    suffix_array_ = new GSA();
    suffix_array_->load(
        lcp_file.c_str(), mcp_file.c_str(), gsa_file.c_str()
    );
    suffix_array_->setSequences(sequence_);
    suffix_array_->setReadCount(num_seqs_);
    is_sfa_built_ = true;
  }
  return;
}

void SFABuild::CountDBSize(void) {
  if(is_size_counted_)  {
    return;
  }
  db_size_MB_ = 0;
  for(int i = 0; i < num_seqs_; ++ i) {
    db_size_MB_ += (double) strlen(sequence_[i]);
  }
  db_size_MB_ /= 1000000;
  is_size_counted_ = true;
  return;
}

double SFABuild::GetDBSizeInMegaBase(void)  {
  if(!is_size_counted_)  {
    this->CountDBSize();
  }
  return db_size_MB_;
}


// the sequences will be stored in global variables "header" and "sequence_"
void SFABuild::LoadSequences(std::string& file_name)  {
  vector<string> files_in;
  files_in.push_back(file_name);
  num_seqs_ = (unsigned int) seq::totalSequenceCount(files_in);
  header_ = new char* [num_seqs_];
  sequence_ = new char* [num_seqs_];
  seq::loadSequences(files_in, header_, sequence_, TAGSEQ);
  for(int i = 0; i < num_seqs_; ++ i) {
    int l = strlen(sequence_[i]);
    // chomp tailing non-characters
    while(!isalpha(sequence_[i][l - 1]))  {
      sequence_[i][l - 1] = '\0';
      -- l;
    }
  }
  is_header_loaded_ = is_sequence_loaded_ = true;
  return;
}

// copy sequence from one array to another
void SFABuild::CopySequences(int num_sequences, char **source, char **target) {
  target = new char* [num_sequences];
  for(int i = 0; i < num_sequences; ++ i) {
    target[i] = new char[strlen(source[i]) + 1];
    strcpy(target[i], source[i]);
  }
  is_sequence_loaded_ = true;
  num_seqs_ = num_sequences;
  return;
}

// reverse the sequences in array "sequence_"
void SFABuild::InPlaceReverse(void)  {
  for(int i = 0; i < num_seqs_; ++ i) {
    int l = strlen(sequence_[i]);
    for(int j = 0; j < floor(l / 2); ++ j) {
      swap(sequence_[i][j], sequence_[i][l - 1 - j]);
    }
  }
  return;
}

// building the suffix array on the entire set of sequences (default)
void SFABuild::BuildSFADefault(void) {
  this->BuildSuffixArray(sequence_, suffix_array_);
  return;
}

void SFABuild::BuildSFAMulti(const SFAIDXTYPE max_size, std::string &dir, std::string &file_stem)  {
  RIDTYPE begin = 0;
  RIDTYPE block_ID = 0;
  SFAIDXTYPE accum_size = 0;
  for(RIDTYPE i = 0; i < num_seqs_; ++ i) {
    accum_size += (SFAIDXTYPE) strlen(sequence_[i]);
    if(accum_size >= max_size)  {
      // build the SFA for the current sequence block
      if(is_sfa_built_) { delete suffix_array_; }
      suffix_array_ = new GSA((char**) (sequence_ + begin), i - begin + 1, true);
      is_sfa_built_ = true;
      // write the suffix array index
      std::string block_stem = file_stem + "." + std::to_string(block_ID);
      DumpSFA(dir, block_stem);
      // update the index trackers
      ++ block_ID;
      begin = i + 1;
      accum_size = 0;
    }
  }
  return;
}

// building the suffix array
void SFABuild::BuildSuffixArray(char** sequences, GSA* suffix_array) {
  suffix_array_ = new GSA(sequence_, num_seqs_, true);
  is_sfa_built_ = true;
  return;
}

// building key array that indicates the maximal-extension suffix for each entry
void SFABuild::BuildKeyArray(GSA* suffix_array, std::vector<SFAIDXTYPE>& key_array)  {
  SFAIDXTYPE n = suffix_array->getSize();
  key_array.resize(n);
  SFAIDXTYPE block_begin = 0;
  for(SFAIDXTYPE i = 0; i < n - 1; ++ i) {
    // if current suffix length greater than LCP with the next suffix
    // then the current suffix is the end of the block
    if(suffix_array->getSuffixLength(i) > suffix_array->getLcpAt(i + 1)) {
      for(SFAIDXTYPE j = block_begin; j <= i; ++ j) {
        key_array[j] = i;
      }
      block_begin = i + 1;
    }
  }
  // the last block
  for(SFAIDXTYPE j = block_begin; j < n; ++ j) {
    key_array[j] = n - 1;
  }
  return;
}

// building the default key array
void SFABuild::BuildKeyArrayDefault(void)  {
  if(!is_sfa_built_)  {
    cout << "Warning: SFABuild::BuildKeyArrayDefault attempt to build key array without suffix array built" << endl;
  }
  BuildKeyArray(suffix_array_, key_array_);
  is_k_array_built_ = true;
  return;
} 

// access a given sequence
std::string SFABuild::GetSequence(int index)  {
  if(!is_sequence_loaded_)  {
    cout << "Warning: SFABuild::GetSequence no sequence loaded" << endl;
    return string("");
  }
  if(index < 0 || index >= num_seqs_)  {
    cout << "Warning: SFABuild::GetSequence sequence index out of range" << endl;
    return string("");
  }
  return string(sequence_[index]);
}

std::string SFABuild::GetHeader(int index)  {
  if(!is_sequence_loaded_)  {
    cout << "Warning: SFABuild::GetSequence no sequence loaded" << endl;
    return string("");
  }
  if(index < 0 || index >= num_seqs_)  {
    cout << "Warning: SFABuild::GetSequence sequence index out of range" << endl;
    return string("");
  }
  return string(header_[index]);
}

// access a given suffix from the suffix array
std::string SFABuild::GetSuffixSFA(int index) {
  if(!is_sfa_built_)  {
    cout << "Warning: SFABuild::GetSuffixSFA suffix array not built" << endl;
    return string("");
  }
  if(index < 0 || index >= suffix_array_->getSize())  {
    cout << "Warning: SFABuild::GetSuffixSFA suffix index out of range" << endl;
    return string("");
  }
  return string(suffix_array_->getSuffix(index));
}

// search the suffix array within the object
std::pair<SFAIDXTYPE, SFAIDXTYPE> SFABuild::SearchSFA(std::string& search_seed) {
  if(!is_sfa_built_)  {
    cout << "Fatal Error: SFABuild::SearchSFA suffix array not built" << endl;
    exit(1);
  }
  return suffix_array_->searchWithLCPs(
      (SFACHARTYPE*) search_seed.c_str(), search_seed.length()
  );
}

// find a list of maximal extension reads within the given range
void SFABuild::GetMaxExtInfoWithinRange(
    std::pair<SFAIDXTYPE, SFAIDXTYPE>& range, 
    std::list<PositionType>& pos_list
)  {
  if(!is_k_array_built_)  {
    cout << "Fatal Error: SFABuild::GetMaxExtRIDWithRange key array not built" << endl;
    exit(1);
  }
  SFAIDXTYPE index = range.first;
  do  {
    index = key_array_[index] + 1;
    PositionType posT;
    posT.rid = (RIDTYPE) suffix_array_->getId(index - 1);
    posT.pos = (LCPTYPE) suffix_array_->getPos(index - 1);
    pos_list.push_back(posT);
  } while (index <= range.second);
  return;
}
