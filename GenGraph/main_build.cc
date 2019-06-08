#include "align_batch.h"
#include "bwt.h"
#include "bwt_search.h"
#include "loader.h"
#include "bio_alphabet.h"
#include "kmer_unitcoder.h"
#include "minimizer_sort.h"
#include "string_graph.h"
#include "sequence_search.h"
#include "kmer_unitcoder.h"
#include "scoring_prot.h"
#include "reduced_alphabet.h"
#include "kmer_filtering.h"
#include "sfa_build.h"
#include "database_index.h"
#include "util_func.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

static string workspace_dir = "index";
static string db_file;
static int extd_len;
static int num_threads = 1;
static int neighbor_score = 11;
static string verbose;
static int scoring_matrix = 0;
static int mer_len = 3;

static int seed_len;
static int alph_id;
static int min_seed_coverage;
static int min_ext_coverage;

void PrintUsage()  {
  cout << "Usage: grasp-build [peptide_db (FASTA)]" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
}

double MyTime (void)
{
    int flag;
    clockid_t cid = CLOCK_REALTIME; // CLOCK_MONOTONE might be better
    timespec tp;
    double timing;
	
    flag = clock_gettime(cid, &tp);
    if (flag == 0) timing = tp.tv_sec + 1.0e-9*tp.tv_nsec;
    else           timing = -17.0;         // If timer failed, return non-valid time
	
    return(timing);
}

void PrintElapsed( double s, double e, const char *task )
{
	double elapsed = e - s ;
	printf ("[%02.0f:%02.0f:%05.2f]\t%s\n", 
			floor(elapsed/3600.0), 
			floor(fmod(elapsed,3600.0)/60.0), 
			fmod(elapsed,60.0),
			task);
	return;
}


int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("db_file", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("work_space", boost::program_options::value<string>(&workspace_dir)->default_value("index"), "working directory for indexing file dump")
      ("seed_coverage", boost::program_options::value<int>(&min_seed_coverage)->default_value(3), "minimum coverage required for seed sequence")
      ("extension_coverage", boost::program_options::value<int>(&min_ext_coverage)->default_value(1), "minimum coverage required for extending assembly")
      ("seed_len", boost::program_options::value<int>(&seed_len)->default_value(6), "length of the seeds")
      ("extension_len", boost::program_options::value<int>(&extd_len)->default_value(10), "minimum overlap length for path extension")
      ("alphabet", boost::program_options::value<int>(&alph_id)->default_value(4), "reduced alphabet to be used for seeding\n  0:ALL20  1:DSSP5  2:DSSP10  3:GBMR4\n  4:GBMR10  5:HSDM5  6:SDM6  7:MURPHY5\n  8:MURPHY10  9:TD5  10:TD10")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("db_file", 1);
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  if(vm.count("help"))  {
    cout << "Usage: graps-build db_file(expecting FASTA-format)\n" << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  boost::filesystem::path abs_workspace = workspace_dir;
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: grasp-build: db_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::is_directory(workspace_dir))  {
    cout << workspace_dir << endl;
    cout << "Error: grasp-build: working space does not exist (please provide full path)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(seed_len < 3 || seed_len > 10)  {
    cout << "Error: grasp-build: seed length out of range (allowed range: 3-10)." << endl;
    exit(0);
  }
  if(extd_len < 6 || seed_len > 20)  {
    cout << "Error: grasp-build: extension length out of range (allowed range: 6-20)." << endl;
    exit(0);
  }
  if(alph_id < 0 || alph_id > 10)  {
    cout << "Error: grasp-build: reduced alphabet not supported (supported alphabets: 0-10)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  cout << "============================================================" << endl;
  double start_time = mytime();
  // load in sequence from file for forward sequences
  UtilFunc util;
  string db_stem = util.GetFileStem(db_file);
  SFABuild db_seq(db_file);
  db_seq.BuildSFAMulti(500000, workspace_dir, db_stem);
  //db_seq.BuildSFAMulti(15000000, workspace_dir, db_stem);
  //db_seq.BuildSFAMulti(5, workspace_dir, db_stem);
  cout << " Index file written." << endl;
  cout << "============================================================" << endl;

  // DEBUG
  //db_seq.PrintAllSuffixes();

  //db_seq.BuildSFADefault();
  //db_seq.BuildKeyArrayDefault();
  //double fw_sfa_build_time = mytime();
  //printElapsed(start_time, fw_sfa_build_time, "GRASPx::Build info: finished building forward suffix array");
  //start_time = mytime();
  // dump the forward version suffix array
  //db_seq.DumpSFA(workspace_dir, db_stem);
  //double fw_sfa_dump_time = mytime();
  //printElapsed(start_time, fw_sfa_dump_time, "GRASPx::Build info: finished writing forward suffix array");
  
  /*  
  start_time = mytime();
  // prepare DatabaseIndex object
  DatabaseIndex db_dump(alph_id, seed_len, extd_len);
  //db_dump.BuildSeedmerMap(db_seq);
  unordered_map<std::string, std::list<std::string> > reduc_alph_map;
  //db_dump.CreateReducedMap(reduc_alph_map);
  string out_file = workspace_dir + "/" + db_stem + ".rdm";
  //db_dump.DumpReducedMap(reduc_alph_map, out_file);
  // create forward seed extensions
  unordered_map<std::string, std::list<PositionType> > fw_seed_ext;
  db_dump.GetSeedExt(db_seq, min_seed_coverage, fw_seed_ext);
  // create forward read extensions
  std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> > fw_read_ext;
  db_dump.CreateReadExtWorker(min_ext_coverage, db_seq, fw_read_ext);
  double fw_ext_time = mytime();
  printElapsed(start_time, fw_ext_time, "GRASPx::Build info: finished resolving forward extensions");
  start_time = mytime();
  // finish using forward suffix array for now, destruct it to save memory
  db_seq.DestructSFA();
  // build the reverse suffix array first, as we can discarded immediately
  // after resolving extensions
  SFABuild db_seq_rev(db_seq);
  db_seq_rev.InPlaceReverse();
  // build reverse suffix array 
  db_seq_rev.BuildSFADefault();
  db_seq_rev.BuildKeyArrayDefault();
  double re_sfa_build_time = mytime();
  printElapsed(start_time, re_sfa_build_time, "GRASPx::Build info: finished building reverse suffix array");
  start_time = mytime();
  // create reverse seed_extensions
  unordered_map<std::string, std::list<PositionType> > re_seed_ext;
  db_dump.GetSeedExtRev(db_seq_rev, min_seed_coverage, re_seed_ext);
  // create reverse read extensions
  std::unordered_map<RIDTYPE, std::list<OVERLAPTYPE> > re_read_ext;
  db_dump.CreateReadExtWorker(min_ext_coverage, db_seq_rev, re_read_ext);
  db_seq_rev.DestructSFA();
  db_seq_rev.DestructSequences();
  double re_ext_time = mytime();
  printElapsed(start_time, re_ext_time, "GRASPx::Build info: finished resolving reverse extensions");
  start_time = mytime();
  // merge seed extension pairs and dump read extension
  // we need to load forward SFA to merge the seed paris
  db_seq.LoadSFA(workspace_dir, db_stem);
  unordered_map<string, list<READPAIRTYPE> > seed_pair_ext;
  db_dump.MatchSeedPair(db_seq, fw_seed_ext, re_seed_ext, seed_pair_ext);
  double pair_time = mytime();
  printElapsed(start_time, pair_time, "GRASPx::Build info: finished bridging seed pairs");
  start_time = mytime();
  //db_dump.CreateSeedExt(min_seed_coverage, db_seq, db_seq_rev, seed_pair_ext);
  string out_file_seed_ext = workspace_dir + "/" + db_stem + ".sxt";
  db_dump.DumpSeedExt(seed_pair_ext, out_file_seed_ext);
  double sxt_dump_time = mytime();
  printElapsed(start_time, sxt_dump_time, "GRASPx::Build info: finished writing seed extensions");
  start_time = mytime();
  // dump read extension
  string out_file_read_ext = workspace_dir + "/" + db_stem + ".rxt";
  db_dump.DumpReadExt(fw_read_ext, re_read_ext, out_file_read_ext);
  double rxt_dump_time = mytime();
  printElapsed(start_time, rxt_dump_time, "GRASPx::Build info: finished writing read extensions");
  start_time = mytime();
  // compute and dump high-scoring k-mers
  unordered_map<string, list<string> > hs_mer;
  //db_dump.CreateHighScoreMer(hs_mer);
  string out_file_hs_mer = workspace_dir + "/" + db_stem + ".hsm";
  //db_dump.DumpHighScoreMer(hs_mer, out_file_hs_mer);
  double hsmer_time = mytime();
  printElapsed(start_time, hsmer_time, "GRASPx::Build info: finished indexing high-scoring k-mer matches");
  cout << "============================================================" << endl;
  */
  
  /*
  // build forward suffix array
  db_seq.BuildSFADefault();
  db_seq.BuildKeyArrayDefault();
  // dump the forward version suffix array
  db_seq.DumpSFA(workspace_dir, db_stem);
  // load in sequence from file for reverse sequences
  
  // create seed-mer map for given database
  DatabaseIndex db_dump(alph_id, seed_len, extd_len);
  
  
  db_dump.BuildSeedmerMap(db_seq);
  
  double sfa_build_time = mytime();
  printElapsed(start_time, sfa_build_time, "build suffix array");
  start_time = mytime();
  
  // create and dump map for reduced sequence to original sequence
  unordered_map<std::string, std::list<std::string> > reduc_alph_map;
  db_dump.CreateReducedMap(reduc_alph_map);
  string out_file = workspace_dir + "/" + db_stem + ".rdm";
  db_dump.DumpReducedMap(reduc_alph_map, out_file);
  
  double rdc_build_time = mytime();
  printElapsed(start_time, rdc_build_time, "build reduced-alphabet mapping");
  start_time = mytime();
  
  // create and dump seed extensions
  unordered_map<string, list<READPAIRTYPE> > seed_pair_ext;
  db_dump.CreateSeedExt(min_seed_coverage, db_seq, db_seq_rev, seed_pair_ext);
  string out_file_seed_ext = workspace_dir + "/" + db_stem + ".sxt";
  db_dump.DumpSeedExt(seed_pair_ext, out_file_seed_ext);
  
  double sde_build_time = mytime();
  printElapsed(start_time, sde_build_time, "build seed extension");
  start_time = mytime();
  
  // create and dump maximal extension reads links
  unordered_map<RIDTYPE, list<OVERLAPTYPE> > fw_read_ext, re_read_ext;
  db_dump.CreateReadExt(min_ext_coverage, db_seq, db_seq_rev, fw_read_ext, re_read_ext);
  string out_file_read_ext = workspace_dir + "/" + db_stem + ".rxt";
  db_dump.DumpReadExt(fw_read_ext, re_read_ext, out_file_read_ext);
  
  double rde_build_time = mytime();
  printElapsed(start_time, rde_build_time, "build read extension");
  start_time = mytime();
  
  unordered_map<string, list<string> > hs_mer;
  db_dump.CreateHighScoreMer(hs_mer);
  string out_file_hs_mer = workspace_dir + "/" + db_stem + ".hsm";
  db_dump.DumpHighScoreMer(hs_mer, out_file_hs_mer);
  
  double hsm_build_time = mytime();
  printElapsed(start_time, hsm_build_time, "build high-score 3-mer map");
  start_time = mytime();
  */
  return 0;
}


/*
int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("db_file", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("index", boost::program_options::value<string>(&workspace_dir)->default_value("index"), "working directory for indexing file dump")
      ("extension_len", boost::program_options::value<int>(&extd_len)->default_value(10), "minimum overlap length for path extension")
      ("neighbor_score", boost::program_options::value<int>(&neighbor_score)->default_value(11), "neighbor score for 3-mer seed matches") 
      ("num_threads", boost::program_options::value<int>(&num_threads)->default_value(1), "maximum number of threads to be used")
      ("verbose", boost::program_options::value<string>(&verbose), "print intermediate information (default true)")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("db_file", 1);
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  if(vm.count("help"))  {
    cout << "Usage: grasp-build [peptide_db (FASTA)]" << endl << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  boost::filesystem::path abs_workspace = workspace_dir;
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: grasp-build: db_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::is_directory(workspace_dir))  {
    cout << workspace_dir << endl;
    cout << "Error: grasp-build: working space does not exist (please provide full path)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(extd_len < 6 || extd_len > 20)  {
    cout << "Error: grasp-build: extension length out of range (allowed range: 6-20)." << endl;
    exit(0);
  }
  bool is_verbose = true;
  if(verbose == "False" || verbose == "false" || verbose == "No" || verbose == "no" || verbose == "0")  {
    is_verbose = false;
  }
  
  
  if(is_verbose)  {
    cout << "============================================================" << endl;
    cout << "GRASP2-Build: Begin of program execution." << endl;
  }
  
  BioAlphabet protein_alphabet(PROT);
  ReducedAlphabet reduced_alphabet((enum Alphabet) 10);
  ScoringProt scoring_function(static_cast<enum MatrixName>(scoring_matrix), -10, -1); 
  
  // Load in the peptide sequences to be searched against
  double start_time = MyTime();
  double check_time;
  Loader pepdb_loader;
  int num_seqs = pepdb_loader.CountFastaNumSeqs(db_file.c_str());
  char **seqs = new char* [num_seqs];
  //num_seqs = pepdb_loader.LoadFasta(protein_alphabet, db_file.c_str(), header, seqs);
  pepdb_loader.LoadFasta(protein_alphabet, db_file.c_str(), seqs);
  string concat_seq; 
  //Concatenator concat_obj(seqs, num_seqs, concat_seq);
  
  //cout << "check seq: " << seqs[823927] << endl;
  //cout << "check seq: " << seqs[7792734] << endl;
  //return 0;

  // sort the reads based on minimizers
  KmerUnitcoder min_sort_mer(protein_alphabet, 6);
  MinimizerSort m_sort;
  int *order = new int [num_seqs];
  m_sort.SortSeqs(min_sort_mer, 10000000, num_seqs, seqs, order);
  
  //for(int i = 0; i < num_seqs; ++ i) {
  //  cout << seqs[i] << endl;
  //}
  //return 0;

  vector<string> concat_seqs;
  vector<int> id_begin;
  Concatenator concat_obj(seqs, num_seqs, 100000000, concat_seqs, id_begin);
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Load peptide database done.";
    PrintElapsed(start_time, check_time, "");
  }
  start_time = MyTime();
  //for(auto it = concat_seqs.begin(); it != concat_seqs.end(); ++ it) {
  //  cout << *it << endl << endl;
  //}

  // construct multiple BWT for each block of sequences
  int seq_idx = 0;
  for(auto it = concat_seqs.begin(); it != concat_seqs.end(); ++ it) {
    string bwt_prefix = workspace_dir + "/" + "bwt." + to_string(seq_idx);
    string rev_bwt_prefix = workspace_dir + "/" + "rev_bwt." + to_string(seq_idx);
    // construct forward sequence BWT and write it on hard disk
    BWT bwt;
    bwt.Construct(protein_alphabet, it->c_str());
    bwt.WriteIndex(bwt_prefix);
    bwt.Purge();
    // construct reverse sequence BWT and write it on hard disk
    //string rev_seq = string(it->rbegin(), it->rend());
    //BWT rev_bwt;
    //rev_bwt.Construct(protein_alphabet, rev_seq.c_str());
    //rev_bwt.WriteIndex(rev_bwt_prefix);
    //rev_bwt.Purge();
    ++ seq_idx;
  }
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Construct Burrows-Wheeler transformation done.";
    PrintElapsed(start_time, check_time, "");
  }
  start_time = MyTime();

  // load BWT one-by-one and conduct the search
  seq_idx = 0;
  StringGraph strG;
  BWTSearch bwt_searcher;
  //std::vector<std::vector<TargetOVERLAPTYPE> > extension;
  vector<vector<TargetOVERLAPTYPE> *> *extension = new vector<vector<TargetOVERLAPTYPE> *>;
  extension->resize(num_seqs);
  for(BWTINT i = 0; i < num_seqs; ++ i) {
    (*extension)[i] = new vector<TargetOVERLAPTYPE>;
  }
  vector<bool> *contained = new vector<bool> (num_seqs, false);
  //TargetOVERLAPTYPE **extension = new TargetOVERLAPTYPE* [num_seqs];
  //int *ext_count = new int [num_seqs];
  //cout << "num_threads: " << num_threads << endl;
  for(auto it = concat_seqs.begin(); it != concat_seqs.end(); ++ it) {
    string bwt_prefix = workspace_dir + "/" + "bwt." + to_string(seq_idx);
    //string rev_bwt_prefix = workspace_dir + "/" + "rev_bwt." + to_string(seq_idx);
    // load forward and backward BWTs from index
    BWT reload_bwt;
    reload_bwt.ConstructFromIndex(protein_alphabet, it->c_str(), bwt_prefix);
    //string rev_seq = string(it->rbegin(), it->rend());
    //BWT reload_rev_bwt;
    //reload_rev_bwt.ConstructFromIndex(protein_alphabet,rev_seq.c_str(), rev_bwt_prefix);
    //cout << "Finish loading index" << endl;
    strG.MultiComputeExtension(
        num_threads, extd_len, num_seqs, seqs, 
        id_begin[seq_idx], reload_bwt, extension, contained
    );
    //cout << "Done computing multiple extension" << endl;
    reload_bwt.Purge();
    it->clear();
    ++ seq_idx;
  }
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Compute read overlap done.";
    PrintElapsed(start_time, check_time, "");
  }
  start_time = MyTime();

  // TODO: handle the read ID that are shuffled by minimizer sorting

  vector<vector<TargetOVERLAPTYPE> *> *rev_extension = new vector<vector<TargetOVERLAPTYPE> *>;
  rev_extension->resize(num_seqs);
  for(BWTINT i = 0; i < num_seqs; ++ i) {
    (*rev_extension)[i] = new vector<TargetOVERLAPTYPE>;
  }
  strG.FillRevExtension(contained, extension, rev_extension);

  string db_stem = GetFileStem(db_file);
  string idx_unitig_file = workspace_dir + "/" + db_stem + ".utg";  

  strG.WriteUnitigsFromExtension(seqs, contained, extension, rev_extension, order, idx_unitig_file);
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Write unitigs done";
    PrintElapsed(start_time, check_time, "");
  }

  // TODO: memory collection for extension and rev_extension
  // Compute the length of the sequence represent by the edge connecting the two sequences
  //strG.ComputeEdgeLen(num_seqs, seqs, extension);
  // remove contained reads
  //vector<bool> contained(num_seqs, false);
  //strG.DetectContained(num_seqs, seqs, extension, contained);
  //strG.RemoveReducibleEdges(num_seqs, contained, extension);
  strG.ImportExtension(num_seqs, contained, extension);
  extension->clear(); delete extension;
  rev_extension->clear(); delete rev_extension;
  contained->clear(); delete contained;

  //return 0;
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Construct string graph done.";
    PrintElapsed(start_time, check_time, "");
  }

  // Post-processing of the string graph
  start_time = MyTime();

  
  strG.CheckGraph(); 
  strG.CondenseGraph(seqs);


  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Post-process string graph done.";
    PrintElapsed(start_time, check_time, "");
  }
  // Write the string graph to hard dist
  start_time = MyTime();
  //string db_stem = GetFileStem(db_file);
  //string idx_unitig_file = workspace_dir + "/" + db_stem + ".utg";  
  strG.WriteGraph(protein_alphabet, seqs, order, idx_unitig_file);
  //strG.WriteGraph(protein_alphabet, seqs, idx_unitig_file);
  strG.Purge();
  // Collect memory
  for(int idm = 0; idm < num_seqs; ++ idm) {
    delete [] seqs[idm]; 
  }
  delete [] seqs;
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Write string graph unitigs done.";
    PrintElapsed(start_time, check_time, "");
  }
  
  start_time = MyTime();
  SequenceSearch seq_search; 

  string idx_neighbor_file = workspace_dir + "/" + db_stem + ".knb"; 
  seq_search.IndexKmerNeighbor(
      mer_len, protein_alphabet, scoring_function, 
      neighbor_score, idx_neighbor_file
  );
  
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Constructing and writing k-mer index done. ";
    PrintElapsed(start_time, check_time, "");
    cout << "GRASP2-Build: End of program execution." << endl;
    cout << "============================================================" << endl;
  }
  
  return 0;
}
*/