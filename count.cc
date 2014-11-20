#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>

#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_iterator.hpp>

using jellyfish::mer_dna;
typedef jellyfish::stream_manager<char**> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>> sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**>> read_parser;
 
/* ========== mer_counter ========== */
class mer_counter : public jellyfish::thread_exec {
  bool            canonical_;
  int             nb_threads_;
  mer_hash&       ary_;
  stream_manager  streams_;
  sequence_parser parser_;

public:
  mer_counter(int nb_threads, mer_hash& ary, char** file_begin, char** file_end, bool canonical) :
    canonical_(canonical),
    nb_threads_(nb_threads),
    ary_(ary),
    streams_(file_begin, file_end, 1), // 1: parse one file at a time
    parser_(mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
  { }

  virtual void start(int thid) {
    for(mer_iterator mers(parser_, canonical_) ; mers; ++mers) {
      ary_.add(*mers, 1);
    }
    ary_.done();
  }

};

/* ========== main ========== */
namespace po = boost::program_options;

int main(int argc, char *argv[]) {

  /* provided by Guillaume */
  const int counter_length = 7;
  const int max_reprobe = 126;
  const bool canonical = false;

  /* default values */
  bool kmer_stats = false;
  int num_threads = 1;

  /* manadatory arguments */
  int kmer_length;
  unsigned long hash_size;
  std::string reads;
  std::string genes;

  po::options_description desc("Clownfish options - writes to STDOUT");
  desc.add_options()
    ("help,h", "Show help message")
    ("kmer-length,k", po::value<int>(&kmer_length)->required(), "Kmer length to use")
    ("hash-size,s", po::value<unsigned long>(&hash_size)->required(), "Initial hash table size for Jellyfish (will grow if full)")
    ("reads,r", po::value<std::string>(&reads)->required(), "Reads from sample, fast(a|q)")
    ("genes,g", po::value<std::string>(&genes)->required(), "Genes file, fasta")
    ("threads,t", po::value<int>(&num_threads)->default_value(1), "Number of threads to use")
    ("kmer-stats", po::bool_switch(&kmer_stats)->default_value(false), "Print max count, total, and unique kmers from Jellyfish (Adds large time cost)");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
      std::cout << desc;
      exit(0);
    } 
    po::notify(vm);
  } catch(po::error &e) {
    std::cerr << "ERROR:" << e.what() << std::endl << desc;
    exit(1);
  }

  std::cerr << "Starting clownfish with: kmer_length=" << kmer_length << 
                                         ",hash_size=" << hash_size <<
                                             ",reads=" << reads <<
                                             ",genes=" << genes <<
                                           ",threads=" << num_threads <<
                                         ",kmer_stats=" << kmer_stats << std::endl;

  mer_dna::k(kmer_length);
  mer_hash hash(hash_size, mer_dna::k() * 2, counter_length, num_threads, max_reprobe);

  std::ifstream genes_fs(genes);
  if(!genes_fs.good()) {
      std::cerr << "Problem opening '" << genes << "' - will now exit" << std::endl;
      exit(1);
  }
  std::cout.precision(20);

  /* ---------- Fill Jellyfish Hashtable ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Filling jellyfish hash ===" << std::endl;
    // jellyfish likes c style strings 'n stuff
    char **reads_c = new char*[1];
    reads_c[0] = (char *) reads.c_str();

    mer_counter counter(num_threads, hash, reads_c, reads_c + 1, canonical);
    counter.exec_join(num_threads);

    delete[] reads_c;
  }

  mer_array* ary = hash.ary();

  /* ---------- Calculate kmer stats ---------- */
  if(kmer_stats) {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    std::cerr << "=== Calcuating kmer stats ===" << std::endl;
    uint64_t max = 0;
    uint64_t total = 0;
    uint64_t distinct = 0;
    uint64_t val = 0;    
    for(auto it = ary->begin(); it != ary->end(); ++it) {
      val = it->second;
      distinct++;
      total += val;
      max = std::max(val, max);
    }
    std::cerr << "distinct:" << distinct << std::endl;
    std::cerr << "max:" << max << std::endl;
    std::cerr << "total:" << total << std::endl;
  }

  /* ---------- Count kmers in genes ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    std::cerr << "=== Counting genes ===" << std::endl;
    
    mer_dna mer;
    std::string buf; 
    // skip the first line
    std::getline(genes_fs, buf); 
    while(!genes_fs.eof()) {
      // "parse" the fasta file, wouldnt be necessary if fasta files didn't allow wrapped lines!
      // TODO: move this into a boost filter or something
      std::string gene = "";
      do {
        std::getline(genes_fs, buf);
        // 2 checks on the same condition, don't know how to combine
        if(buf[0] != '>') {
          gene += buf; 
        }
      } while(buf[0] != '>' && !genes_fs.eof());

     uint64_t sum = 0;
     uint64_t val;
     int length = gene.size() - kmer_length + 1;
      for(int i = 0; i < length; ++i) {
        mer = gene.substr(i, kmer_length);
        if(!ary->get_val_for_key(mer, &val)) {
          val = 0;
        }
        sum += val;
      }
      std::cout << (double) sum / length << std::endl;
    }
  }

  //std::cerr << "=== Exiting clownfish count ===" << std::endl;
  return 0;
}
