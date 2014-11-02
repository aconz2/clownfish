#include <string>
#include <fstream>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>

using jellyfish::mer_dna;
typedef jellyfish::stream_manager<char**> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > read_parser;

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
    for(mer_iterator mers(parser_, canonical_) ; mers; ++mers)
      ary_.add(*mers, 1);
    ary_.done();
  }
};

int main(int argc, char *argv[]) {
  if(argc != 6) {
    std::cerr << "Error: Wrong number of arguments\n"
              << "Usage: " << argv[0] << " <k_mer_len> <canonical> <hash_size> <nb_threads> <reads_file> < <genes_file>" << std::endl;
    exit(1);
  }

  const int kmer_length = std::stoi(argv[1]);
  mer_dna::k(kmer_length);
  const bool canonical = strcmp(argv[2], "true") == 0;
  const size_t initial_hash_size = std::stoul(argv[3]);
  const int    nb_threads        = std::stoi(argv[4]);

  std::cerr << "Starting to fill jelly hash" << std::endl;
  // Parse the input file.
  // 7 is the length in bits of the counter field
  // 126 is the maximum reprobe value
  mer_hash ary(initial_hash_size, mer_dna::k() * 2, 7, nb_threads, 126);
  {
    // Count the k-mers in the reads and store them in the hash array 'ary'
    mer_counter counter(nb_threads, ary, argv + 5, argv + 6, canonical);
    counter.exec_join(nb_threads);
  }

  std::cerr << "Done filling jelly fish" << std::endl; 
  std::cerr << "Starting to count genes" << std::endl;
  mer_dna mer;
  uint64_t val = 0;
  auto hash = ary.ary();
 
  for(std::string line; std::getline(std::cin, line); ) {
    for(int i = 0; i < (int) line.size() - kmer_length + 1; ++i) {
      mer = line.substr(i, kmer_length);
      if(!hash->get_val_for_key(mer, &val)) {
        val = 0;
      }
      std::cout << val << ' ';
    }
    std::cout << '\n';
  }

  // Query the database for all the mers on the command line
/* commented out for test
  mer_dna  m;
  uint64_t val = 0;
  auto hash = ary.ary();
  for(int i = 5; i < argc; ++i) {
    m = argv[i];
    m.canonicalize();
    if(!hash->get_val_for_key(m, &val))
      val = 0;
    std::cout << m << ' ' << val << '\n';
  }
*/

  return 0;
}
