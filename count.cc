#include <string>
#include <fstream>

#include <boost/timer/timer.hpp>

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
  if(argc != 8) {
    std::cerr << "Error: Wrong number of arguments\n"
              << "Usage: " << argv[0] << " <k_mer_len> <canonical> <hash_size> <nb_threads> <reads_file> <genes_file> <output_file>" << std::endl;
    exit(1);
  }

  const int kmer_length = std::stoi(argv[1]);
  mer_dna::k(kmer_length);
  const bool canonical = strcmp(argv[2], "true") == 0;
  const size_t initial_hash_size = std::stoul(argv[3]);
  const int    nb_threads        = std::stoi(argv[4]);
  std::ifstream genes_fs(argv[6]);
  std::ofstream output_fs(argv[7], std::ios::binary);

  std::cerr << "=== Filling jellyfish hash ===" << std::endl;
  // Parse the input file.
  // 7 is the length in bits of the counter field
  // 126 is the maximum reprobe value
  mer_hash ary(initial_hash_size, mer_dna::k() * 2, 7, nb_threads, 126);
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    // Count the k-mers in the reads and store them in the hash array 'ary'
    mer_counter counter(nb_threads, ary, argv + 5, argv + 6, canonical);
    counter.exec_join(nb_threads);
  }

  auto hash = ary.ary();

  {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    std::cerr << "=== Calcuating kmer stats ===" << std::endl;
    uint64_t max = 0;
    uint64_t total = 0;
    uint64_t distinct = 0;
    for(auto it = hash->begin(); it != hash->end(); ++it) {
      auto pair = *it;
      //std::cerr << pair.first << ':'  << pair.second <<  std::endl;
      distinct++;
      total += pair.second;
      if(pair.second > max) {
        max = pair.second;
      }
    }
    std::cerr << "distinct:" << distinct << std::endl;
    std::cerr << "max:" << max << std::endl;
    std::cerr << "total:" << total << std::endl;
  }

  {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    std::cerr << "=== Counting genes ===" << std::endl;

    mer_dna mer;
    uint64_t val = 0;
    uint64_t val_ = 0;
    std::string buf; 
    std::string gene;
    
    // skip the first line
    std::getline(genes_fs, buf); 
    while(!genes_fs.eof()) {
      // "parse" the fasta file, wouldnt be necessary if fasta files 
      // didn't allow wrapped lines!
      gene = "";
      do {
        std::getline(genes_fs, buf);
        // 2 checks on the same condition, don't know how to combine
        if(buf[0] != '>') {
          gene += buf; 
        }
      } while(buf[0] != '>' && !genes_fs.eof());

      // count the kmers
      unsigned int length = (int) gene.size() - kmer_length + 1;
      output_fs.write(reinterpret_cast<const char *>(&length), sizeof(uint32_t));
      for(unsigned int i = 0; i < length; ++i) {
        mer = gene.substr(i, kmer_length);
        if(!hash->get_val_for_key(mer, &val)) {
          val = 0;
        }
        // when canonical flag, set the value to be the max of it and its reverse complement
        if(canonical) {
          mer.reverse_complement();
          if(!hash->get_val_for_key(mer, &val_)) {
            val_ = 0;
          }
          if(val_ > val) {
            val = val_;
          }
        } 

        //std::cerr << val << ' ';
        output_fs.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
      }
      //std::cerr << std::endl;
    }
  }

  return 0;
}
