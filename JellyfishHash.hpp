#include <string>

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

class jelly_hash {
  bool canonical_ = false;
  int max_reprobe_ = 126;
  int counter_length_ = 7; // length in bits
  int kmer_length_;
  int num_threads_;
  std::string file_;
  mer_hash hash_;
  mer_array* ary_;

public:
  jelly_hash(unsigned long hash_size, int kmer_length, int num_threads, std::string& file) :
    kmer_length_(kmer_length),
    num_threads_(num_threads),
    file_(file),
    hash_(hash_size, mer_dna::k(kmer_length) * 2, counter_length_, max_reprobe_)
    {
      ary_ = hash_.ary();
    }

  void fill() {
    // pointer hacks to make jellyfish happy
    int length = file_.length() + 1;
    char *file_cpy = new char [length];
    memcpy(file_cpy, file_.c_str(), sizeof(char) * length);
    char *file_cpy_end = file_cpy + sizeof(char) * length;
    //std::cerr << file_begin << std::endl;
    mer_counter counter(num_threads_, hash_, &file_cpy, &file_cpy_end, canonical_);
    counter.exec_join(num_threads_);
  }

  uint64_t get(const std::string& kmer_s) {
    // TODO: probably use a mer_iterator instead of passing all these strings
    // These two statements not joinable
    mer_dna kmer;
    kmer = kmer_s;
    uint64_t val = 0;

    if(!ary_->get_val_for_key(kmer, &val)) {
      val = 0;
    }

    return val; 
  }

  void stats(uint64_t& distinct, uint64_t& max, uint64_t& total) {
    uint64_t max_ = 0;
    uint64_t total_ = 0;
    uint64_t distinct_ = 0;
    uint64_t val = 0;    
    for(auto it = ary_->begin(); it != ary_->end(); ++it) {
      val = it->second;
      distinct_++;
      total_ += val;
      max_ = std::max(val, max);
    }
    distinct = distinct_;
    max = max_;
    total = total_;
  }
 
};
