#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>

#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>

#include "JellyfishHash.hpp"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {

  bool kmer_stats = false;
  int num_threads = 1;

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
    ("debug,d", po::value<bool>(&kmer_stats)->default_value(false), "Print max count, total, and unique kmers from Jellyfish");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
      std::cout << desc;
      exit(0);
    } 
    po::notify(vm);
  } catch(po::error &e) {
    std::cerr << "ERROR:" << e.what() << std::endl;
    exit(1);
  }

  mer_dna::k(kmer_length);
  jelly_hash hash(hash_size, kmer_length, num_threads, reads);

  std::ifstream genes_fs(genes);
  std::cout.precision(20);

  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Filling jellyfish hash ===" << std::endl;
    hash.fill();
  }

  if(kmer_stats) {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    std::cerr << "=== Calcuating kmer stats ===" << std::endl;
    uint64_t distinct;
    uint64_t max;
    uint64_t total;
    hash.stats(distinct, max, total);
    std::cerr << "distinct:" << distinct << std::endl;
    std::cerr << "max:" << max << std::endl;
    std::cerr << "total:" << total << std::endl;
  }

  {
    boost::timer::auto_cpu_timer t(std::cerr, 2); 
    std::cerr << "=== Counting genes ===" << std::endl;
    
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
     int length = gene.size() - kmer_length + 1;
      for(int i = 0; i < length; ++i) {
        std::string kmer = gene.substr(i, kmer_length);
        sum += hash.get(kmer);
        //std::cerr << val << ' ';
        //output_fs.write(reinterpret_cast<const char *>(&val), sizeof(uint32_t));
      }
      std::cout << (double) sum / length << std::endl;
    }
  }

  std::cerr << "=== Exiting clownfish count ===" << std::endl;
  return 0;
}
