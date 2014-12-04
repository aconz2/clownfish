# Notes #
    - IGC.fa has 10M genes
    - Metahit has 3.9M
    - mhgc longest gene is 88086
    - igc longest gene is 88230

# TODO:
  - move the per-gene scoring into c++, each gene score will always be the mean, so do that to write to the file
  - clownfish count does not handle resizing of jellyfish hash well i.e. memory explodes. everything is fine so long as the --size given is safely over the number of distinct kmers. may want to move to taking a .jf file as input or using a .jf file as an intermediate.

  - 20 simulated samples, abundances random from 0-20, readlength 70, no mistakes, 30M reads, at least 4 samples must contain each gene for the canopy to work
  - run to kmer size 35
  - install canopy

# jellyfish -C flag test
input: ATCGCGGTA

without -C
TCG 1
ATC 1
CGG 1
GTA 1
GCG 1
GGT 1
CGC 1

with -C
CGA 1
ATC 1
GTA 1
ACC 1
CGC 2
CCG 1

so CGC collapsed with GCG to become 2 and not 1 each


