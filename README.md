# Notes #
    - IGC.fa has 10M genes
    - Metahit has 3.9M
    - mhgc longest gene is 88086
    - igc longest gene is 88230
    - using jellyfish.QueryMerDNA was convenient, but very slow. I believe it does a binary search of the file, and (maybe) mmaps it. only got through 500K genes in 3 hours. the jellyfish counting only took 2 mins of that! Now trying to load it into the python dict (but still reading the jellyfish database binary file directly)
    - loading into dict was real slow! and used more memory than jellyfish did
    - now have a cpp program using jellyfish as a lib to do the counting, then querying for each gene, outputs to stdout a space seperated list of counts

# TODO:
clownfish count does not handle resizing of jellyfish hash well i.e. memory explodes. everything is fine so long as the --size given is safely over the number of distinct kmers. may want to move to taking a .jf file as input or using a .jf file as an intermediate.

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


