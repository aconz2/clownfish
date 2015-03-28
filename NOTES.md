# Notes #
    - IGC.fa has 10M genes
    - Metahit has 3.9M
    - mhgc longest gene is 88086
    - igc longest gene is 88230

# TODO:

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


