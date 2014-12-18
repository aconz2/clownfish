# Clownfish #
Clownfish is a gene abundance estimator. Given a set of reference genes and a
set of reads, we wish to assign an abundance to each gene. Clownfish works at
the k-mer level and (currently) calculates the abundance to be the mean of the
abundance of its component k-mers.

Clownfish is under active development and new methods are currently be developed and evaluated.

# Requirements #

## C libraries ##
- jellyfish >= 2.0
- boost = 1.54 (may work with earlier versions, not tested)
- gcc >= 4.8.3 (may work with earlier versions, not tested)

## Python Modules ##
- python3 = 3.3.2 (may work with other versions, not tested)
- Biopython
- pandas
- numpy

# Installation #
There is no autoconf or cmake script at this time.

```
git clone https://github.com/aconz2/clownfish.git
make
export PATH="`pwd`/bin:$PATH"
```

Make sure that jellyfish is in your `PKG_CONFIG_PATH` and `LD_LIBRARY_PATH`, something like
- `PKG_CONFIG_PATH` /some/prefix/dir/jellyfish-2.1.4/build/lib/pkgconfig/ 
- `LD_LIBRARY_PATH` /some/prefix/dir/jellyfish-2.1.4/build/lib/ 

And that `BOOST_ROOT` contains the root to your boost installation, something like
- `BOOST_ROOT` /some/prefix/dir/boost-1.54.0

# Usage #

## clownfish ##

```
Clownfish options - writes to STDOUT:
  -h [ --help ]             Show help message
  -k [ --kmer-length ] arg  Kmer length to use
  -s [ --hash-size ] arg    Initial hash table size for Jellyfish (will grow if
                            full)
  -r [ --reads ] arg        Reads from sample, fast(a|q)
  -g [ --genes ] arg        Genes file, fasta
  -t [ --threads ] arg (=1) Number of threads to use
  --counter-length arg (=7) Length (in bits) of counter field in Jellyfish hash
  --kmer-stats              Print max count, total, and unique kmers from 
                            Jellyfish (Adds large time cost)
  --raw                     Print tab-seperated counts for each gene, one per 
                            line. ie. don't take the mean
```

- `--hash-size` is passed to jellyfish to presize its hash table. This serves
  as an estimate and will grow itself when necessary.
- `--counter-length` is passed to jellyfish. This is usually not necessary
  unless using a smaller k, usually < 14. The problem arises with low k because
the kmers appear much more frequently and quickly use the default 7 bits to
count with. Jellyfish "solves" this by then stealing another entry in the hash
table to use as a counter. This quickly fills the size of the table, even with
a large hash size.
- `--reads` and `--genes` should be uncompressed. You may pass a process
  substitution to uncompress on the fly. Ex: `<(zcat genes.fasta.gz)`
- `--kmer-stats` will print the number of unique and distinct kmers from the
  reads, as well as the max-count. This does add an appreciable amount of time.
- `--raw` the output on each line is a list of numbers, tab separated. The
  `j`th number on the `i`th line will be the frequencty in which the `j`th kmer
from the `i`th gene appears in the set of reads.

`clownfish` writes to STDOUT and logs to STDERR.

## cf_merge ## 

```
usage: cf_merge [-h] --genes GENES cf_files [cf_files ...]

positional arguments:
  cf_files       Files to merge

optional arguments:
  -h, --help     show this help message and exit
  --genes GENES  Genes file
```

Once clownfish has been run on one or more set of reads and you have captured
STDOUT to say, `/some/prefix/dir/sample1.cf` and `/some/prefix/dir/sample2.cf`
you can run `cf_merge` to combine these into a matrix with the gene id as row
labels and sample id as column label. Ex.

```
cf_merge --genes <same gene file used in counting step> /some/prefix/dir/sample1.cf /some/prefix/dir/sample2.cf
```

This will give the labels of `sample1` and `sample2` for columns 1 and 2.

`cf_merge` writes to STDOUT in a tab seperated format
