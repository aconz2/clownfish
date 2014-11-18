#!/usr/bin/env python3

"""Count occurrences of each kmer in each gene against a sample"""

import os
import re
import sys
import logging
import argparse
import resource
import functools
import subprocess
import struct

from Bio import SeqIO
import pandas
import numpy

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

def memory_usage():
    """Returns in kilobytes, the total (or maybe peak) memory usage"""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def file(f):
    if os.path.exists(f):
        return f
    else:
        logging.error('file "%s" does not exist, cannot continue', f)
        sys.exit(1)

def trim_and_clamp(arr,
                   left_trim=None,
                   right_trim=None,
                   minimum=None,
                   maximum=None):
    """arr should be a numpy array,
    the first `left_trim` values in the array will be dropped
    the last `right_trim` values in the array will be dropped
    values in the array < `minimum` will be set to 0
    values in the array > `maximum` will be set to 0
    (uses kwargs for easy partial application, all parameters
    may be given as None)"""
    if left_trim is not None:
        arr = arr[left_trim:]

    if right_trim is not None:
        arr = arr[0:-right_trim]

    if minimum is not None:
        mins = numpy.ones(arr.size) * minimum
        arr = numpy.multiply(arr, arr >= mins, arr)

    if maximum is not None:
        maxs = numpy.ones(arr.size) * maximum
        arr = numpy.multiply(arr, arr <= maxs, arr)

    return arr

def parse_cf(handle):
    """Yields each numpy.array of integers found in the cf binary file
    format is `((uint32 length)(uint32,..., length of them))...` and
    uint32 is little endian"""
    # this is scary!
    while True:
        length_b = handle.read(4)
        if length_b == b'':
            return 
        # convert uint32 little endian to integer
        (length, ) = struct.unpack('<I', length_b)
        yield numpy.fromstring(handle.read(4 * length), dtype=numpy.uint32)

def score_cf(handle, method,
                     left_trim=None, 
                     right_trim=None,
                     maximum=None,
                     minimum=None):
    """Score a cf file with parameters. Return a float list"""
    trim_and_clamper = functools.partial(trim_and_clamp,
                                         left_trim=left_trim,
                                         right_trim=right_trim,
                                         minimum=minimum,
                                         maximum=maximum)
    arrs = parse_cf(handle)
    return [method(trim_and_clamper(arr)) for arr in arrs]

def score_cfs(cf_files, gene_labels, method,
                                     left_trim=None, 
                                     right_trim=None,
                                     maximum=None,
                                     minimum=None):
    """Convenience method for scoring multiple cf_files. Returns a pandas
    DataFrame. kwargs are the same as `score_cf`"""
    samples = {} 
    for cf_file in cf_files:
        with open(cf_file, 'rb') as cf_handle:
            logging.debug('Scoring %s...', cf_file)
            scores = score_cf(cf_handle, method,
                                         left_trim=left_trim, 
                                         right_trim=right_trim,
                                         maximum=maximum,
                                         minimum=minimum)

            samples[re.sub('\.cf$', '', cf_file.split('/')[-1])] = scores

    return pandas.DataFrame(data=samples, index=gene_labels)
 
def count_main():
    """Count kmers in `reads` file, then list the corresponding count for 
    each gene in `genes`"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genes', required=True, type=file,
                        help='FASTA file with genes')
    parser.add_argument('--samples', required=True,
                        type=file, help='FASTA or FASTQ samples')
    parser.add_argument('--canonical', action='store_true', default=False,
                        help='Consider each count to be the max of it and it\'s reverse complement')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to run jellyfish counter with')
    parser.add_argument('--size', type=str, required=True,
                        help='Initial size of jellyfish hash table') 
    parser.add_argument('--kmer-length', type=int, required=True, 
                        help='Length of kmer to use')
    parser.add_argument('--kmer-stats', action='store_true', default=False,
                        help='Output distinct, max count, and total kmers from sample')
    parser.add_argument('--output', type=str, default='occur.cf',
                        help='Name of output file to write to')
    args = parser.parse_args(sys.argv[2:])
    logging.debug(args)

    if os.path.exists(args.output):
        logging.error('Output file "%s" already exists, will not overwrite it', args.output)
        sys.exit(1)

    # lower case it!
    args.canonical = 'true' if args.canonical else 'false'
    args.kmer_stats = 'true' if args.kmer_stats else 'false'
    # <k_mer_len> <canonical> <hash_size> <nb_threads> <reads_file> <genes_file> <output>
    child_args = [os.path.join(os.path.dirname(__file__), 'count'),
                 args.kmer_length, args.canonical,  args.size, args.threads,
                 args.samples, args.genes, args.output, args.kmer_stats]
    child_args = list(map(str, child_args))

    logging.debug("Launching child process `count`")
    logging.debug(child_args)
    subprocess.check_call(child_args,
                          stdout=sys.stdout, stderr=sys.stderr,
                          close_fds=False)

def score_main():
    """Scores one or more cf files. This collapses each sample to be a
    vector of scores, one for each gene. The output will have each sample
    as a column and each gene as a row. Can be output to csv or a pickled
    pandas data_frame"""
    score_methods = {'mean': numpy.mean,
                     'median': numpy.median}

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--left-trim', type=int, default=None,
                        help='Drop the first <left> numbers')
    parser.add_argument('--right-trim', type=int, default=None,
                        help='Drop the last <right> numbers')
    parser.add_argument('--min', type=int, default=None,
                        help='Consider values < <min> to be 0')
    parser.add_argument('--max', type=int, default=None,
                        help='Consider values > <max> to be 0')
    parser.add_argument('--method', choices=list(score_methods.keys()),
                        required=True, type=str,
                        help='Scoring method to use for kmer counts')
    parser.add_argument('cf_file', nargs='+', type=file,
                        help='Clownfish binary file to score')
    parser.add_argument('--output', required=True, type=str,
                        help='File to write to')
    parser.add_argument('--genes', required=True, type=file,
                        help='Genes counted, used for labels')
    parser.add_argument('--format', default='pickle', choices=['pickle', 'tsv'],
                        help='Output format')
    args = parser.parse_args(sys.argv[2:])
    logging.debug(args)

    if os.path.exists(args.output):
        print('File `{}` already exists, will not overwrite'.format(args.output))
        sys.exit(1)

    scoring_method = score_methods[args.method]

    logging.debug('Reading gene labels')
    with open(args.genes) as genes_handle:
        gene_labels = [record.id for record in SeqIO.parse(genes_handle, 'fasta')]

    data_frame = score_cfs(args.cf_file, gene_labels,
                                         scoring_method,
                                         left_trim=args.left_trim,
                                         right_trim=args.right_trim,
                                         minimum=args.min,
                                         maximum=args.max)

    logging.debug('Starting to write output file')
    if args.format == 'pickle':
        data_frame.to_pickle(args.output)
    else:
        data_frame.to_csv(args.output, sep='\t')

def dump_main():
    """Convert the cf binary file to ascii"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cf_file', type=file, help='Clownfish binary file to dump')
    args = parser.parse_args(sys.argv[2:])
    logging.debug(args)
   
    try: 
        with open(args.cf_file, 'rb') as cf_handle:
            arrs = parse_cf(cf_handle)
            str_arrs = ('\t'.join(map(str, arr)) for arr in arrs)
            sys.stdout.writelines(str_arr + '\n' for str_arr in str_arrs)
    except BrokenPipeError:
        pass

def main():
    """Estimate relative gene abundance using kmers"""
    dir(file)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['count', 'score', 'dump'])
    args = parser.parse_args(sys.argv[1:2])
    if args.command == 'count':
        count_main()
    elif args.command == 'score':
        score_main()
    elif args.command == 'dump':
        dump_main()
    else:
        print("Invalid choice")
        sys.exit(1)
    
if __name__ == "__main__":
    main()
