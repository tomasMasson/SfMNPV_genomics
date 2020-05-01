#!/usr/bin/env python3

import argparse
import pysam

def get_softclipped_reads(bam_file):
    '''
    Filter a .bam file and returns the reads with soft-clipped regions
    bigger than 20 nucleotides.
    The output file is in BAM format.
    '''
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    clipped_reads = pysam.AlignmentFile('soft-clipped_reads.bam',
                                        'wb', template=samfile)

    for read in samfile:
        cigar = read.cigarstring
        if cigar != '101M' and cigar is not None:
            clips = [item for item in read.cigartuples
                     if item[0] == 4 and item[1] > 20]
            if clips:
                clipped_reads.write(read)

def argument_parser():
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='BAM file')
    args = parser.parse_args()
    bam_file = args.bam
    return bam_file

def main():
    bam_file = argument_parser()
    get_softclipped_reads(bam_file)

if __name__ == '__main__':
    main()
