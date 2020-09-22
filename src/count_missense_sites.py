#!/usr/bin/env python3

import argparse
import numpy as np
from Bio import SeqIO


def count_segregating_sites(alignment):
    '''
    This script takes an alignment in fasta format and
    returns the number of polymorphic sites.
    '''

    align = list(SeqIO.parse(alignment, format='fasta'))
    gene = alignment.split('/')[2].split('.')[0]
    align_array = np.array([list(rec)
                           for rec in align])
    sites = 0
    for site in align_array.T:
        if len(set(site)) > 1 and '-' not in site:
            sites += 1
    return f'{gene},{sites}'


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('alignment',
                        help='Fasta with proteome sequences')
    args = parser.parse_args()
    return args.alignment


def main():
    alignment = argument_parser()
    print(count_segregating_sites(alignment))


if __name__ == '__main__':
    main()
