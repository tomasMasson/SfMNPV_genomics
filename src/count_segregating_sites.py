#!/usr/bin/env python3

import argparse
import numpy as np
from Bio import SeqIO


def count_segregating_sites(alignment):
    '''
    This script takes a blast output file in xml format and
    returns a multifasta for each orthogroup (group of sequences
    matching against the same query).
    Sequences are extracted from a multifasta comprising all
    proteins present in all SfMNPV isolates.
    '''

    align = list(SeqIO.parse(alignment, format='fasta'))
    gene = "".join(align[0].id.split('_')[1:])
    align_array = np.array([list(rec)
                           for rec in align])
    sites = 0
    for site in align_array.T:
        if len(set(site)) > 1 and '-' not in site:
            sites += 1
    return f'{sites},{gene}'


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
