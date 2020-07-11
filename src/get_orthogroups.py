#!/usr/bin/env python3

import argparse
from Bio import SearchIO
from Bio import SeqIO


def get_orthogroups(blast_xml, multifasta):
    '''
    This script takes a blast output file in xml format and
    returns a multifasta for each orthogroup (group of sequences
    matching against the same query).
    Sequences are extracted from a multifasta comprising all
    proteins present in all SfMNPV isolates.
    '''

    blast_handle = SearchIO.parse(blast_xml, 'blast-xml')
    orthogroups = {}
    for record in blast_handle:
        key = record.id.split('_')[0]
        orthogroups[key] = [hit.id for hit in record]
    sequences = list(SeqIO.parse(multifasta, 'fasta'))
    for key in orthogroups:
        with open(f'{key}.faa', 'w') as f:
            for sequence in sequences:
                if sequence.id in orthogroups[key]:
                    f.write(f'>{sequence.id}\n{sequence.seq}\n')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('blast',
                        help='input BLAST results file')
    parser.add_argument('multifasta',
                        help='Fasta with proteome sequences')
    args = parser.parse_args()
    return args.blast, args.multifasta


def main():
    blast, proteome = argument_parser()
    get_orthogroups(blast, proteome)


if __name__ == '__main__':
    main()
