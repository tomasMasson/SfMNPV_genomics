#!/usr/bin/env python3

from Bio import SeqIO
import argparse


def rearrange_sequence(input_file):
    '''
    Takes the output file from megahit and keep the contigs larger
    than 10 Kbp.
    Contigs sequence is fixed to match the granulin gene as start
    position.
    '''

    granulin = 'ATGTATACTCGTTACAGCTATAACCCATCTTTGGGTCGCACCTACGTGTA'
    sequences = SeqIO.parse(input_file, 'fasta')
    for sequence in sequences:
        if len(sequence) > 100000:
            genome_seq = sequence.seq.reverse_complement()
            start = genome_seq.find(granulin)
            region1 = genome_seq[start:]
            region2 = genome_seq[:start]
            sequence.seq = region1 + region2
            print(f'>sfmnpv_argentina\n{sequence.seq}')

            
def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='File with the contigs from megahit')
    args = parser.parse_args()
    return args.input


def main():
    input_file = argument_parser()
    rearrange_sequence(input_file)


if __name__ == '__main__':
    main()
