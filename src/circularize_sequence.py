#!/usr/bin/env python3

from Bio import SeqIO
import argparse


def circularize_sequence(input_file):
    '''
    Takes the output file from megahit and keep the contigs larger
    than 100 Kbp.
    Contigs sequence is fixed to match polyhedrin ATG codon as start
    position.
    Additionally, we add a missing sequence from hrs5 that was
    confirmed using SPAdes assembler and Tablet for visualisation
    of reads alignment at the region 50170.
    '''

    polyhedrin = 'ATGTATACTCGTTACAGCTATAACCCATCTTTGGGTCGCACCTACGTGTA'
    hr5 = 'ACAATCTTTGCTTTCGGTGAAGTGTTTCGCTGAAAGCAAACTTTGATAAAATGACGCAATAAAATGATAAAATTATTGTGCAATAAAGTCTTCAATGTTTGCTTTCGGCAAAGTGTTTCGCTGAAAGCAAAGATTGCGATTATTGCACAATGA'
    sequences = SeqIO.parse(input_file, 'fasta')
    for sequence in sequences:
        if len(sequence) > 100000 and polyhedrin in sequence.seq:
            genome_seq = sequence.seq
            start = genome_seq.find(polyhedrin)
            region1 = genome_seq[start:]
            region2 = genome_seq[:start]
            sequence.seq = region1 + hr5 + region2
        elif len(sequence) > 100000:
            genome_seq = sequence.seq.reverse_complement()
            start = genome_seq.find(polyhedrin)
            region1 = genome_seq[start:]
            region2 = genome_seq[:start]
            sequence.seq = region1 + hr5 + region2
    print(f'>genome_draft\n{sequence.seq}')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='File with the contigs from megahit')
    args = parser.parse_args()
    return args.input


def main():
    input_file = argument_parser()
    circularize_sequence(input_file)


if __name__ == '__main__':
    main()
