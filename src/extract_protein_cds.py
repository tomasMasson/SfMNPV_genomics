#!/usr/bin/env python3

'''
Extract all protein coding sequences from a genome sequence in
Fasta format using as guide an annotation table in GTF format.
'''

from Bio import SeqIO
import argparse


def extract_protein_cds(genome, annotation):
    '''
    Extract protein coding sequences from a genome in Fasta format.
    Coding regions are specified using an annotation table in GTF
    format.
    '''
    sequence = SeqIO.read(genome, 'fasta')
    with open(annotation, 'r') as f:
        for index, line in enumerate(f):
            data = line.split()
            gene = data[10]
            start = int(data[3]) - 1
            end = int(data[4])
            seq = sequence.seq
            polarity = data[6]
            if polarity == '+':
                cds = seq[start:end-3].translate()
                print(f'>cds{index+1}_{gene}\n{cds}\n')
            else:
                cds = seq[start+3:end].reverse_complement().translate()
                print(f'>cds{index+1}_{gene}\n{cds}\n')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('Genome',
                        help='Genome sequence in Fasta format')
    parser.add_argument('Annotation',
                        help='Annotation file with coding seqs')
    args = parser.parse_args()
    return args.Genome, args.Annotation


def main():
    genome, annotation = argument_parser()
    extract_protein_cds(genome, annotation)


if __name__ == '__main__':
    main()
