#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import MutableSeq
import vcf
import argparse


def get_alternative_alleles(vcf_file):
    '''
    Takes a .vcf file and returns a dictionary showing
    the sites were alternative allele frequency is greater
    than reference allele.
    This output can be used to mutate the reference sequence.
    '''

    vcf_handle = vcf.Reader(open(vcf_file, 'r'))
    alternative_alleles = {}
    for record in vcf_handle:
        if record.INFO['AF'] >= 0.5:
            alternative_alleles[record.POS] = [record.REF,
                                               str(record.ALT[0]),
                                               record.INFO['AF']]
    return alternative_alleles


def mutate_sequence(sequence, mutations):
    '''
    Uses a dictionary (created by get_alternative_alleles)
    to select positions in a reference sequence and introduce
    specific point mutation.
    '''

    seq = SeqIO.read(sequence, 'fasta')
    mutable_seq = seq.seq.tomutable()
    for key in mutations:
        mutable_seq[key-1] = mutations[key][1]
    new_seq = mutable_seq.toseq()
    print(f'>genome_assembly\n{new_seq}')


def argument_parser():
    ''' Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence',
                        help='Faste file to be mutated')
    parser.add_argument('vcf',
                        help='.vcf file containing the mutations')
    args = parser.parse_args()
    return args.sequence, args.vcf


def main():
    sequence, vcf = argument_parser()
    mutations = get_alternative_alleles(vcf)
    mutate_sequence(sequence, mutations)


if __name__ == '__main__':
    main()
