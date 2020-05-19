#!/usr/bin/env python3

from Bio import AlignIO
import vcf
import argparse
import csv

def extract_variant(vcf_file):
    '''
    Extract variants position for SfMNPV from a .vcf file.
    '''
    dic = {}
    positions = []
    references = []
    variants = []
    vcf_handle = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_handle:
        for item in record:
            if item.sample == 'sfmnpv_argentina' and item.gt_type:
                positions.append(record.POS)    
                references.append(record.REF)
                variants.append(str(record.ALT))
    dic['positions'] = positions
    dic['references'] = references
    dic['Variants'] = variants
    return dic

def correct_positions(align_file, variants):
    '''
    Correct variants positions extracted from a sequence alignment so that
    they coincide with sfmnpv_argentina positions.
    '''
    # Load sequence alignment
    align = AlignIO.read(align_file, 'fasta')
    for seq in align:
        if seq.id == 'sfmnpv_argentina':
            seq = seq.seq
    # Load positions from dictionary created from extract_variant()
    positions = variants['positions'] 

    # For each variants search gapped positions in the alignment and
    # substract it from the position number
    positions_updated = []
    for pos in positions:
        initial = int(pos)
        gaps = 0
        for index, base in enumerate(seq):
            if index < pos and base == '-':
                gaps += 1
        final = initial - gaps
        positions_updated.append(final)

    # Update variants dictionary
    variants['positions'] = positions_updated
    # Write outputs to a .csv file
    print(variants.keys())
    with open('sfmnpv_argentina_variants.csv', 'w') as fh:
        writer = csv.writer(fh)
        writer.writerows(variants.keys())
        writer.writerows(variants.items())
    
def arguments_parser():
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf',
                        help='VCF file created by snp-sites')
    parser.add_argument('align',
                        help='Fasta sequence alignment')
    args = parser.parse_args()
    return args.vcf, args.align

def main():
    vcf_file, align_file = arguments_parser()
    variants = extract_variant(vcf_file)
    correct_positions(align_file, variants)

if __name__ == '__main__':
    main()
