#!/usr/bin/env python3

import argparse
import vcf


def make_snv_report(vcf_file):
    '''
    Makes a .csv report from a vcf annotated file.
    Report includes information about Position, Reference,
    Alternative, Frequence, Gene and Mutation of the variant.
    '''
    print('Position,Reference,Alternative,Frequence,Type,Gene,AA_variant')
    vcf_data = vcf.Reader(open(vcf_file))
    for record in vcf_data:
        # Extract data fields
        position = record.POS
        reference = record.REF
        variant = record.ALT[0]
        frequency = record.INFO['AF']
        variant_type = record.INFO['ANN'][0].split('|')[1]
        gene = '-'
        protein_variant = '-'
        if variant_type != 'upstream_gene_variant':
            gene = record.INFO['ANN'][0].split('|')[4]
            protein_variant = record.INFO['ANN'][0].split('|')[10]
        # Print to standard output
        print(f'{position},{reference},{variant},{frequency},{variant_type},{gene},{protein_variant}')


def command():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', help='VCF annotated file')
    args = parser.parse_args()
    return args.vcf


def main():
    vcf_file = command()
    make_snv_report(vcf_file)


if __name__ == '__main__':
    main()
