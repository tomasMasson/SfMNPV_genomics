#!/usr/bin/env python3

import vcf
import argparse

def get_mutation_gene(mutations, annotations):
    '''
    Retrieves the gene affected by each mutation in a VCF file.
    '''
    vcf_reader = vcf.Reader(open(mutations, 'r'))
    mut_list = [int(record.POS) for record in vcf_reader]

    with open(annotations, 'r') as fh:
        genes = [(record.split()[3],
                  record.split()[4],
                  record.strip().split('\t')[8])
                  for record in fh]
    for mut in mut_list:
        for gene in genes:
            if mut > int(gene[0]) and mut < int(gene[1]):
                print(f'{mut} affect {gene[2]}')
            else:
                print(f'{mut} affect {gene[2]}')

def argument_parser(): 
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('mut', help='Mutations VCF file')
    parser.add_argument('ann', help='Mutations VCF file')
    args = parser.parse_args()
    return args.mut, args.ann

def main():
    mutations, annotations = argument_parser()
    get_mutation_gene(mutations, annotations)

if __name__ == '__main__':
    main()
