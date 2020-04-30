#!/usr/bin/env python3

from Bio import AlignIO
import argparse
import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def count_variants_per_window(alignment, variants, width, step):
    '''
    Takes a multiple sequence alignment in .fasta format and a .vcf file
    created by snp-sites.
    Return a .csv files containing the number of SNV per region in the
    genome. The width and step of the window are provided as arguments.
    '''

    align = AlignIO.read(alignment, 'fasta')
    vcf_handle = vcf.Reader(open(variants, 'r'))
    variants_record = [int(record.POS) for record in vcf_handle]
    variants_per_window = {}
    for i in range(0, len(align[0]) - int(width) + 1, int(step)):
        count = 0
        for variant in variants_record:     
            if variant >= i and variant < (i+int(width)):
                count += 1
        variants_per_window[i] = count
    
    return variants_per_window

def plot_variants(data):
    '''
    Creates a scatter plot of the number of SNV across genome windows
    '''
    
    df = pd.DataFrame(data, index=[0])
    df = df.T
    df.columns = ['#SNV']
    plot = sns.scatterplot(x=df.index, y=df['#SNV'], data=df)
    plot.figure.savefig('isolates_diversity/SNV_distribution.png')

def arguments_parser():
    '''Command line argument parser.'''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment', help='Alignment file (.fasta)')
    parser.add_argument('variants', help='Variants file (.vcf)')
    parser.add_argument('width', help='Window size (int)')
    parser.add_argument('step', help='Step size (int)')
    args = parser.parse_args()

    align = args.alignment
    variants = args.variants
    width = args.width
    step = args.step
    return align, variants, width, step

def main():
    align, variants, width, step = arguments_parser() 
    results = count_variants_per_window(align, variants, width, step)
    plot_variants(results)

if __name__ == '__main__':
    main()
