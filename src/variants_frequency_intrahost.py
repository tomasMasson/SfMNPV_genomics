#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def variants_frequency(variants):
    '''
    Takes a multiple sequence alignment in .fasta format and a .vcf file
    created by snp-sites.
    Return a .csv files containing the number of SNV per region in the
    genome. The width and step of the window are provided as arguments.
    '''

    vcf_handle = list(vcf.Reader(open(variants, 'r')))
    variants_record = [int(record.POS) for record in vcf_handle]
    variants_freq = [float(record.INFO['AF']) for record in vcf_handle]
    variants_dic = {'Position': variants_record,
                           'Frequency': variants_freq}
    return variants_dic

def plot_variants(data):
    '''
    Creates a scatter plot of the number of SNV across genome windows
    '''
    
    df = pd.DataFrame(data) 
    plt.figure(figsize=(8,4))
    ax = sns.scatterplot(x='Position', y='Frequency', data=df,
                         hue='Frequency', legend=False, 
                         x_jitter=True, y_jitter=True)
    ax.set(xlabel='Genome Position',
           ylabel='Frequency',
           title='SfMNPV Genomic Diversity')
#    plt.show()
    ax.figure.savefig('SNV_frequence_distribution.svg')

def arguments_parser():
    '''Command line argument parser.'''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('variants', help='Variants file (.vcf)')
    args = parser.parse_args()

    variants = args.variants
    return variants

def main():
    variants = arguments_parser() 
    results = variants_frequency(variants)
    plot_variants(results)

if __name__ == '__main__':
    main()
