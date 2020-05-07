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
    discard = ['splice_region_variant&stop_retained_variant',
               'stop_gained',
               'start_lost']
    positions = [int(record.POS) for record in vcf_handle
                 if record.INFO['ANN'][0].split('|')[1] not in discard]
    frequencies = [float(record.INFO['AF']) for record in vcf_handle
                   if record.INFO['ANN'][0].split('|')[1] not in discard]
    # iList to discard rare events (only 4 stop codon gain/loss)
    annotations = [record.INFO['ANN'][0].split('|')[1]
                   for record in vcf_handle
                   if record.INFO['ANN'][0].split('|')[1] not in discard]
    variants_dic = {'Position':positions, 
                    'Frequency': frequencies,
                    'Annotation': annotations}
    return variants_dic

def scatterplot_variants(data):
    '''
    Creates a scatter plot of the number of SNV across genome windows
    '''
    
    df = pd.DataFrame(data) 
    sns.set_style('ticks')
    plt.figure(figsize=(18, 6))
    sns.scatterplot(x='Position',
                    y='Frequency',
                    data=df,
                    hue='Annotation',
                    palette='deep',
                    x_jitter=True,
                    y_jitter=True)
    plt.xlabel('Genome Position')
    plt.ylabel('SNV Frequency')
    plt.savefig('SNV_frequency_distribution.svg')
    plt.close()

def frequency_plus_counts_variants(data):
    '''
    Creates a boxplot showing the frequency of each variant class 
    (synonymous, missense and non_coding) at the right panel.
    In the left panel creates a barplot showing the couny for each class
    '''
    fig, axs = plt.subplots(1,2, figsize=(10,6)) 
    df = pd.DataFrame(data) 
    sns.set_style('ticks')
    sns.boxplot(x='Annotation',
                y='Frequency',
                data=df,
                palette='deep',
                ax=axs[1])

    classes = ['Synonymous', 'Missense', 'Non Coding']
    counts = df['Annotation'].value_counts()
    sns.barplot(x=classes, y=counts, ax=axs[0],
                palette='deep', linewidth=1.5, edgecolor='.2')
    plt.savefig('SNV_count.svg')

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
    scatterplot_variants(results)
    frequency_plus_counts_variants(results)

if __name__ == '__main__':
    main()
