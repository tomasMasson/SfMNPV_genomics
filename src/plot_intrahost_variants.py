#!/usr/bin/env python3

import argparse
import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def variants_frequency(variants):
    '''
    Takes a VCF file and return a dictionary containing position,
    frequency and annotation information.
    Additionally, discard rare variantes affecting start or stop
    codons (there are only 2 in the dataset).
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
    variants_dic = {'Position': positions,
                    'Frequency': frequencies,
                    'Annotation': annotations}
    return variants_dic


def scatterplot_variants(data):
    '''
    Creates a scatter plot of variants frequency.
    '''

    df = pd.DataFrame(data)
    sns.set()
    plt.figure(figsize=(12, 6))
    color_mapping = {'missense_variant': '#beaed4',
                     'synonymous_variant': '#ffff99',
                     'upstream_gene_variant': '#fdc086'}
    sns.scatterplot(x='Position',
                    y='Frequency',
                    data=df,
                    hue='Annotation',
                    palette=color_mapping,
                    edgecolor='black',
                    legend=False)
    plt.xlabel('Genome Position')
    plt.ylabel('Frequency')
    plt.savefig('snv_distribution.svg')
    plt.close()


def frequency_plus_counts_variants(data):
    '''
    Creates a boxplot showing the frequency of each variant class
    (synonymous, missense and non_coding) at the right panel.
    In the left panel creates a barplot showing the couny for each class
    '''
    df = pd.DataFrame(data)
    df = df.replace({'missense_variant': 'Missense',
                     'synonymous_variant': 'Synonymous',
                     'upstream_gene_variant': 'Intergenic'})
    color_mapping = {'Missense': '#beaed4',
                     'Synonymous': '#ffff99',
                     'Intergenic': '#fdc086'}

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    sns.set()
    order_list = ['Synonymous', 'Missense', 'Intergenic']

    sns.countplot(x='Annotation',
                  data=df,
                  palette=color_mapping,
                  edgecolor='gray',
                  linewidth=1,
                  order=order_list,
                  ax=axs[0])
    axs[0].set(xlabel='', ylabel='Count')

    sns.boxplot(x='Annotation',
                y='Frequency',
                data=df,
                palette=color_mapping,
                linewidth=1,
                order=order_list,
                ax=axs[1])
    axs[1].set(xlabel='', ylabel='Frequency')

    plt.savefig('snv_classification.svg')


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
