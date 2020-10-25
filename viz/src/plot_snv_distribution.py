#!/usr/bin/env python3

import argparse
import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_variants_frequency(data):
    '''
    Takes a VCF file and return a dictionary
    containing position, frequency and annotation
    for each variant.
    Additionally, discard rare variants affecting
    start or stop codons (there are only 2 in the
    dataset).
    '''

    # Read variants from VCF file
    variants = list(vcf.Reader(open(data, 'r')))
    # Set variants classes to discard
    discard = ['splice_region_variant&stop_retained_variant',
               'stop_gained',
               'start_lost']
    # Store variants positions
    positions = [int(rec.POS) for rec in variants
                 if rec.INFO['ANN'][0].split('|')[1] not in discard]
    # Store variants frequencies
    freqs = [float(rec.INFO['AF']) for rec in variants
             if rec.INFO['ANN'][0].split('|')[1] not in discard]
    # Store variants annotations
    annotations = [rec.INFO['ANN'][0].split('|')[1]
                   for rec in variants
                   if rec.INFO['ANN'][0].split('|')[1] not in discard]
    # Create a dictionary with filtered variants
    filtered_variants = {'Position': positions,
                         'Frequency': freqs,
                         'Annotation': annotations}

    return filtered_variants


def plot_snv_distribution(data):
    '''

    '''

    # Convert data to a pandas DataFrame
    df = pd.DataFrame(data)
    # Set default settings for Seaborn
    sns.set()
    # Set figure size
    plt.figure(figsize=(12, 12))
    # Rename variant classes
    df = df.replace(
            {'missense_variant': 'Nonsynonymous',
             'synonymous_variant': 'Synonymous',
             'upstream_gene_variant': 'Intergenic'})
    # Set colors for variant classes
    color_mapping = {'Nonsynonymous': '#beaed4',
                     'Synonymous': '#ffff99',
                     'Intergenic': '#fdc086'}
    # Set variant classes order
    order_list = ['Synonymous',
                  'Nonsynonymous',
                  'Intergenic']
    # Create figure panels (axis)
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2, 2), (1, 0))
    ax3 = plt.subplot2grid((2, 2), (1, 1))

    # Scatterplot of variants position vs frequency
    sns.scatterplot(x='Position',
                    y='Frequency',
                    data=df,
                    hue='Annotation',
                    palette=color_mapping,
                    edgecolor='black',
                    legend=False,
                    ax=ax1)
    ax1.set(xlabel='Genome Position',
            ylabel='Allele Frequency')

    # Barplot with the number of variants per class
    sns.countplot(x='Annotation',
                  data=df,
                  palette=color_mapping,
                  edgecolor='gray',
                  linewidth=1,
                  order=order_list,
                  ax=ax2)
    ax2.set(xlabel='',
            ylabel='SNVs')

    # Boxplot with variant frequency per class
    sns.boxplot(x='Annotation',
                y='Frequency',
                data=df,
                palette=color_mapping,
                linewidth=1,
                order=order_list,
                ax=ax3)
    ax3.set(xlabel='',
            ylabel='Allele Frequency')

    # Save figure
    plt.savefig('snv_distribution.svg')


def main():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('variants',
                        help='VCF variants file')
    args = parser.parse_args()
    variants = args.variants
    results = get_variants_frequency(variants)
    plot_snv_distribution(results)


if __name__ == '__main__':
    main()
