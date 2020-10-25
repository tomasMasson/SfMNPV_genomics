#!/usr/bin/env python3

import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def plot_snp_distributions(data):
    '''
    Creates a histogram displaying the
    distribution of SNP across genome position
    and barplot for SNPs between isolates.
    '''

    # Read SNPs from VCF file
    snps = list(vcf.Reader(open(data, 'r')))
    # Create a list with the positions of the SNPs
    positions = [snp.start for snp in snps]
    # Group SNPs according to its isolate
    snp_groups = {}
    for record in snps:
        # Extract isolate name
        isolate = [field['GT'] for field in record.samples]
        isolate = ''.join(isolate)
        # Add isolate if its the first time
        if isolate not in snp_groups:
            snp_groups[isolate] = 1
        # Increase SNP count
        elif isolate in snp_groups:
            count = snp_groups[isolate] + 1
            snp_groups[isolate] = count

    # Set defaults for Seaborn
    sns.set()
    # Create DataFrame from isolates SNPs
    df = pd.DataFrame.from_dict(snp_groups,
                                orient='index',
                                columns=['SNPs'])
    # Rename isolates names
    isolates_names = {
            '011111': '3AP2',
            '000001': 'ARG-M',
            '001010': '19+Colombian',
            '001000': '19',
            '000010': 'Colombian',
            '010101': 'ARG-M+NicB+NicG',
            '010100': 'NicB+NicG',
            '010111': '3AP2+19',
            }
    df = df.rename(isolates_names)

    # Create two vertical figures
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,12))
    # Plot histogram of SNPs distribution
    sns.histplot(data=positions,
                 bins=13,
                 color='grey',
                 ax=ax1)
    ax1.set(xlabel='Genome Position',
            ylabel='SNPs')
    # Set isolates order for barplot
    order = ['Colombian', '19', '3AP2', 'NicB+NicG',
             'ARG-M+NicB+NicG', 'ARG-M',
             '19+Colombian', '3AP2+19']
    # Plot the number of SNPs per isolate
    sns.barplot(data=df,
                x='SNPs',
                y=df.index,
                order=order,
                color='gray',
                ax=ax2)

    # Save figure
    fig.savefig('snp_distribution.svg')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('vcf',
                        help='VCF file created with SNP-sites')
    args = parser.parse_args()
    return args.vcf


def main():
    data = argument_parser()
    plot_snp_distributions(data)


if __name__ == '__main__':
    main()
