#!/usr/bin/env python3

import vcf
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def filter_snv_sites(filename):
    '''
    Takes a VCF file created by snp-sites and return a
    dictionary ready to be inyected into pandas.
    Variants are displayed only for mutant sites.
    '''

    vcf_reader = vcf.Reader(open(filename, 'r'))
    positions = []
    isolate = []
    for record in vcf_reader:
        for item in record:
            if item.gt_type != 0:
                positions.append(record.POS)
                isolate.append(item.sample)
    data = {'Genome Position': positions,
            'Isolate': isolate}
    return data


def plot_snv_distribution(data):
    '''
    Creates a categorical plot displaying the distribution
    of SNV for each isolate.
    '''

    df = pd.DataFrame(data)
    sns.stripplot(x='Genome Position', y='Isolate', data=df)
    fig = plt.gcf()
    fig.set_size_inches(16, 8)
    fig.savefig('isolates_snv_distribution.svg')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('vcf',
                        help='VCF file created with SNP-sites')
    args = parser.parse_args()
    return args.vcf


def main():
    filename = argument_parser()
    data = filter_snv_sites(filename)
    plot_snv_distribution(data)


if __name__ == '__main__':
    main()
