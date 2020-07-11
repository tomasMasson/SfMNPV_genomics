#!/usr/bin/env python3

import argparse


def create_snpeff_gtf(gtf_file):
    '''
    Add gene products names to a GTF file.
    '''

    with open(gtf_file, 'r') as fh:
        for index, line in enumerate(fh):
            record = line.split()
            fld1 = record[0]
            fld2 = record[1]
            fld3 = 'CDS'
            fld4 = record[3]
            fld5 = record[4]
            fld6 = '.'
            fld7 = record[6]
            fld8 = '0'
            fld9 = f'gene_id "cds{index+1}"'
            print(f'{fld1}\t{fld2}\t{fld3}\t{fld4}\t{fld5}\t{fld6}\t{fld7}\t{fld8}\t{fld9};')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='GTF file')
    args = parser.parse_args()
    return args.gtf


def main():
    gft = argument_parser()
    create_snpeff_gtf(gft)


if __name__ == '__main__':
    main()
