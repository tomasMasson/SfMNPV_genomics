#!/usr/bin/env python3

import argparse


def add_names_gtf(gtf_file, names_files):
    '''
    Add gene products names to a GTF file.
    '''

    with open(names_files, 'r') as fh:
        names = [(record.split()[0], record.split()[1])
                 for record in fh]

    with open(gtf_file, 'r') as fh:
        for line in fh:
            record = line.split()
            fld1 = record[0]
            fld2 = record[1]
            fld3 = 'CDS'
            fld4 = record[3]
            fld5 = record[4]
            fld6 = '.'
            fld7 = record[6]
            fld8 = '0'
            for item in names:
                if item[0] == record[10]:
                    fld9 = f'gene_id "{item[1]}"'
            print(f'{fld1}\t{fld2}\t{fld3}\t{fld4}\t{fld5}\t{fld6}\t{fld7}\t{fld8}\t{fld9};')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='GTF file')
    parser.add_argument('names', help='Names file')
    args = parser.parse_args()
    return args.gtf, args.names


def main():
    gft, names = argument_parser()
    add_names_gtf(gft, names)


if __name__ == '__main__':
    main()
