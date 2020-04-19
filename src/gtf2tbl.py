#!/usr/bin/env python3

import argparse

def gtf2tbl(gtf_file):
    '''
    Takes a .gff annotation file and returns an NCBI feature table file.
    '''

    with open(gtf_file, 'r') as fh:
        print(f'>Feature sfmnpv_argentina')
        for line in fh:
            data = line.strip().split('\t')
            if data[6] == '+':
                cds = f'{data[3]}\t{data[4]}\tCDS\n\t\t\tproduct\t{data[8]}'
                print(f'{cds}')
            else:
                cds = f'{data[4]}\t{data[3]}\tCDS\n\t\t\tproduct\t{data[8]}'
                print(f'{cds}')

def main():
    ''' Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf',
                        help='.gtf file containing the annotations')
    args = parser.parse_args()
    gtf = args.gtf
    gtf2tbl(gtf)    

if __name__ == '__main__':
    main()
