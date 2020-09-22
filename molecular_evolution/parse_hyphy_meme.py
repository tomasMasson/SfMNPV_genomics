#!/usr/bin/env python3

from Bio import AlignIO
import csv
import json
import pandas as pd
import argparse


def parse_hyphy(hyphy_file):
    '''
    Takes an hyphy output file and returns two lists, for the data and headers.
    '''
    with open(hyphy_file, 'r') as fh:
        dic = json.load(fh)
    data = dic['MLE']['content']['0']
    headers = [header[0] for header in dic['MLE']['headers']]
    df = pd.DataFrame(data)
    df = df.loc[:, (df != 0).any(axis=0)]
    df.columns = headers
    return df


def filter_dataframe(df, algn_file):
    '''
    Takes a data and headers lists from parse_hyphy and filter the positions
    corresponding to the Human sequence. Returns a pandas dataframe.
    '''
    algn = AlignIO.read(algn_file, 'fasta')
    for sequence in algn:
        if 'YP_001036321.1' in sequence.id:
            residues = [item for item in sequence.seq]

    df['Residue'] = residues
    df = df[df['Residue'] != '-']
    df.index = range(1, len(df)+1, 1)
    return df


def assign_selective_pressure(df):
    '''
    Takes a datafame with hyphy results and assign a column displaying
    the type of selective pressure in each residue.
    '''
    episodic_selection = []
    for index, row in df.iterrows():
        if row['p-value'] <= 0.05:
            episodic_selection.append('Episodic Diversifiyng')
        else:
            episodic_selection.append('Not significant')
    df['Episodic Selection'] = episodic_selection
    return df


def main():
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('hyphy',
                        help='File with hyphy results in .json format')
    parser.add_argument('Alignment',
                        help='Sequence alignment used to filter codons corresponding to human sequence')
    args = parser.parse_args()
    hyphy_file, alignment = args.hyphy, args.Alignment
    output_file = f'{hyphy_file.split(".json")[0]}.csv'
    df = parse_hyphy(hyphy_file)
    df = filter_dataframe(df, alignment)
    df = assign_selective_pressure(df)
    df.to_csv(output_file, index=True)


if __name__ == '__main__':
    main()
