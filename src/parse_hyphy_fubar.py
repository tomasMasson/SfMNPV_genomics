#!/usr/bin/env python3

from Bio import AlignIO
import json
import pandas as pd
import argparse


def parse_hyphy(hyphy_file):
    '''
    Takes an hyphy output file and returns two lists,
    for the data and headers.
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
    Takes a data and headers lists from parse_hyphy
    and filter the positions corresponding to the Human sequence.
    Returns a pandas dataframe'''

    algn = AlignIO.read(algn_file, 'fasta')
    for sequence in algn:
        if 'YP_001036321.1' in sequence.id:
            residues = [item for item in sequence.seq]

    df['Omega'] = df['beta'] / df['alpha']
    df['Residue'] = residues
    df = df[df['Residue'] != '-']
    df.index = range(len(df))
    return df


def assign_selective_pressure(df):
    '''
    Takes a datafame with hyphy results and assign a column
    displaying the type of selective pressure in each residue.
    '''

    selection_type = []
    for index, row in df.iterrows():
        if row['Prob[alpha>beta]'] >= 0.9:
            selection_type.append('Purifiyng')
        elif row['Prob[alpha<beta]'] >= 0.9:
            selection_type.append('Diversifiyng')
        else:
            selection_type.append('Neutral')
    df['Selection_type'] = selection_type
    return df


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('hyphy',
                        help='Hyphy results in .json format')
    parser.add_argument('alignment',
                        help='Alignment used to filter codons corresponding to human sequence')
    args = parser.parse_args()
    return args.hyphy, args.alignment


def main():
    hyphy_file, alignment = argument_parser()
    df = parse_hyphy(hyphy_file)
    df = filter_dataframe(df, alignment)
    df = assign_selective_pressure(df)
    output_file = f'{hyphy_file.split(".json")[0]}.csv'
    df.to_csv(output_file, index=False)


if __name__ == '__main__':
    main()
