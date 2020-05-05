#!/usr/bin/env python3

from Bio import AlignIO
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def sequence_similarity(seq1, seq2):
    '''
    Calculates the pairwise sequence identity between two sequences.
    '''
    length = len(seq1)
    diff = 0
    for i, j in zip(seq1, seq2):
        if i != j:
            diff += 1
    similarity = ((length - diff) * 100) / length
    return similarity

def get_similarity_profile(seq1, seq2, step: int, window: int):
    '''
    Computes the sequence similarity score for two sequences using a
    sliding window approach with a custom step and window. 
    To calculate similarity, this function uses sequences_similarity().
    '''
    length = len(seq1)
    profile = []
    for index in range(0, length, step):
        subseq1 = seq1[index:index+window]
        subseq2 = seq2[index:index+window]
        profile.append(sequence_similarity(subseq1, subseq2))

    return profile

def alignment_similarity_profile(align_file, step: int, window: int):
    '''
    Takes a multiple sequence alignment and extract the individual
    sequences. Also set a reference sequence. 
    '''
    align = AlignIO.read(align_file, 'fasta')
    sequences = []
    for seq in align:
        if seq.id == 'sfmnpv_argentina':
            reference = seq
        else:
            sequences.append(seq)
    profiles = {}
    for seq in sequences:
        name = seq.id
        profiles[name] = get_similarity_profile(seq.seq, reference.seq, step, window)
    profiles['positions'] = np.arange(0, len(reference), step) 
    return profiles

def plot_similarity(profiles):
    '''Plot sequences pairs similarities. '''

    df = pd.DataFrame(profiles)
    sns.lineplot(x='positions', y='HM595733.1', data=df)
    sns.lineplot(x='positions', y='EU258200.1', data=df)
    sns.lineplot(x='positions', y='JF899325.1', data=df)
    sns.lineplot(x='positions', y='NC_009011.2', data=df)
    sns.lineplot(x='positions', y='KF891883.1', data=df)
    plt.legend(['Nicaraguan', '19', 'Defective', '3AP2', 'Colombian'])
    plt.xlabel('Genome Alignment Position')
    plt.ylabel('Sequence Similarity (%)')
    plt.savefig('sequence_similarity_plot.svg')

def argument_parser():
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('align',
                        help='Sequence alignment file in Fasta format')
    parser.add_argument('step', type=int, help='Step size')
    parser.add_argument('width', type=int, help='Windows size')
    args = parser.parse_args()
    return args.align, args.step, args.width

def main():
    align, step, width = argument_parser()
    profiles = (alignment_similarity_profile(align, step, width))
    plot_similarity(profiles)

if __name__ == '__main__':
    main()
