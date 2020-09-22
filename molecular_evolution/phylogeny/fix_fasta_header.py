#!/usr/bin/env python3

import argparse

def fix_fasta_header(filename):
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line.lstrip('>').split()[0]
                specie_name = line.split('[')[1][:-2]
                name = '_'.join(specie_name.split())
                print(f'>{protein_id}_{name}')
            else:
                print(line.strip())


def command():
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    arg = parser.parse_args()
    return arg.file


if __name__ == '__main__':
    fix_fasta_header(command())
