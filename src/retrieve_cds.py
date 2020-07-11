#!/usr/bin/env python3

import argparse
from Bio import Entrez
from Bio import SearchIO


def retrieve_cds(blast_xml):
    '''
    Retrieve coding sequences from for the significant hits of a
    BLAST search output file (in XML format).
    '''

    blast_handle = SearchIO.read(blast_xml, 'blast-xml')
    identifiers = [hit.id for hit in blast_handle]

    Entrez.email = 'tomasmasson0@gmail.com'
    search_handle = Entrez.epost(
            db='protein',
            id=','.join(identifiers),
            usehistory='y'
            )
    results_handle = Entrez.read(search_handle)
    search_handle.close()
    webenv = results_handle['WebEnv']
    query_key = results_handle['QueryKey']

    fetch_handle = Entrez.efetch(
            db='protein',
            rettype='fasta_cds_na',
            retmode='text',
            webenv=webenv,
            query_key=query_key
            )
    print(fetch_handle.read())


def argument_parser():
    '''Command line argument parser.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('blast_xml',
                        help='BLAST results file in .xml')
    args = parser.parse_args()
    return args.blast_xml


def main():
    blast_xml = argument_parser()
    retrieve_cds(blast_xml)


if __name__ == '__main__':
    main()