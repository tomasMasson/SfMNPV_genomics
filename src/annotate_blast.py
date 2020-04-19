#!/usr/bin/env python3

import argparse
from Bio import SearchIO
from Bio import Entrez
from Bio import SeqIO

def parse_blast_result(blast_results):
    '''
    This script takes a blast output file in xml format and returns
    a table displaying the query id, the best hit id and the annotation
    retrieved from Entrez database
    '''

    blast_handle = SearchIO.parse(blast_results, 'blast-xml')
    results = {}
    for record in blast_handle:
        for hit in record:
            results[record.id] = f'YP_{hit.id.split("_")[4][:-2]}'
    return results

def fetch_annotations(results):
    '''Fetch Entrez annotations for each hit in a BLAST search.'''

    # BLAST hits annotation retrieval from NCBI Entrez
    Entrez.email = 'tomas.masson@biol.unlp.edu.ar'
    entrez_query = [query for query in results.values()]
    entrez_search = Entrez.efetch('protein',
                                  id=entrez_query,
                                  rettype='gb', retmode='text')
    entrez_handle = SeqIO.parse(entrez_search, 'gb')
    hits_annotation = {}
    for record in entrez_handle:
        hits_annotation[record.name] = record.description 

    # Query identifier and annotation matching 
    queries_annotation = {}
    for record in results:
        for hit in hits_annotation:
            if results[record] == hit:
                queries_annotation[record] = [hit, hits_annotation[hit]]

    return queries_annotation

def get_features_table(annotation_file):
    '''
    Convert a raw annotation into a feature table similar to .gtf format.
    '''     

    unsorted_table = {}
    for number, record in enumerate(annotation_file):
        name  = number
        start = int(record.split(':')[1]) + 1
        end = int(record.split(':')[2]) + 1
        annotation = ' '.join(word for word in annotation_file[record])
        if end > start:
            strand = '+'
            unsorted_table[name] = [start, end, strand, annotation]
        else:
            strand = '-'
            unsorted_table[name] = [end, start, strand, annotation]
    sorted_table = {key: value for key, value in sorted(unsorted_table.items(), key=lambda item: item[1])}

    for index, key in enumerate(sorted_table, 1):
        data = sorted_table[key]
#        print(f'{index}\t{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\n')
        print(f'sfmnpv\tBLASTp\tCDS\t{data[0]}\t{data[1]}\t.\t{data[2]}\t0\tSimilar to {data[3]}')

def argument_parser():
    '''Execute blast_parser main function.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='<BLAST results in xml format>',
                        type=str, help='input BLAST results file')
    args = parser.parse_args()
    blast_results = args.input

    return blast_results

def main():
    raw_blast = argument_parser()
    blast_results = parse_blast_result(raw_blast)
    raw_annotations = fetch_annotations(blast_results)
    features_table = get_features_table(raw_annotations)

if __name__ == '__main__':
    main()
