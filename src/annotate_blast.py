#!/usr/bin/env python3

import argparse
from Bio import SearchIO
from Bio import Entrez
from Bio import SeqIO


def parse_blast_result(blast_xml):
    '''
    This script takes a blast output file in xml format and
    returns a table displaying the query id, the best hit id
    and the annotation retrieved from Entrez database
    '''

    # Results keep protein homologue from SfMNPV 3AP2, because it
    # has a standarize identifier (YP_XXX or NP_XXX)
    blast_handle = SearchIO.parse(blast_xml, 'blast-xml')
    results = {}
    for record in blast_handle:
        hits = [hit.id for hit in record]
        for hit in hits:
            if 'NC_009011' in hit:
                pseudoname = hit.split('prot_')[1]
                hit_name = pseudoname.split('.')[0]
                results[record.id] = hit_name
                break
            else:
                pseudoname = hit.split('prot_')[1]
                hit_name = pseudoname.split('.')[0]
                results[record.id] = hit_name
    return results


def fetch_annotations(blast_results):
    '''
    Fetch Entrez annotations for each hit in a BLAST search.
    '''

    # BLAST hits annotation retrieval from NCBI Entrez using
    # the protein identifier from SfMNPV 3AP2
    Entrez.email = 'tomas.masson@biol.unlp.edu.ar'
    queries = [query for query in blast_results.values()]
    entrez_search = Entrez.efetch('protein',
                                  id=queries,
                                  rettype='gb',
                                  retmode='text')
    search_handle = SeqIO.parse(entrez_search, 'gb')
    queries_annotation = {}
    for record in search_handle:
        queries_annotation[record.name] = record.description

    # Create annotation using BLAST records (proteins) and the
    # annotations retrieved from NCBI
    annotation = {}
    for record in blast_results:
        for query in queries_annotation:
            if blast_results[record] == query:
                annotation[record] = [query, queries_annotation[query]]
    return annotation


def get_features_table(raw_annotation):
    '''
    Convert a raw annotation into a feature table similar
    to .gtf format.
    '''

    # Extract gene coordinates from fasta header and creates an
    # unsorted feature table
    unsorted_annotation = {}
    for number, record in enumerate(raw_annotation):
        name = number
        start = int(record.split(':')[1]) + 1
        end = int(record.split(':')[2]) + 1
        annotation = ' '.join(word for word in raw_annotation[record])
        if end > start:
            strand = '+'
            unsorted_annotation[name] = [start, end, strand, annotation]
        else:
            strand = '-'
            unsorted_annotation[name] = [end, start, strand, annotation]
    # Sort feature table
    sorted_table = {key: value for key, value in sorted(unsorted_annotation.items(), key=lambda item: item[1])}

    # Create a dictionary with common protein names
#    with open(protein_names, 'r') as f:
#        names = {line.split()[0]: line.split()[1] for line in f}

    # Print to standard output
    for index, key in enumerate(sorted_table, 1):
        data = sorted_table[key]
        print(f'genome_assembly\tBLASTp\tCDS\t{data[0]}\t{data[1]}\t.\t{data[2]}\t0\t{data[3]}')


def argument_parser():
    '''Command line argument parser.'''

    parser = argparse.ArgumentParser()
    parser.add_argument('blast',
                        help='input BLAST results file')
    args = parser.parse_args()
    return args.blast


def main():
    raw_blast = argument_parser()
    blast_results = parse_blast_result(raw_blast)
    raw_annotations = fetch_annotations(blast_results)
    get_features_table(raw_annotations)


if __name__ == '__main__':
    main()
