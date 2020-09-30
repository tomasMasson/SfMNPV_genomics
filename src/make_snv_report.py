#!/usr/bin/env python3

import argparse
import vcf


def get_gene_class(gene):
    """
    Return the class of a gene.
    Available options are Core or Non-core.
    """
    core_genes = ['pif-0', 'pif-1', 'pif-2',
                  'pif-3', 'pif-4', 'pif-5',
                  'pif-6', 'pif-7', 'pif-8',
                  'dnapol', 'desmoplakin',
                  'dnahel', 'alk-exo', 'p33',
                  'lef-1', 'lef-2', 'lef-4',
                  'lef-5', 'p18',
                  'lef-8', 'lef-9', 'p47',
                  'vlf-1', '38k', 'ac53-like',
                  'p40', 'p6.9', 'vp1054',
                  'vp39', 'ac81-like', 'gp41',
                  'odv-e18', 'odv-e25', 'p45',
                  'odv-ec43', 'p49',
                  'ac78-like', 'odv-ec27'
                  ]
    if gene in core_genes:
        return "Core gene"
    else:
        return "Non-core gene"


def change_polarity(mutation):
    """
    Return + if the protein mutation introduce
    a change in the polarity of the residue.
    """

    groups = {
            "Ala": "Non polar",
            "Asn": "Polar",
            "Asp": "Acidic",
            "Arg": "Basic",
            "His": "Basic",
            "Cys": "Cysteine",
            "Gln": "Polar",
            "Glu": "Acidic",
            "Gly": "Glicine",
            "Pro": "Proline",
            "Leu": "Non polar",
            "Lys": "Basic",
            "Ile": "Non polar",
            "Met": "Non polar",
            "Phe": "Non polar",
            "Ser": "Polar",
            "Thr": "Polar",
            "Tyr": "Polar",
            "Trp": "Non polar",
            "Val": "Non polar",



            }
    if "Ter" in mutation:
        return "-"
    if "*" in mutation:
        return "-"
    if "?" in mutation:
        return "-"

    polarity_wt = groups[mutation[2:5]]
    polarity_mutant = groups[mutation[-3:]]
    if polarity_wt != polarity_mutant:
        return '+'
    else:
        return '-'


def make_snv_report(vcf_file):
    '''
    Makes a .csv report from a vcf annotated file.
    Report includes information about Position, Reference,
    Alternative, Frequence, Gene and Mutation of the variant.
    '''
    print('Position,Reference,Alternative,Frequence,Variant_class,Gene,Gene_class,AA_variant,polarity_change')
    vcf_data = vcf.Reader(open(vcf_file))
    for record in vcf_data:
        # Extract data fields
        position = record.POS
        reference = record.REF
        variant = record.ALT[0]
        frequency = record.INFO['AF']
        variant_class = record.INFO['ANN'][0].split('|')[1]
        gene = '-'
        gene_class = '-'
        aa_variant = '-'
        polarity_change = '-'
        if variant_class != 'upstream_gene_variant':
            gene = record.INFO['ANN'][0].split('|')[4]
            gene_class = get_gene_class(gene) 
            aa_variant = record.INFO['ANN'][0].split('|')[10]
            polarity_change = change_polarity(aa_variant)
        # Print to standard output
        print(f'{position},{reference},{variant},{frequency},{variant_class},{gene},{gene_class},{aa_variant},{polarity_change}')


def command():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', help='VCF annotated file')
    args = parser.parse_args()
    return args.vcf


def main():
    vcf_file = command()
    make_snv_report(vcf_file)


if __name__ == '__main__':
    main()
