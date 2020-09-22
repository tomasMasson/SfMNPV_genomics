#!/usr/bin/env python3


def create_ft(annotation):
    """
    Annotation model:

    CDS     1..741
    gene="polh"
    note="ORF1; Similar to SfMNPV-3AP2 sf1, AcMNPV-C6 ac8"
    codon_start=1
    product="Polyhedrin"
    """
    # Source header definitions
    organism = "Spodoptera frugiperda multiple nucleopolyhedrovirus"
    mol_type = "genomic DNA"
    isolate = "ARG"
    db_xref = "taxon:10455"
    country = "Argentina"
    note = "SfMNPV-ARG"
    collection_date = ""
    collected_by = ""

    print(f'>Features {note}')

    with open(annotation, "r") as fh:
        next(fh)
        for line in fh:
            fields = line.strip().split(',')
            polarity = fields[3]
            if polarity == "+":
                start = fields[1]
                stop = fields[2]
            else:
                start = fields[2]
                stop = fields[1]
            gene = fields[4]
            product = fields[5]
            if fields[7] != '-':
                note = f"{fields[0]}; Similar to SfMNPV-3AP2 {fields[6]}, AcMNPV-C6 {fields[7]}"
            else:
                note = f"{fields[0]}; Similar to SfMNPV-3AP2 {fields[6]}"
            coordinates = f'{start}\t{stop}\tCDS'
            gene_rec = f'\t\t\tgene\t{gene}'
            codon_rec = f'\t\t\tcodon_start\t1'
            prod_rec = f'\t\t\tproduct\t{product}'
            note_rec = f'\t\t\tnote\t{note}'
            print(coordinates)
            print(gene_rec)
            print(codon_rec)
            print(prod_rec)
            print(note_rec)


create_ft("annotation_table.csv")
