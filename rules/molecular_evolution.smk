rule all:
    input:
        "sf29.fubar.json"

#rule blast_search:
#    input:
#        "sequence.faa"
#    output:
#        "blast_search.xml"
#    params:
#        db="database",
#        taxids="baculoviridae.txids"
#    shell:
#        """
#        blastp -query {input} -db {params.db} -taxids {params.taxids} -outfmt 5
#        """

rule retrieve_sequences:
    input:
        "{gene}_blast.xml"
    output:
        "{gene}.cds.fna",
        "{gene}.cds.faa"
    shell:
        """
        python src/retrieve_cds.py {input}
        """

rule protein_alignment:
    input:
        "{gene}.cds.faa"
    output:
        "{gene}.aln.faa"
    shell:
        """
        mafft --localpair --maxiterate 1000 {input} >> {output}
        """

rule codon_alignment:
    input:
        "{gene}.aln.faa",
        "{gene}.cds.fna"
    output:
        "{gene}.aln.fna"
    shell:
        """
        pal2nal.pl {input} -output fasta >> {output}
        """

rule phylogeny_inference:
    input:
        "{gene}.aln.faa"
    output:
        "{gene}.treefile"
    params:
        "{gene}"
    shell:
        """
        iqtree -s {input} --prefix {params} -bb 1000 
        """

#rule gard_recombination_test:
#    input:
#        "{gene}.aln.fna"
#    output:
#        "{gene}.gard.json"
#    shell:
#        """
#        hyphy gard --alignment {input} --output {output} 
#        """
#
#rule busted_selection_test:
#    input:
#        "{gene}.aln.fna",
#        "{gene}.treefile"
#    output:
#        "{gene}.busted.json"
#    shell:
#        """
#        hyphy busted --alignment {input[0]} --tree {input[1]} --output {output}
#        """
        
rule fubar_selection_test:
    input:
        "{gene}.aln.fna",
        "{gene}.treefile"
    output:
        "{gene}.fubar.json"
    shell:
        """
        hyphy fubar --alignment {input[0]} --tree {input[1]} --output {output}
        """
