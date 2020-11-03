rule retrieve_sequences:
    input:
        "data/{gene}_blast.xml"
    output:
        "molecular_evolution/{gene}.cds.fna",
        "molecular_evolution/{gene}.cds.faa"
    shell:
        """
        python src/retrieve_cds.py {input}
        """

rule protein_alignment:
    input:
        "molecular_evolution/{gene}.cds.faa"
    output:
        "molecular_evolution/{gene}.aln.faa"
    shell:
        """
        mafft --localpair --maxiterate 1000 {input} >> {output}
        """

rule codon_alignment:
    input:
        "molecular_evolution/{gene}.aln.faa",
        "molecular_evolution/{gene}.cds.fna"
    output:
        "molecular_evolution/{gene}.aln.fna"
    shell:
        """
        pal2nal.pl {input} -output fasta >> {output}
        """

rule phylogeny_inference:
    input:
        "molecular_evolution/{gene}.aln.faa"
    output:
        "molecular_evolution/{gene}.treefile"
    params:
        "molecular_evolution/{gene}"
    shell:
        """
        iqtree -s {input} --prefix {params} -bb 1000 -T AUTO
        """

rule fel_selection_test:
    input:
        "{gene}.aln.fna",
        "{gene}.treefile"
    output:
        "{gene}.fel.json"
    shell:
        """
        hyphy fel --alignment {input[0]} --tree {input[1]} --output {output}
        """

rule meme_selection_test:
    input:
        "{gene}.aln.fna",
        "{gene}.treefile"
    output:
        "{gene}.meme.json"
    shell:
        """
        hyphy meme --alignment {input[0]} --tree {input[1]} --output {output}
        """

rule absrel_selection_test:
    input:
        "{gene}.aln.fna",
        "{gene}.treefile"
    output:
        "{gene}.absrel.json"
    shell:
        """
        hyphy absrel --alignment {input[0]} --tree {input[1]} --output {output}
        """
