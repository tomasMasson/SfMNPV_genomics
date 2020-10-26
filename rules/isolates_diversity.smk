rule concatenate_genomes:
    input:
        "data/isolates_genomes.fna",
        "assembly/genome.fa"
    output:
        "isolates_diversity/genomes.fna"
    shell:
        """
        cat {input} > {output}
        """

rule align_genomes:
    input:
        "isolates_diversity/genomes.fna"
    output:
        "isolates_diversity/genomes.aln.fna"
    shell:
        """
        mafft {input} >> {output}
        """

rule build_phylogeny:
    input:
        "isolates_diversity/genomes.aln.fna"
    output:
        "isolates_diversity/phylogeny.treefile"
    params:
        "isolates_diversity/phylogeny"
    shell:
        """
        iqtree -s {input} -nt 4 -pre {params} -bb 1000 -redo
        """

rule add_root_alignment:
    input:
        root="data/semnpv.fna",
        msa="isolates_diversity/genomes.aln.fna"
    output:
        "isolates_diversity/rooted_phylogeny/genomes.aln.fna"
    shell:
        """
        mafft --add {input.root} {input.msa} >> {output}
        """

rule build_rooted_phylogeny:
    input:
        "isolates_diversity/rooted_phylogeny/genomes.aln.fna"
    output:
        "isolates_diversity/rooted_phylogeny/phylogeny.treefile"
    params:
        "isolates_diversity/rooted_phylogeny/phylogeny"
    shell:
        """
        iqtree -s {input} -o NC_002169.1 -nt 4 -pre {params} -bb 1000 -redo
        """

rule extract_cds:
    input:
        "assembly/genome.fa",
        "assembly/annotation.gtf"
    output:
        "isolates_diversity/proteome.faa",
    shell:
        """
        python src/extract_protein_cds.py {input} > {output}
        """

rule cluster_orthologous_proteins:
    input:
        "isolates_diversity/proteome.faa",
        "data/isolates_proteomes.faa"
    output:
        "isolates_diversity/blast.xml"
    params:
        "isolates_diversity/proteomes.db"
    shell:
        """
        cat {input} > {params} && \
        makeblastdb -in {params} -dbtype prot && \
        blastp -query {input[0]} -db {params} -evalue 0.0001 -outfmt 5 -max_target_seqs 6 > {output}
        """

rule count_missense_sites:
    input:
        "isolates_diversity/blast.xml"
    output:
        "isolates_diversity/gene_missense_sites.csv"
    params:
        "isolates_diversity/proteomes.db"
    shell:
        """
        mkdir -p isolates_diversity/orthogroups && \
        python src/get_orthogroups.py {input} {params} && \
        mv cds* isolates_diversity/orthogroups &&\
        for file in isolates_diversity/orthogroups/*; do mafft --localpair --maxiterate 1000 $file > $file.aln;done &&\
        for file in isolates_diversity/orthogroups/*.aln; do src/count_missense_sites.py $file >> {output};done
        """

rule extract_snp:
    input:
        "isolates_diversity/genomes.aln.fna"
    output:
        "isolates_diversity/snp.vcf"
    shell:
        """
        snp-sites -cv {input} > {output}
        """
