SAMPLES = ["sfmnpv"]

rule all:
    input:
        "annotation/sfmnpv.gb",
        "isolates_diversity/tree_plot.svg",
        "isolates_diversity/rooted_phylogeny/tree_plot.svg",
        "isolates_diversity/sequence_similarity_plot.svg"

rule quality_filter:
    input:
        "data/{sample}_1.fq.gz",
        "data/{sample}_2.fq.gz"
    output:
        "genome_assembly/{sample}_sample1.fq",
        "genome_assembly/{sample}_sample2.fq",
        html = "genome_assembly/{sample}_sample.html",
        json = "genome_assembly/{sample}_sample.json"
    shell:
        """
        fastp\
        --in1 {input[0]} --in2 {input[1]}\
        --out1 {output[0]} --out2 {output[1]}\
        --html {output.html} --json {output.json}\
        --reads_to_process 1000000
        """

rule genome_assembly:
    input:
        "genome_assembly/{sample}_sample1.fq",
        "genome_assembly/{sample}_sample2.fq"
    output:
        "genome_assembly/megahit/{sample}.contigs.fa"
    params:
        name="{sample}"
    shell:
        """
        megahit\
        -1 {input[0]} -2 {input[1]}\
        --out-dir megahit --out-prefix '{params.name}'\
        && mv megahit genome_assembly/
        """

rule rearrange_sequence:
    input:
        "genome_assembly/megahit/sfmnpv.contigs.fa"
    output:
        "genome_assembly/sfmnpv_draft.fa"
    shell:
        """
        python src/rearrange_sequence.py {input} |\
        fold -w 60 > {output}
        """

rule bwa_alignment:
    input:
        "genome_assembly/sfmnpv_draft.fa",
        "genome_assembly/sfmnpv_sample1.fq",
        "genome_assembly/sfmnpv_sample2.fq"
    output:
        "genome_assembly/reads.bam"
    shell:
        """
        bwa index {input[0]} &&\
        bwa mem {input} | samtools view -Sb - > {output}
        """

rule samtool_sort:
    input:
        "genome_assembly/reads.bam"
    output:
        "genome_assembly/reads_sort.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule lofreq_snv:
    input:
        genome="genome_assembly/sfmnpv_draft.fa",
        bam="genome_assembly/reads_sort.bam"
    output:
        "genome_assembly/snv.vcf"
    shell:
        """
        lofreq call -f {input.genome} -o {output} {input.bam}
        """

rule correct_alleles:
    input:
        genome="genome_assembly/sfmnpv_draft.fa",
        vcf="genome_assembly/snv.vcf"
    output:
        "genome_assembly/sfmnpv.fa"
    shell:
        """
        python src/alternative_alleles.py {input.genome} {input.vcf} |\
        fold -w 60 > {output}
        """

rule orf_detection:
    input:
        "genome_assembly/sfmnpv.fa"
    output:
        "annotation/sfmnpv_orf.fa"
    shell:
        """
        ./src/ORFfinder -in {input} -out {output} -c t -s 0 -ml 150
        """

rule annotation_blast_search:
    input:
        "annotation/sfmnpv_orf.fa"
    output:
        "annotation/orf_blastp.xml"
    params:
        "data/spodoptera_sp_db.faa",
        "annotation/blast_db"
    shell:
        """
        makeblastdb -in {params[0]} -out {params[1]} -dbtype prot &&\
        blastp -db {params[1]} -query {input}\
               -evalue 0.0001 -max_target_seqs 1 -outfmt 5\
               -out {output}
        """

rule orf_annotation:
    input:
        "annotation/orf_blastp.xml"
    output:
        "annotation/sfmnpv_annotation.gtf"
    shell:
        """
        python src/annotate_blast.py {input} > {output}
        """

rule build_feature_table:
    input:
        "annotation/sfmnpv_annotation.gtf"
    output:
        "annotation/sfmnpv.tbl"
    shell:
        """
        python src/gtf2tbl.py {input} > {output}
        """

rule create_asn:
    input:
        "data/template.sbt",
        "annotation/sfmnpv.tbl"
    output:
        "annotation/sfmnpv.sqn"
    shell:
        """
        cp genome_assembly/sfmnpv.fa annotation/sfmnpv.fsa &&\
        src/tbl2asn -t {input[0]} -i annotation/sfmnpv.fsa -f {input[1]} 
        """

rule create_genbank:
    input:
        "annotation/sfmnpv.sqn"
    output:
        "annotation/sfmnpv.gb"
    shell:
        """
        src/asn2gb -i {input} > {output} 
        """

rule join_sfmnpv_genomes:
    input:
        "data/sfmnpv_genomes_ncbi.fna",
        "genome_assembly/sfmnpv.fa"
    output:
        "isolates_diversity/sfmnpv_genomes.fna"
    shell:
        """
        cat {input} > {output}
        """
rule align_sfmnpv_genomes:
    input:
        "isolates_diversity/sfmnpv_genomes.fna"
    output:
        "isolates_diversity/sfmnpv_genomes_aln.fna"
    shell:
        """
        mafft {input} >> {output}
        """

rule build_sfmnpv_phylogeny:
    input:
        "isolates_diversity/sfmnpv_genomes_aln.fna"
    output:
        "isolates_diversity/sfmnpv_genomes.treefile"
    params:
        "isolates_diversity/sfmnpv_genomes"
    shell:
        """
        iqtree -s {input} -pre {params} -bb 1000 -redo
        """

rule plot_sfmnpv_phylogeny:
    input:
        "isolates_diversity/sfmnpv_genomes.treefile"
    output:
        "isolates_diversity/tree_plot.svg"
    params:
        "KF891883.1"
    shell:
        """
        python src/draw_tree.py {input} {params} &&\
        mv tree_plot.svg {output}
        """

rule plot_sequence_similarity:
    input:
        "isolates_diversity/sfmnpv_genomes_aln.fna"
    output:
        "isolates_diversity/sequence_similarity_plot.svg"
    shell:
        """
        python src/sequence_similarity_profile.py {input} 200 2000 &&\
        mv sequence_similarity_plot.svg {output}
        """

rule join_sfmnpv_root_genomes:
    input:
        "data/semnpv_root.fna",
        "isolates_diversity/sfmnpv_genomes.fna"
    output:
        "isolates_diversity/rooted_phylogeny/sequences.fna"
    shell:
        """
        cat {input} > {output}
        """
rule align_root_genomes:
    input:
        "isolates_diversity/rooted_phylogeny/sequences.fna"
    output:
        "isolates_diversity/rooted_phylogeny/sequences_aln.fna"
    shell:
        """
        mafft {input} >> {output}
        """

rule build_rooted_phylogeny:
    input:
        "isolates_diversity/rooted_phylogeny/sequences_aln.fna"
    output:
        "isolates_diversity/rooted_phylogeny/rooted_phylogeny.treefile"
    params:
        "isolates_diversity/rooted_phylogeny/rooted_phylogeny"
    shell:
        """
        iqtree -s {input} -pre {params} -bb 1000 -redo
        """

rule plot_rooted_phylogeny:
    input:
        "isolates_diversity/rooted_phylogeny/rooted_phylogeny.treefile"
    output:
        "isolates_diversity/rooted_phylogeny/tree_plot.svg"
    params:
        "NC_002169.1"
    shell:
        """
        python src/draw_tree.py {input} {params} && \
        mv tree_plot.svg {output}
        """

rule plot_isolates_snv_distribution:
    input:
        "isolates_diversity/sfmnpv_genomes_snv.vcf"
    output:
        "isolates_diversity/isolates_snv_distribution.svg"
    shell:
        """
        python src/isolates_variants_distribution.py {input} && \
        mv isolates_snv_distribution.svg {output}
        """
