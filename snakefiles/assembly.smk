rule all:
    input:
        "assembly/reads1.qc.fq",
        "assembly/reads2.qc.fq",
        "assembly/reads.qc.html",
        "assembly/reads.qc.json",
        "assembly/megahit/genome.contigs.fa",
        "assembly/genome_draft.fa",
        "assembly/genome.fa",
        "assembly/predicted_orfs.fa",
        "assembly/annotation.gtf"

rule quality_filter:
    input:
        "data/reads1.fq.gz",
        "data/reads2.fq.gz"
    output:
        "assembly/reads1.qc.fq",
        "assembly/reads2.qc.fq",
        html = "assembly/reads.qc.html",
        json = "assembly/reads.qc.json"
    shell:
        """
        fastp\
        --in1 {input[0]} --in2 {input[1]}\
        --out1 {output[0]} --out2 {output[1]}\
        --html {output.html} --json {output.json}\
        --reads_to_process 2000000
        """

rule genome_assembly:
    input:
        "assembly/reads1.qc.fq",
        "assembly/reads2.qc.fq"
    output:
        "assembly/megahit/genome.contigs.fa"
    params:
        name="genome"
    shell:
        """
        megahit\
        -1 {input[0]} -2 {input[1]}\
        --out-dir megahit --out-prefix '{params.name}'\
        && mv megahit assembly/
        """

rule circularize_sequence:
    input:
        "assembly/megahit/genome.contigs.fa"
    output:
        "assembly/genome_draft.fa"
    shell:
        """
        python src/circularize_sequence.py {input} |\
        fold -w 60 > {output}
        """

rule bwa_alignment:
    input:
        "assembly/genome_draft.fa",
        "assembly/reads1.qc.fq",
        "assembly/reads2.qc.fq"
    output:
        "assembly/reads_mapped.bam"
    shell:
        """
        bwa index {input[0]} &&\
        bwa mem {input} | samtools view -Sb - > {output}
        """

rule samtool_sort:
    input:
        "assembly/reads_mapped.bam"
    output:
        "assembly/reads_sorted.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule lofreq_snv:
    input:
        genome="assembly/genome_draft.fa",
        bam="assembly/reads_sorted.bam"
    output:
        "assembly/snv.vcf"
    shell:
        """
        lofreq call\
        -f {input.genome} -o {output} {input.bam}
        """

rule correct_variants:
    input:
        genome="assembly/genome_draft.fa",
        vcf="assembly/snv.vcf"
    output:
        "assembly/genome.fa"
    shell:
        """
        python src/correct_variants.py {input.genome} {input.vcf} |\
        fold -w 60 > {output}
        """

rule orf_detection:
    input:
        "assembly/genome.fa"
    output:
        "assembly/predicted_orfs.fa"
    shell:
        """
        ./src/ORFfinder -in {input} -out {output} -c t -s 0 -ml 150
        """

rule orf_homology_search:
    input:
        "assembly/predicted_orfs.fa"
    output:
        "assembly/orf_blastp.xml"
    params:
        "data/reference_proteome.faa",
        "assembly/blast_db"
    shell:
        """
        makeblastdb -in {params[0]} -out {params[1]} -dbtype prot &&\
        blastp -db {params[1]} -query {input}\
               -evalue 0.0001 -max_target_seqs 1 -outfmt 5\
               -out {output}
        """

rule orf_annotation:
    input:
        "assembly/orf_blastp.xml",
        "data/protein_names.csv"
    output:
        "assembly/annotation.gtf"
    shell:
        """
        python src/annotate_blast.py {input} > {output}
        """
