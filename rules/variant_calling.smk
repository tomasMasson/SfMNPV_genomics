rule reads_quality_control:
    input:
        "data/reads1.fq.gz",
        "data/reads2.fq.gz"
    output:
        "variants_calling/reads1.qc.fq",
        "variants_calling/reads2.qc.fq",
        html = "variants_calling/reads.html",
        json = "variants_calling/reads.json"
    shell:
        """
        fastp\
        --in1 {input[0]} --in2 {input[1]}\
        --out1 {output[0]} --out2 {output[1]}\
        --html {output.html} --json {output.json}
        """

rule bwa_index:
    input:
        "assembly/genome.fa"
    output:
        "variants_calling/genome.fa",
        "variants_calling/genome.fa.amb",
        "variants_calling/genome.fa.ann",
        "variants_calling/genome.fa.bwt",
        "variants_calling/genome.fa.pac",
        "variants_calling/genome.fa.sa"
    params:
        "variants_calling/genome.fa"
    shell:
        """
        cp {input} {params} &&\
        bwa index {params}\
        """

rule bwa_alignment:
    input:
        "variants_calling/genome.fa",
        "variants_calling/reads1.qc.fq",
        "variants_calling/reads2.qc.fq"
    output:
        "variants_calling/reads_mapped.bam"
    params:
        rg=r"@RG\tID:genome_assembly\tSM:1\tLB:lib1\tPL:Illumina\tPU:unit1"
    threads: 4
    shell:
        """
        bwa mem -R '{params.rg}' -t {threads} {input} |\
        samtools view -Sb > {output}
        """

rule samtool_sort:
    input:
        "variants_calling/reads_mapped.bam"
    output:
        "variants_calling/reads_sorted.bam",
        "variants_calling/reads_sorted.bam.bai",
    threads: 4
    shell:
        """
        samtools sort -@ {threads} {input} > {output[0]} && \
        samtools index -@ {threads} {output[0]}
        """

rule lofreq_snv_calling:
    input:
        genome="variants_calling/genome.fa",
        bam="variants_calling/reads_sorted.bam"
    output:
        "variants_calling/snv.vcf"
    threads: 4
    shell:
        """
        samtools faidx {input.genome} &&\
        lofreq call-parallel --pp-threads {threads}\
        -f {input.genome} -o {output} {input.bam}
        """

rule snpEff_snv_annotation:
    input:
        "variants_calling/snv.vcf"
    output:
        "variants_calling/snv_ann.vcf"
    params:
        config="data/snpEff.config",
        db="genome_assembly",
        seq="variants_calling/snpEff_db/genome_assembly/sequences.fa",
        gtf="variants_calling/snpEff_db/genome_assembly/genes.gtf"
    shell:
        """
        mkdir -p variants_calling/snpEff_db/genome_assembly && \
        cp assembly/genome.fa {params.seq} && \
        src/create_snpeff_gtf.py assembly/annotation.gtf > {params.gtf} && \
        snpEff build -gtf22 -c {params.config} {params.db} && \
        snpEff -c {params.config} {params.db} {input} > {output} &&\
        rm snpEff_*
        """

rule delly_sv_calling:
    input:
        genome="variants_calling/genome.fa",
        bam="variants_calling/reads_sorted.bam"
    output:
        "variants_calling/delly_sv.bcf"
    shell:
        """
        delly call -g {input.genome} -o {output} {input.bam}
        """
