#!/bin/bash

bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" variants_calling/genome.fa variants_calling/reads1.qc.fq variants_calling/reads2.qc.fq | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -Sb - > variants_calling/lumpy-sv.bam

samtools view -b -F 1294 variants_calling/lumpy-sv.bam > variants_calling/discordant_lumpy-sv_unsorted.bam 

samtools view -h variants_calling/lumpy-sv.bam | ~/miniconda3/envs/lumpy/bin/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > variants_calling/splitters_lumpy-sv_unsorted.bam

samtools sort variants_calling/discordant_lumpy-sv_unsorted.bam -o variants_calling/discordant_lumpy-sv.bam

samtools sort variants_calling/splitters_lumpy-sv_unsorted.bam -o variants_calling/splitters_lumpy-sv.bam

samtools view variants_calling/lumpy-sv.bam | tail -n 100000 | ~/miniconda3/envs/lumpy/share/lumpy-sv-0.3.0-3/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o variants_calling/lib1.histo

lumpy -mw 4 -tt 0 -pe id:sfmnpv,bam_file:variants_calling/discordant_lumpy-sv.bam,histo_file:variants_calling/lib1.histo,mean:439,stdev:118,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:sfmnpv,bam_file:variants_calling/splitters_lumpy-sv.bam,back_distance:10,weight:1,min_mapping_threshold:20 > variants_calling/lumpy-sv.vcf
