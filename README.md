# Genomic Diversity in a Population of Spodoptera frugiperda Nucleopolyhedrovirus

<img src="thumbnail.png" width="400" style="float:right"/>

This project explored the genomic diversity inside an argentinian (ARG-M) isolate of *Spodoptera frugiperda* multiple nucleopolyhedrovirus (SfMNPV).

## SfMNPV genomic diversity (assembly/ and variants_calling/)

A consensus genome sequence for SfMNPV ARGM-M was assembled (assembly/genome.fa) with [Megahit](https://github.com/voutcn/megahit) from a subset of Illumina reads. Single nucleotide and structural variants were discovered through Illumina deep sequencing together with [LoFreq](http://csb5.github.io/lofreq/), [Delly](https://github.com/dellytools/delly) and [Lumpy](https://github.com/arq5x/lumpy-sv) as variant callers (variants_calling/snv_report, variants_calling/delly_sv.bcf and variants_calling/lumpy-sv.vcf).

## Genetic diversity present in SfMNPV isolates (isolates_diversity/)

All SfMNPV sequenced isolates were aligned using [MAFFT](https://mafft.cbrc.jp/alignment/software/) with default settings (isolates_diversity/genomes.aln.fna) and single polymorphisms (isolates_diversity/snp.vcf) were extracted with [SNP-sites](https://github.com/sanger-pathogens/snp-sites). Additionally, a maximum likelihood phylogeny (isolates_diversity/phylogeny.treefile) was reconstructed with [IQ-TREE](http://www.iqtree.org/).

Coding sequences genetic diversity was compared (viz/panels/cds_diversity.svg) using the missense mutations within and between SfMNPV isolates (viz/circos/gene_missense_mutations.csv).

## Molecular evolution of *sf29* (molecular_evolution/)

BLASTp search results (molecular_evolution/sf29_blast.xml), phylogenetic reconstruction (molecular_evolution/sf29.aln.faa.treefile) and HyPhy evolutionary rate inferences (FEL, MEME and aBSREL) are provided for the putative collagenase *sf29*.
