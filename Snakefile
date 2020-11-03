rule all:
    input:
      "assembly/annotation.gtf",
      "variants_calling/snv_ann.vcf",
      "variants_calling/delly_sv.bcf",
      "isolates_diversity/snp.vcf",
      "isolates_diversity/gene_missense_sites.csv",
      "isolates_diversity/phylogeny.treefile",
      "isolates_diversity/rooted_phylogeny/phylogeny.treefile",
      "molecular_evolution/sf29.fel.json",
      "molecular_evolution/sf29.meme.json",
      "molecular_evolution/sf29.absrel.json"

include: "rules/assembly.smk"
include: "rules/variant_calling.smk"
include: "rules/isolates_diversity.smk"
include: "rules/molecular_evolution.smk"
