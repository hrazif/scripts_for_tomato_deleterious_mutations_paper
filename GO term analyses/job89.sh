#!/bin/bash 

# bsub -q gpu -W 72:00 -n 1 -R rusage[mem=4000] -R span[hosts=1] -o job89.out -e job89.err ./job89.sh

deleterious_vcf=/project/uma_ana_caicedo/Hamid/20191118_deleterious_alleles/43_deleterious_mutations_shared_by_SIFT_and_PROVEAN_splittinSLL_wImp/deleterious_SNPs.vcf.gz

#to convert to dadi format
zcat $deleterious_vcf > deleterious.vcf

perl ~/software/convert_vcf_to_dadi_input.pl deleterious.vcf dadi_pop.txt


module load R/3.5.0
R CMD BATCH finding_derived_deleterious.r
	

#annotate
R CMD BATCH GO_term_enrichment_all_topGO.r



