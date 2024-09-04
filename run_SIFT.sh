#!/bin/bash


#nohup ./job24.sh >> job24.out &


input_vcf=/home/caicedo/Hamid/20200406/2_Imputation_including_outgroup/LinkImputeOutput_all.vcf.gz

#to extract Pennelli
~/usr/bin/vcftools --gzvcf $input_vcf --remove-indels --indv "EA00585" --max-missing 1.0 --recode --stdout|gzip -c > EA00585.vcf.gz

#to remove header and heterozgous sites in SPENN
zcat EA00585.vcf.gz|grep -v '##'|sed s/'#CHROM'/'CHROM'/g| grep -v '0/1' > EA00585_x_het.txt

#run to_prepare_input_VCF_for_SITF.r, then

R CMD BATCH to_prepare_input_VCF_for_SITF.r

cat EA00585_x_het_wANC_DER.txt|awk '{print $1"\t"$2"\t"$3"\t"$11"\t"$12"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}'|sed 1d > EA00585_x_het_wANC_DER_no_colnames.txt

cat header.txt EA00585_x_het_wANC_DER_no_colnames.txt|sed s'/SL2.50ch0//'g |sed s'/SL2.50ch//'g > EA00585_x_het_wANC_DER.vcf

#java -jar <Path to SIFT4G_Annotator> -c -i <Path to input vcf file> -d <Path to SIFT4G database directory> -r <Path to your results folder> -t
annotator=/home/caicedo/software/sift4g/SIFT4G_Annotator_v2.4.jar
database=/home/caicedo/Hamid/20191104/2_SIFT/tomato_database/scripts_to_build_SIFT_db/tomato_2.50/SL2.50
java -jar $annotator -c -i EA00585_x_het_wANC_DER.vcf -d $database -r ./ -t