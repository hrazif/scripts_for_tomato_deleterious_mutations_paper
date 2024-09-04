#!/bin/bash 

#run on Cornell bioHPC (cbsumm14)

#nohup ./run_manta.sh>>run_manta.output&

#to configure manta for each accession

mkdir manta
cd manta

bams=/local/workdir/hr355/20191118_deleterious_alleles/20220425_SV_calling_LA2093_shortread

genome=/local/workdir/hr355/20191118_deleterious_alleles/20220425_SV_calling_LA2093_shortread/S_lycopersicum_chromosomes.4.00.fa.bgz

#for i in `ls $bams/*.bam|sed s@/home/FrankLab/hr355/Analyses/20211115_variant_calling_w_new_assembly/@@g|sed s@.bam@@g`; do
for i in SRR12039813;do 

	mkdir "$i"; cd "$i"

	~/Software/manta-1.6.0.centos6_x86_64/bin/configManta.py --bam=$bams/"$i".bam --referenceFasta=$genome

	cd MantaWorkflow/
	
	./runWorkflow.py
	
	cd ..
	cd ..
done


#@!cd manta

#@!mkdir manta_output_vcfs

#@!for i in `ls $bams/*.bam|sed s@/home/FrankLab/hr355/Analyses/20211115_variant_calling_w_new_assembly/@@g|sed s@.bam@@g`; do
#@!for i in SRR12039813; do

#@!mv "$i"/MantaWorkflow/results/variants/diploidSV.vcf.gz ./manta_output_vcfs/"$i".vcf.gz
#@!	rm -r "$i"
#@!done
