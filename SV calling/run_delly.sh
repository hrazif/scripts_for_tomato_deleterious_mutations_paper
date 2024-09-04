#!/bin/bash 

#run on Cornell bioHPC (cbsumm20)

#nohup ./run_delly.sh>>run_delly.output&





#reference genome
#@!cp /home/FrankLab/shared/genomes/eggplant/20211029_organized_genome_linear.fasta ./


#@!bgzip -c 20211029_organized_genome_linear.fasta > 20211029_organized_genome_linear.fasta.bgz

#to index the genome
#@!bwa index 20211029_organized_genome_linear.fasta.bgz


#@!mkdir delly; 

cd delly


#to run delly
export PATH=/programs/delly-0.8.7:$PATH

bam_folder=/local/workdir/hr355/20191118_deleterious_alleles/20220425_SV_calling_LA2093_shortread

#for i in `ls $bam_folder/*.bam|sed s@/home/FrankLab/hr355/Analyses/20211115_variant_calling_w_new_assembly/@@g|sed s@.bam@@g`; do
#for i in MM01491 MM01493 MM01506 MM01536 MM01543 MM01545 MM01547 MM01553 MM01565 MM01572 MM01582 MM01584 MM01588 MM01592 MM01597 MM01654 MM01658 MM01659 MM01676 MM01677 MM01678 MM01712 MM01790 MM01791 MM01792 MM01808 MM01826 MM01828 MM01831 MM01838 MM10417 MM10438 MM10439 MM12137 MM12437 PI105347 PI116953 PI140456 PI141970 PI143402 PI163268 PI164581 PI164710 PI166362 PI179997 PI180343 PI200854 PI200856 PI222267 PI279872 PI462370; do

genome=/local/workdir/hr355/20191118_deleterious_alleles/20220425_SV_calling_LA2093_shortread/S_lycopersicum_chromosomes.4.00.fa.bgz
for i in SRR12039813; do

delly call -g $genome $bam_folder/"$i".bam -o "$i".bcf

bcftools view "$i".bcf > "$i".vcf

done



