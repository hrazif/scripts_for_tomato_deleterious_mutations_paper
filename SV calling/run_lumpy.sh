#!/bin/bash 

#nohup ./run_lumpy.sh>>run_lumpy.output&

#run on Cornell bioHPC 




#reference genome, Heinz 2.5
#wget ftp://ftp.solgenomics.net/tomato_genome/Heinz1706/assembly/build_2.50/S_lycopersicum_chromosomes.2.50.fa.gz

#to index the genome
#bwa index S_lycopersicum_chromosomes.2.50.fa.gz

#@!cp /local/workdir/hr355/20211228_SV_calling_short_read_data/20211029_organized_genome_linear.* ./


#@!mkdir lumpy

cd lumpy 

bams=/local/workdir/hr355/20191118_deleterious_alleles/20220425_SV_calling_LA2093_shortread


#@!for i in `ls $bams/*.bam|sed s@/home/FrankLab/hr355/Analyses/20211115_variant_calling_w_new_assembly/@@g|sed s@.bam@@g`; do
for i in SRR12039813; do


# Extract the discordant paired-end alignments.
samtools view --threads 20 -b -F 1294 $bams/"$i".bam > "$i".discordants.unsorted.bam

# Extract the split-read alignments
samtools view --threads 20 -h $bams/"$i".bam \
    | /programs/lumpy-sv-0.3.0/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view --threads 20 -Sb - \
    > "$i".splitters.unsorted.bam

# Sort both alignments
samtools sort --threads 20 "$i".discordants.unsorted.bam -o "$i".discordants.bam
samtools sort "$i".splitters.unsorted.bam -o "$i".splitters.bam


samtools view $bams/"$i".bam \
    | tail -n+100000 \
    | /programs/lumpy-sv-0.3.0/scripts/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 10000 \
    -o "$i".lib1.histo

#get mean and sd from above, then
/home/hr355/Software/lumpy-sv/bin/lumpy \
   -mw 4 \
   -tt 0 \
    -pe id:"$i",bam_file:"$i".discordants.bam,histo_file:"$i".lib1.histo,mean:457,stdev:20,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:"$i",bam_file:"$i".splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    > "$i".vcf
	
	

done	

bams=/local/workdir/hr355/20191118_deleterious_alleles/20220425_SV_calling_LA2093_shortread

#for SV typer
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
export PYTHONPATH=/programs/svtyper-0.7.1/lib/python2.7/site-packages:/programs/svtyper-0.7.1/lib64/python2.7/site-packages
export PATH=/programs/svtyper-0.7.1/bin:$PATH

for i in SRR12039813; do
	
	svtyper -B $bams/"$i".bam  -l "$i".json
	svtyper -B $bams/"$i".bam -i "$i".vcf -o "$i".gt.vcf
	
done	