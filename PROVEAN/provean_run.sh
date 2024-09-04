#!/bin/bash 

#bsub -q gpu -W 72:00 -n 20 -R rusage[mem=1000] -R span[hosts=1] -ojob41.out -ejob41.err ./job41.sh
#@!module load R/3.6.1	


#@!grep 'NONSYNONYMOUS' /project/uma_ana_caicedo/Hamid/20191118_deleterious_alleles/43_deleterious_mutations_shared_by_SIFT_and_PROVEAN_splittinSLL_wImp/EA00585_x_het_wANC_DER_SIFTannotations.xls > Nonsynonymous_SNPs_SIFT.txt

#@!cut -f5 Nonsynonymous_SNPs_SIFT.txt|uniq > protein_list.txt

#@!wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG2.4_release/ITAG2.4_proteins.fasta

#@!sed -e 's/\(^>.*$\)/#\1#/' ITAG2.4_proteins.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > ITAG2.4_proteins_linear.fasta

#@!while read protein_ID; do grep -A1 $protein_ID ITAG2.4_proteins_linear.fasta> $protein_ID.fasta; done < protein_list.txt 

#@!R CMD BATCH making_ancestral_protein_seqs.r

#compare resulting fasta files

#@!while read protein_ID; do comm -23 "$protein_ID".fasta  "$protein_ID"_ancestral.fasta >> fasta_comparison.txt;done < protein_list.txt

#@!mkdir fasta_files

#@!mv *.fasta ./fasta_files/

#@!R CMD BATCH making_variant_files.r



#@!module load provean/1.1.5
#@!module load cdhit/4.8.1 ncbi-blast/2.4.0+


#@!protein_list=(`ls ./fasta_files/*var|sed s@./fasta_files/@@|sed s/'.var'//`)
#@!mkdir provean_output
#@!for i in {1..21}; do
	#from original vcf, make smaller VCFs
#@!	a=`expr $i - 1` 
#@!	b=`expr $a \* 1372`
#@!	c=`expr $i \* 1372`
#@!	d=`expr $c - 1` 
#@!	for j in `seq $b $d`; do
	
#@!		~/software/provean/bin/provean.sh --num_threads 1 -q ./fasta_files/${protein_list[j]}_ancestral.fasta -v ./fasta_files/${protein_list[j]}.var > ./provean_output/${protein_list[j]}.provean_output

			
#@!	done &
#@!done

#bsub -q short -W 4:00 -n 1 -R rusage[mem=5000] -R span[hosts=1] -ojob41.out -ejob41.err ./job41.sh

#to organize the output
output_list=(`ls ./provean_output`)
num_output_files=`ls ./provean_output/*|wc -l|cut -f1 -d ' '`

for i in `seq 1 $num_output_files`; do
	a=`expr $i - 1`
	
	cat ./provean_output/${output_list[$a]}|grep ','|tr ',' '\t' > temp1
	echo ${output_list[$a]}|sed s/'.provean_output'// > temp2
	
	num_variants=`wc -l temp1|cut -f1 -d ' '`
	while read line; do for i in `seq 1 $num_variants`; do echo "$line"; done; done < temp2 > temp3
	
	paste temp1 temp3 > temp4
	
	cat temp4 >> all.provean_output
	
done	
	