protein_list=read.table("protein_list.txt", header=F)
protein_list=as.character(protein_list[,1])


Nonsynonymous_SNPs=read.table("Nonsynonymous_SNPs_SIFT.txt", header=F, sep="\t", check.names=F)

for(i in 1:length(protein_list)){

	protein=read.table(paste0(protein_list[i],".fasta"), header=T, sep="\t", check.names=F)

	temp=as.character(protein[,1])

	temp2=Nonsynonymous_SNPs[Nonsynonymous_SNPs[,5]==protein_list[i],]

	for (j in 1:nrow(temp2)){

		pos=as.numeric(temp2[j,12])
		
		#to double check this is working

		if (substr(temp,pos,pos)==as.character(temp2[j,10])) {
			print(paste0(pos,"_",protein_list[i]))	
		}
		
		substr(temp,pos,pos)=as.character(temp2[j,10])

	}
	
	protein[,1]=temp
	write.table(protein, paste0(protein_list[i], "_ancestral.fasta"), quote=F, row.names=F, col.names=T)
}
