protein_list=read.table("protein_list.txt", header=F)
protein_list=as.character(protein_list[,1])


Nonsynonymous_SNPs=read.table("Nonsynonymous_SNPs_SIFT.txt", header=F, sep="\t", check.names=F)

for(i in 1:length(protein_list)){

	temp=Nonsynonymous_SNPs[Nonsynonymous_SNPs[,5]==protein_list[i],]
	
	temp2=temp[,c(12,10,11)]

	write.table(temp2, paste0("./fasta_files/",protein_list[i], ".var"), quote=F, row.names=F, col.names=F, sep=",")
}
