dadi_file=read.table("input.vcf.data", header=T, check.names=F)
pops=c("SP", "SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX_CA_NSA", "SLC_MEX", "SLL_Americas","SLL_modern")


ref_alleles_counts=dadi_file[,4:(length(pops)+3)]
#to sort columns by populations
ref_alleles_counts=ref_alleles_counts[,pops]


#counts of alternate alleles
alt_alleles_counts=dadi_file[,(length(pops)+5):((length(pops)+5)+length(pops)-1)]
#to sort columns by populations
alt_alleles_counts=alt_alleles_counts[,pops]


#to find private reference alleles
summary_pops_num_derived=data.frame(pops=pops, num_sites_w_derived_alleles=0, DAF=0)
for (i in 1:length(pops)) {
	
	temp=cbind(ref_alleles_counts[i],alt_alleles_counts[i])
	
	temp$Ref=dadi_file$Ref
	temp$OUT=dadi_file$OUT

	
	#count of derived alleles in each site
	temp$x=0
	temp[temp$Ref==temp$OUT,]$x=temp[temp$Ref==temp$OUT,2]
	temp[temp$Ref!=temp$OUT,]$x=temp[temp$Ref!=temp$OUT,1]

	temp$chromosome=dadi_file$Gene
	temp$position=dadi_file$Postion
	

		
	num_sites_w_derived_alleles_all_chr=data.frame()
	for(j in c("01","02","03","04","05","06","07","08","09","10","11","12")){
		temp2=temp[temp$chromosome==paste0("SL2.50ch",j),]
		
		last_window=ceiling(temp2[nrow(temp2),"position"]/100000)
		
		num_sites_w_derived_alleles=data.frame(chromosome=j, windows=1:last_window, num_sites_w_derived_alleles=0)
		for(k in 1:last_window){
			temp3=temp2[temp2$position > (k-1)*100000 & temp2$position <= (k)*100000 ,]
						
			if(nrow(temp3) > 0) {
				#finding the number of sites with derived alleles
				num_sites_w_derived_alleles$num_sites_w_derived_alleles[k]=nrow(temp3[temp3$x > 0,])
				
				if(is.na(nrow(temp3[temp3$x > 0,]))) {
				num_sites_w_derived_alleles$num_sites_w_derived_alleles[k]=0
				}
			}
			
		}
			
		num_sites_w_derived_alleles_all_chr=rbind(num_sites_w_derived_alleles_all_chr, num_sites_w_derived_alleles)
		
	}
		
	write.table(num_sites_w_derived_alleles_all_chr, paste0(pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), quote=F, row.names=F, col.names=T)	


	#to check progress
	print(pops[i])
	
	
}	
