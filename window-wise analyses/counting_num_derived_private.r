individuals=c("BGV006208","BGV007339","BGV007151","BGV007366","BGV007198","TS-418","TS-440","TS-417","TS-79","BGV006336","BGV006457","TR00028","BGV006363","BGV006478","TS-435","TS-415","TS-434","EA00676","TS-421","TS-439","TS-422","TS-425","TS-77","TS-429","BGV006327","BGV006454","TS-419","BGV006347","BGV006353","BGV015382","TS-420","BGV015380","PAS014479","TS-438","BGV006370","TS-416","TS-156","TS-411","BGV007109","BGV007111","BGV006775","BGV007194","BGV007149","BGV007169","BGV007152","BGV007161","BGV007158","BGV007181","BGV006230","PI129033","BGV006828","BGV006881","BGV012625","TS-71","BGV007992","TS-148","TR00027","TS-158","BGV006865","BGV007023","BGV006867","BGV006225","BGV006232","BGV012639","TS-53","BGV006806","BGV006859","BGV006896","BGV006904","BGV006235","BGV006229","BGV006910","BGV006906","BGV006907","BGV006753","BGV006231","BGV006899","BGV006901","BGV007015","BGV007017","BGV006175","PI487625","BGV006234","BGV005912","BGV006927","BGV006931","BGV006934","BGV008058","TS-299","PI129026","BGV006777","BGV006852","BGV006792","BGV006768","BGV006779","BGV008095","BGV008096","BGV013161","TS-56","BGV007990","BGV007989","TS-273","BGV006825","BGV013945","BGV008225","TS-57","TS-231","TS-258","BGV007981","BGV008189","BGV014508","TS-256","BGV008100","TS-129","BGV015730","PI406890","BGV008037","BGV008041","TS-247","BGV008077","BGV004584","BGV012640","BGV008036","BGV008061","TS-304","BGV008065","TS-184","BGV014522","BGV015726","LA2309","BGV008098","BGV014516","BGV014514","BGV014515","BGV014518","BGV015734","BGV014519","BGV015727","BGV008354","BGV012627","PI129088","BGV008348","BGV008345","BGV008218","LA1712","BGV008108","LA2697","BGV013175","BGV007927","BGV007931","TS-436","BGV007933","BGV007934","BGV007935","TS-302","TS-131","BGV005895","TS-280","BGV008106","CATIE-11106-1","BGV008219","BGV008221","BGV007900","BGV007899","BGV007894","TS-229","BGV012614","BGV008067","BGV008223","BGV007920","TS-154","TS-165","BGV008051","BGV007902","BGV013134","BGV007901","BGV008070","BGV007909","BGV007908","BGV007918","BGV007911","BGV007910","BGV007921","BGV008224","BGV007871","BGV007878","BGV007895","BGV007875","BGV007857","BGV007865","BGV007876","BGV007936","BGV007872","TS-249","BGV007860","BGV007864","BGV007867","BGV007870","EA01640","EA00371","TS-152","BGV007862","BGV007863","EA04243","TS-75","EA05701","TR00019","Tegucigalpa","TS-261","TS-10","TS-141","EA00940","EA03701","EA01835","TS-44","EA01854","EA04861","TS-282","TS-163","TR00022","TS-73","TS-193","TR00018","TS-251","TR00021","TS-173","TS-86","AlisaCraig","TS-142","EA02617","TS-43","Moneymaker","TS-95","EA00892","EA04828","EA00465","TR00023","EA01155","EA03221","TS-132","TS-192","EA00990","TS-191","EA01037","TR00020","EA00448","EA05581","TS-197","EA00157","EA01019","EA01049", "EA01088")

all_genotypes=read.table("private_del.GT.FORMAT", header=T, stringsAsFactors = F, check.names=F)


output=data.frame(indivdual=individuals, num_derived=0)
for (i in 1:length(individuals)) {
	
	genotypes=all_genotypes[,c(individuals[i], "EA00585")]
	
	#to remove sites with missing or heterozygous ancestral state
	genotypes=genotypes[genotypes[,1]!="." | genotypes[,2]!="." | genotypes[,2]!="0/1", ]
	
		
	genotypes$derived=0
	
	for(j in 1:nrow(genotypes)){
	
		if(genotypes[j,1]!=genotypes[j,2]) {
		
		#to adjust for homozygosity
		
			if(genotypes[j,1]=="0/1"){
				genotypes$derived[j]=0.5
			}else{
				genotypes$derived[j]=1
			}	
		}
	}
	
	output$num_derived[i]=sum(genotypes$derived)
	}


populations=read.table("populations_splittinSLL.txt", header=T, stringsAsFactors=F)

output$population="x"
for (k in 1:nrow(output)){

	output$population[k]=populations[populations$accession==output$indivdual[k], "population"]


}

write.table(output, "num_deleterious_inds.txt", quote=F, sep="\t", col.names=T, row.names=F)


#for plotting
pop_names=c("SP", "SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX_CA_NSA", "SLC_MEX", "SLL_Americas", "SLL_modern")

pop_size=vector()
for (j in pop_names) {
  temp=dim(populations[populations$population==j,])[1]
  pop_size=c(pop_size, temp)
}

box_width=sqrt(pop_size)
space=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1)

pop_names_simple=names=c("SP","SLC ECU","SLC PER","SLC San Martin","SLC MEX-CA-NSA","SLC MEX","SLL Americas", "SLL modern")
#for combined pdfs
pdf("num_deleterious_inds.pdf", height = 4, width = 4)
      
  #to order by population
  output$population<-ordered(output$population, levels=pop_names)
  
  par(mfrow=c(1,1), mai=c(1.7, 0.9, 0.3, 0.05)*0.75)
  at.x=seq(0.1,0.7*length(pop_names),0.7)
  plot (num_derived~population, data=output, main="", col=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3"), varwidth=F, xlab="", ylab="",  xaxt="n", ann = FALSE, width=box_width, at=at.x)
  axis(1, at=at.x, labels=FALSE)
  text(x=at.x, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=pop_names_simple, srt=45, adj=1, xpd=TRUE, col=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3"))
  
  #mtext(side = 2, text = paste(units[i]), line = 2.2)
  

dev.off()
