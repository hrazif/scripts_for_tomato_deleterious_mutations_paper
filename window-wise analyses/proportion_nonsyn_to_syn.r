
#nonsyn=read.table("/home/caicedo/Hamid/20191118_deleterious_alleles/42_der_syn_nonsyn_noncoding_all_pops_per_ind_splittingSLL_bigInput/num_nonsynonymous_inds.txt", header=T, stringsAsFactors=F)

#syn=read.table("/home/caicedo/Hamid/20191118_deleterious_alleles/42_der_syn_nonsyn_noncoding_all_pops_per_ind_splittingSLL_bigInput/num_synonymous_inds.txt", header=T, stringsAsFactors=F)


#nonsyn$num_syn=syn$num_derived

#nonsyn$nonsyn_syn_ratio=nonsyn$num_derived/nonsyn$num_syn

nonsyn=read.table("nonsyn_syn_ratio.txt", header=T)

#kruskal test to test if populations have any effect on this ratio; if significant, the answer is yes
nonsyn$population <- as.factor(nonsyn$population)

kruskal.test(nonsyn$nonsyn_syn_ratio ~ nonsyn$population)


#install.packages("multcompView")
#install.packages("PMCMR")
require(multcompView)
require(PMCMR)

#dunn's test
out <- posthoc.kruskal.dunn.test(nonsyn$nonsyn_syn_ratio ~ nonsyn$population, p.adjust="bonf")


write.table (out$p.value, "dunn_test_results.txt", quote = F, sep = "\t")


#get p-values
out.p <- get.pvalues(out)
#get letters
out.mcV <- multcompLetters(out.p, threshold=0.05, reversed = T)

comp_letters=t(as.data.frame(out.mcV$Letters))


#sort
comp_letters=comp_letters[,pop_names]



#for plotting
pop_names=c("SP", "SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX_CA_NSA", "SLC_MEX", "SLL_Americas", "SLL_modern")

pop_size=vector()
for (j in pop_names) {
  temp=dim(nonsyn[nonsyn$population==j,])[1]
  pop_size=c(pop_size, temp)
}

box_width=sqrt(pop_size)
space=rep(0.1,length(pop_names))

pop_names_simple=names=c("SP", "SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX_CA_NSA", "SLC_MEX", "SLL_Americas", "SLL_modern")
#for combined pdfs

pdf("nonsyn_syn_ratio_inds.pdf", height = 6, width = 6)
      
  #to order by population
  nonsyn$population<-ordered(nonsyn$population, levels=pop_names)
  
  par(mfrow=c(1,1), mai=c(1.8, 0.9, 0.3, 0.05)*0.75)
  at.x=seq(0.1,0.8*length(pop_names),0.8)
  plot (nonsyn_syn_ratio~population, data=nonsyn, main="", col=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3"), varwidth=F, xlab="", ylab="",  xaxt="n", ann = FALSE, width=box_width, at=at.x, boxcol="red", whiskcol="red", medcol="black", staplecol="red", outcol="red", outpch="")
 
  stripchart(nonsyn_syn_ratio~population, data=nonsyn, vertical = T,  
             method = "jitter", add = TRUE, pch = 20, col = 'darkgray', at=at.x)
  
  
   axis(1, at=at.x, labels=FALSE)
  text(x=at.x, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=pop_names_simple, srt=45, adj=1, xpd=TRUE, col=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3"))
  
  #to add compact letters
  text(x=at.x, y=par()$usr[3]+1.03*(par()$usr[4]-par()$usr[3]), labels=as.character(comp_letters), xpd=TRUE)
  
  
  #mtext(side = 2, text = paste(units[i]), line = 2.2)
  

dev.off()


#write.table(nonsyn, "nonsyn_syn_ratio.txt", quote=F, sep="\t", col.names=T, row.names=F)