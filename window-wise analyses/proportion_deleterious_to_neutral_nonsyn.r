
#deleterious=read.table("/home/caicedo/Hamid/20191118_deleterious_alleles/43_deleterious_mutations_shared_by_SIFT_and_PROVEAN_splittinSLL_bigInput_wImp/num_deleterious_inds.txt", header=T, stringsAsFactors=F)

#neutral_nonsyn=read.table("/home/caicedo/Hamid/20191118_deleterious_alleles/45_num_neutral_nonsynonymous_based_on_SIFT_and_PROVEAN_splittingSLL_bigInput_wImp/num_neutral_nonsyn_inds.txt", header=T, stringsAsFactors=F)


#deleterious$num_neutral_nonsyn=neutral_nonsyn$num_derived

#deleterious$deleterious_neutral_nonsyn_ratio=deleterious$num_derived/deleterious$num_neutral_nonsyn



del_to_neutr_nosyn=read.table("deleterious_neutral_nonsyn_ratio.txt", header = T, stringsAsFactors = F, check.names = F)

pop_names=c("SP","SLC_ECU","SLC_PER","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX","SLL_Americas","SLL_modern")


#kruskal test to test if populations have any effect on this ratio; if significant, the answer is yes
del_to_neutr_nosyn$population <- as.factor(del_to_neutr_nosyn$population)

kruskal.test(del_to_neutr_nosyn$deleterious_neutral_nonsyn_ratio ~ del_to_neutr_nosyn$population)


#install.packages("multcompView")
#install.packages("PMCMR")
require(multcompView)
require(PMCMR)

#dunn's test
out <- posthoc.kruskal.dunn.test(del_to_neutr_nosyn$deleterious_neutral_nonsyn_ratio ~ del_to_neutr_nosyn$population, p.adjust="bonf")


write.table (out$p.value, "dunn_test_results.txt", quote = F, sep = "\t")


#get p-values
out.p <- get.pvalues(out)
#get letters
out.mcV <- multcompLetters(out.p, threshold=0.05, reversed = T)

comp_letters=t(as.data.frame(out.mcV$Letters))


#sort
comp_letters=comp_letters[,pop_names]


#for plotting

pop_size=vector()
for (j in pop_names) {
  temp=dim(del_to_neutr_nosyn[del_to_neutr_nosyn$population==j,])[1]
  pop_size=c(pop_size, temp)
}

box_width=sqrt(pop_size)
space=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)

pop_names_simple=c("SP","SLC ECU","SLC PER","SLC San Martin","SLC MEX-CA-NSA","SLC MEX","SLL Americas", "SLL modern")
#for combined pdfs
pdf("deleterious_neutral_nonsyn_ratio_inds_2.pdf", height = 6, width = 6)
      
  #to order by population
del_to_neutr_nosyn$population<-ordered(del_to_neutr_nosyn$population, levels=pop_names)
  
  par(mfrow=c(1,1), mai=c(1.8, 0.9,0.5, 0.05)*0.75)
  at.x=seq(0.1,length(pop_names)*0.8,0.8)
  plot (deleterious_neutral_nonsyn_ratio~population, data=del_to_neutr_nosyn, main="", col=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3"), varwidth=F, xlab="", ylab="",  xaxt="n", ann = FALSE, width=box_width, at=at.x,boxcol="red", whiskcol="red", medcol="black", staplecol="red", outcol="red", outpch="")
  axis(1, at=at.x, labels=FALSE)
 
  stripchart(deleterious_neutral_nonsyn_ratio~population, data=del_to_neutr_nosyn, vertical = T,  
             method = "jitter", add = TRUE, pch = 20, col = 'darkgray', at=at.x)
  
  
   text(x=at.x, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=pop_names_simple, srt=45, adj=1, xpd=TRUE, col=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3"))
  
  #to add compact letters
  text(x=at.x, y=par()$usr[3]+1.03*(par()$usr[4]-par()$usr[3]), labels=as.character(comp_letters), xpd=TRUE)
  
  #mtext(side = 2, text = paste(units[i]), line = 2.2)
  

dev.off()


#write.table(deleterious, "deleterious_neutral_nonsyn_ratio.txt", quote=F, sep="\t", col.names=T, row.names=F)




#install.packages("multcomp")

