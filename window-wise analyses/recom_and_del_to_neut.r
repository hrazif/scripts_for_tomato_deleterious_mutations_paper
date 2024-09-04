
pops=c("SP","SLC_ECU","SLC_PER","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX","SLL_Americas","SLL_modern")



pdf("recom_and_prop_del_neutr.pdf", width = 4, height = 4)

wilcoxon_test_output=data.frame(pop=pops, p_value=0, mean_prop_del_neutr_high_recom=0, mean_prop_del_neutr_low_recom=0,median_prop_del_neutr_high_recom=0, median_prop_del_neutr_low_recom=0)


recom_file=read.table("/project/uma_ana_caicedo/Hamid/20191118_deleterious_alleles/8_recombination_rate_100kb/recom_rate_all_chr.txt", header = T, stringsAsFactors = F)
  
for(i in 1: length(pops)){

 

  #deleterious
  
  deleterious_file=read.table(paste0("/project/uma_ana_caicedo/Hamid/20191118_deleterious_alleles/72_num_deleterious_100kb_bigInput_wImp/",pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), header = T, stringsAsFactors = F)

  
  #neutral_nonsyn
  neutr_nonsyn_file=read.table(paste0("/project/uma_ana_caicedo/Hamid/20191118_deleterious_alleles/73_num_neutral_nonsyn_100kb_bigInput_wImp/",pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), header = T, stringsAsFactors = F)

  
  
  #to merge
  
  recom_del=merge(recom_file, deleterious_file, by=c("chromosome","windows"), all = T)
  
  recom_del_neutr=merge(recom_del, neutr_nonsyn_file, by=c("chromosome","windows"), all = T)
  
  colnames(recom_del_neutr)[4:5]=c("num_del","num_neutr_nonsyn")
  

  
  #to remove windows with recombination rate=NA
  recom_del_neutr=recom_del_neutr[!is.na(recom_del_neutr$recome_rate),]
  
  #to remove windows with inf remobination rate
  recom_del_neutr=recom_del_neutr[!is.infinite(recom_del_neutr$recome_rate),]
 
  #prop del to neutr
  recom_del_neutr$prop_del_neutr=recom_del_neutr$num_del/recom_del_neutr$num_neutr_nonsyn
  
  
  #to remove NaN (0/0) values for prop del to neutr
  recom_del_neutr=recom_del_neutr[!is.nan(recom_del_neutr$prop_del_neutr),]
  
  
  #to avoid Infs
  recom_del_neutr$num_neutr_nonsyn[recom_del_neutr$num_neutr_nonsyn==0]=1
  
   
  #to recalculate prop del to neutr
  recom_del_neutr$prop_del_neutr=recom_del_neutr$num_del/recom_del_neutr$num_neutr_nonsyn
  
  
  #to order by recombination rate
  recom_del_neutr=recom_del_neutr[order(recom_del_neutr$recome_rate, decreasing = T),]
  
  
  recom_del_neutr$type="low recom"
  
  top_windows=round(nrow(recom_del_neutr)/2)
  recom_del_neutr[1:top_windows,"type"]="high recom"
  
  
  recom_del_neutr$type=ordered(recom_del_neutr$type, c("low recom", "high recom"))
  
  boxplot(recom_del_neutr$prop_del_neutr ~ recom_del_neutr$type, main=pops[i], xlab = "recombination rate", ylab = "Prop (del/neutr)")
  
  
  #wilcoxon test
  
  temp1=recom_del_neutr[recom_del_neutr$type=="low recom",]
  
  temp2=recom_del_neutr[recom_del_neutr$type=="high recom",]
  
  wilcoxon_test=wilcox.test(temp1$prop_del_neutr, temp2$prop_del_neutr, paired = F)
  
  wilcoxon_test_output$p_value[i]= wilcoxon_test$p.value 
  
  
  
  #means and medians
  
  wilcoxon_test_output$mean_prop_del_neutr_high_recom[i]=mean(temp2$prop_del_neutr, na.rm = T)
  
  wilcoxon_test_output$mean_prop_del_neutr_low_recom[i]=mean(temp1$prop_del_neutr, na.rm = T)
  
  wilcoxon_test_output$median_prop_del_neutr_high_recom[i]=median(temp2$prop_del_neutr, na.rm = T)
  
  wilcoxon_test_output$median_prop_del_neutr_low_recom[i]=median(temp1$prop_del_neutr, na.rm = T)
  
  
  

}
dev.off()

wilcoxon_test_output$significance="no"

wilcoxon_test_output[wilcoxon_test_output$p_value < 0.05/length(pops), "significance"]="yes"

write.table(wilcoxon_test_output, "wilcoxon_test_output.txt", col.names = T, row.names = F, sep="\t", quote = F)

wilcoxon_test_output

