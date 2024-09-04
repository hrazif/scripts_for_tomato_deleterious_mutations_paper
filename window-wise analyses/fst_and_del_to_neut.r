setwd("C:/BackUps/Analyses/20191118_deleterious_alleles/82_correlation_Fst_prop_del_neutr_100kb")

pops=c("SLC_ECU","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX","SLL_Americas","SLL_modern")

comparisons=c("SP_to_SLC_ECU",
       "SLC_PER_to_SLC_San_Martin",
       "SLC_ECU_PER_SM_to_SLC_MEX_CA_NSA",
       "SLC_ECU_PER_SM_to_SLC_MEX",
       "SLC_MEX_to_SLL_Americas",
       "SLC_MEX_to_SLL_modern")

pdf("Fst_and_prop_del_neutr.pdf", width = 4, height = 4)


wilcoxon_test_output=data.frame(pop=pops, p_value=0, mean_fst_high_load=0, mean_fst_low_load=0,median_fst_high_load=0, median_fst_low_load=0)

for(i in 1: length(pops)){

  fst_file=read.table(paste0("C:/BackUps/Analyses/20200406/4_sweeps_w_Imputation_100kb/1_Fst/",comparisons[i], ".windowed.weir.fst"), header = T, stringsAsFactors = F)
  
  fst_file$CHROM=gsub("SL2.50ch","", gsub("SL2.50ch0","", fst_file$CHROM))
  
  colnames(fst_file)[1]="chromosome"
  
  #deleterious
  
  deleterious_file=read.table(paste0("C:/BackUps/Analyses/20191118_deleterious_alleles/72_num_deleterious_100kb_bigInput_wImp/",pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), header = T, stringsAsFactors = F)
  
  deleterious_file$windows=(deleterious_file$windows-1)*100000+1
  
  colnames(deleterious_file)[2]="BIN_START"
  
  
  
  #neutral_nonsyn
  neutr_nonsyn_file=read.table(paste0("C:/BackUps/Analyses/20191118_deleterious_alleles/73_num_neutral_nonsyn_100kb_bigInput_wImp/",pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), header = T, stringsAsFactors = F)
  
  neutr_nonsyn_file$windows=(neutr_nonsyn_file$windows-1)*100000+1
  
  colnames(neutr_nonsyn_file)[2]="BIN_START"
  
  
  #to merge
  
  fst_del=merge(fst_file, deleterious_file, by=c("chromosome","BIN_START"), all = T)
  
  fst_del_neutr=merge(fst_del, neutr_nonsyn_file, by=c("chromosome","BIN_START"), all = T)
  
  colnames(fst_del_neutr)[7:8]=c("num_del","num_neutr_nonsyn")
  
  
  #to fill in Nas for del and neutr
  fst_del_neutr$num_del[is.na( fst_del_neutr$num_del)]=0
  
  fst_del_neutr$num_neutr_nonsyn[is.na( fst_del_neutr$num_neutr_nonsyn)]=0
  
  
  #to remove windows with Fst=NA
  fst_del_neutr=fst_del_neutr[!is.na(fst_del_neutr$MEAN_FST),]
  
  
  #to avoid Infs
  fst_del_neutr$num_neutr_nonsyn[fst_del_neutr$num_neutr_nonsyn==0]=1
  
  #prop del to neutr
  fst_del_neutr$prop_del_neutr=fst_del_neutr$num_del/fst_del_neutr$num_neutr_nonsyn
  
  
  #to order by Prop del to neutr
  fst_del_neutr=fst_del_neutr[order(fst_del_neutr$prop_del_neutr, decreasing = T),]
  
  
  fst_del_neutr$type="low load"
  
  top_windows=0.05*nrow(fst_del_neutr)
  fst_del_neutr[1:top_windows,"type"]="high load"
  
  
  fst_del_neutr$type=ordered(fst_del_neutr$type, c("low load", "high load"))
  
  boxplot(fst_del_neutr$MEAN_FST ~ fst_del_neutr$type, main=comparisons[i], xlab = "", ylab = "Fst")
  
  
  #wilcoxon test
  
  temp1=fst_del_neutr[fst_del_neutr$type=="low load",]
  
  temp2=fst_del_neutr[fst_del_neutr$type=="high load",]
  
  wilcoxon_test=wilcox.test(temp1$MEAN_FST, temp2$MEAN_FST, paired = F)
  
  wilcoxon_test_output$p_value[i]= wilcoxon_test$p.value 
  
  
  
  #means and medians
  
  wilcoxon_test_output$mean_fst_high_load[i]=mean(temp2$MEAN_FST, na.rm = T)
  
  wilcoxon_test_output$mean_fst_low_load[i]=mean(temp1$MEAN_FST, na.rm = T)
  
  wilcoxon_test_output$median_fst_high_load[i]=median(temp2$MEAN_FST, na.rm = T)
  
  wilcoxon_test_output$median_fst_low_load[i]=median(temp1$MEAN_FST, na.rm = T)
  
  

}
dev.off()

wilcoxon_test_output$significance="no"

wilcoxon_test_output[wilcoxon_test_output$p_value < 0.05/length(pops), "significance"]="yes"

write.table(wilcoxon_test_output, "wilcoxon_test_output.txt", col.names = T, row.names = F, sep="\t", quote = F)

wilcoxon_test_output
