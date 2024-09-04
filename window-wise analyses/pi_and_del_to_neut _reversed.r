setwd("C:/BackUps/Analyses/20191118_deleterious_alleles/83_correlation_Pi_prop_del_neutr_100kb")

pops=c("SLC_ECU","SLC_San_Martin","SLC_MEX_CA_NSA","SLC_MEX","SLL_Americas","SLL_modern")



pdf("pi_and_prop_del_neutr_reversed.pdf", width = 4, height = 4)

wilcoxon_test_output=data.frame(pop=pops, p_value=0, mean_prop_del_neutr_high_PI=0, mean_prop_del_neutr_low_PI=0,median_prop_del_neutr_high_PI=0, median_prop_del_neutr_low_PI=0)


for(i in 1: length(pops)){

  pi_file=read.table(paste0("C:/BackUps/Analyses/20200406/4_sweeps_w_Imputation_100kb/2_pi_ratio/",pops[i], ".windowed.pi"), header = T, stringsAsFactors = F)
  
  colnames(pi_file)[1]="chromosome"
  pi_file$chromosome=gsub("SL2.50ch","", gsub("SL2.50ch0","", pi_file$chromosome))
  

  #deleterious
  
  deleterious_file=read.table(paste0("C:/BackUps/Analyses/20191118_deleterious_alleles/72_num_deleterious_100kb_bigInput_wImp/",pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), header = T, stringsAsFactors = F)
  
  deleterious_file$windows=(deleterious_file$windows-1)*100000+1
  
  colnames(deleterious_file)[2]="BIN_START"
  
  
  
  #neutral_nonsyn
  neutr_nonsyn_file=read.table(paste0("C:/BackUps/Analyses/20191118_deleterious_alleles/73_num_neutral_nonsyn_100kb_bigInput_wImp/",pops[i],"_num_sites_w_derived_alleles_all_chr.txt"), header = T, stringsAsFactors = F)
  
  neutr_nonsyn_file$windows=(neutr_nonsyn_file$windows-1)*100000+1
  
  colnames(neutr_nonsyn_file)[2]="BIN_START"
  

  
  
  #to merge
  
  pi_del=merge(pi_file, deleterious_file, by=c("chromosome","BIN_START"), all = T)
  
  pi_del_neutr=merge(pi_del, neutr_nonsyn_file, by=c("chromosome","BIN_START"), all = T)
  
  colnames(pi_del_neutr)[6:7]=c("num_del","num_neutr_nonsyn")
  

  
  #to remove windows with pi=NA
  pi_del_neutr=pi_del_neutr[!is.na(pi_del_neutr$PI),]
  
  

  #prop del to neutr
  pi_del_neutr$prop_del_neutr=pi_del_neutr$num_del/pi_del_neutr$num_neutr_nonsyn
  
  
  #to remove NaN (0/0) values for prop del to neutr
  pi_del_neutr=pi_del_neutr[!is.nan(pi_del_neutr$prop_del_neutr),]
  
  
  #to avoid Infs
  pi_del_neutr$num_neutr_nonsyn[pi_del_neutr$num_neutr_nonsyn==0]=1
  
  
  #to recalculate prop del to neutr
  pi_del_neutr$prop_del_neutr=pi_del_neutr$num_del/pi_del_neutr$num_neutr_nonsyn
  
  
  #to order by PI
  pi_del_neutr=pi_del_neutr[order(pi_del_neutr$PI, decreasing = T),]
  
  
  pi_del_neutr$type="low PI"
  
  top_windows=round(nrow(pi_del_neutr)/2)
  pi_del_neutr[1:top_windows,"type"]="high PI"
  
  
  pi_del_neutr$type=ordered(pi_del_neutr$type, c("low PI", "high PI"))
  
  boxplot(pi_del_neutr$prop_del_neutr ~ pi_del_neutr$type, main=pops[i], xlab = "PI", ylab = "Prop (del/neutr)")
  
  
  #wilcoxon test
  
  temp1=pi_del_neutr[pi_del_neutr$type=="low PI",]
  
  temp2=pi_del_neutr[pi_del_neutr$type=="high PI",]
  
  wilcoxon_test=wilcox.test(temp1$prop_del_neutr, temp2$prop_del_neutr, paired = F)
  
  wilcoxon_test_output$p_value[i]= wilcoxon_test$p.value 
  
  
  
  #means and medians
  
  wilcoxon_test_output$mean_prop_del_neutr_high_PI[i]=mean(temp2$prop_del_neutr, na.rm = T)
  
  wilcoxon_test_output$mean_prop_del_neutr_low_PI[i]=mean(temp1$prop_del_neutr, na.rm = T)
  
  wilcoxon_test_output$median_prop_del_neutr_high_PI[i]=median(temp2$prop_del_neutr, na.rm = T)
  
  wilcoxon_test_output$median_prop_del_neutr_low_PI[i]=median(temp1$prop_del_neutr, na.rm = T)
  
  
  

}
dev.off()

wilcoxon_test_output$significance="no"

wilcoxon_test_output[wilcoxon_test_output$p_value < 0.05/length(pops), "significance"]="yes"

write.table(wilcoxon_test_output, "wilcoxon_test_output_reversed.txt", col.names = T, row.names = F, sep="\t", quote = F)

wilcoxon_test_output

