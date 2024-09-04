df=read.table("private_and_nonprivate_del_neutr.txt", header = T, stringsAsFactors = F)

library(reshape2)
library(ggplot2)

df_ = melt(df, id.vars=c("population"))

df_$population <- factor(df_$population,levels = c("SP_SECU", "SP_PER", "SP_NECU", "SLC_ECU", "SLC_PER", "SLC_San_Martin", "SLC_MEX-CA-NSA", "SLC_MEX", "SLL", "SLL_ex._Admixed"))

pdf("private_and_nonprivate_del_neutr.pdf", height = 4, width = 8)
ggplot(df_, aes(population, value, fill=variable)) +geom_bar(stat='Identity',position=position_dodge())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

