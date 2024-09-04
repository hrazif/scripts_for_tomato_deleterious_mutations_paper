pop_colors=c("purple", "lightsalmon1", "darkgreen", "magenta", "orange", "blue", "red","firebrick3")

del=read.table("num_deleterious_inds.txt", header=TRUE)
colnames(del)[2]="num_deleterious"

syn=read.table("num_synonymous_inds.txt", header=TRUE)
syn=syn[,1:2]
colnames(syn)[2]="num_synonymous"


nonsyn=read.table("num_non_synonymous_inds.txt", header=TRUE)
nonsyn=nonsyn[,1:2]
colnames(nonsyn)[2]="num_nonsynonymous"


noncoding=read.table("num_noncoding_inds.txt", header=TRUE)
noncoding=noncoding[,1:2]
colnames(noncoding)[2]="num_noncoding"


neutral_nonsyn=read.table("num_neutral_nonsyn_inds.txt", header=TRUE)
neutral_nonsyn=neutral_nonsyn[,1:2]
colnames(neutral_nonsyn)[2]="num_neutral_nonsyn"



all_mutations=merge(del, syn, by="indivdual", all=TRUE)
all_mutations=merge(all_mutations, nonsyn, by="indivdual", all=TRUE)
all_mutations=merge(all_mutations, noncoding, by="indivdual", all=TRUE)
all_mutations=merge(all_mutations, neutral_nonsyn, by="indivdual", all=TRUE)



pi_s=read.table("mean_PIs_pops.txt", header=TRUE)

average_pi_mutations_per_pop=data.frame(pop=pi_s$pops, av_pi=pi_s$mean_PIs, av_del=0, av_syn=0, av_nonsyn=0,av_neutr_nonsyn=0, av_noncoding=0)

pops=pi_s$pops

for (i in 1:length(pops)){
	
	temp1=all_mutations[all_mutations$population==pops[i],]
	average_pi_mutations_per_pop[i,"av_del"]=mean(temp1$num_deleterious, na.rm=TRUE)
	average_pi_mutations_per_pop[i,"av_syn"]=mean(temp1$num_synonymous, na.rm=TRUE)
	average_pi_mutations_per_pop[i,"av_nonsyn"]=mean(temp1$num_synonymous, na.rm=TRUE)
	average_pi_mutations_per_pop[i,"av_noncoding"]=mean(temp1$num_noncoding, na.rm=TRUE)
	average_pi_mutations_per_pop[i,"av_neutr_nonsyn"]=mean(temp1$num_neutral_nonsyn, na.rm=TRUE)
	
}

#plotting correlation analyses

library(tidyverse)
library(ggpmisc)
library(ggplot2)




average_pi_mutations_per_pop$color=pop_colors
names(pop_colors)=average_pi_mutations_per_pop$pop

pdf("20220607_corr_av_pi_vs_av_del.pdf", height = 3, width = 3)


my.formula = y ~ x
ggplot(average_pi_mutations_per_pop, aes(x = av_pi, y = av_del, color = pop, group = 1)) +
  xlab("average pi")+ylab("average # deleterious mutations")+
  geom_point() + scale_color_manual(values = pop_colors) +
  geom_smooth(method="lm", formula = y ~ x ) + 
  stat_poly_eq() + 
  theme_classic() + theme(legend.position="none")

dev.off()






pdf("20220607_corr_av_pi_vs_av_syn.pdf", height = 3, width = 3)

my.formula = y ~ x
ggplot(average_pi_mutations_per_pop, aes(x = av_pi, y = av_syn, color = pop, group = 1)) +
  xlab("average pi")+ylab("average # synonymous mutations")+
  geom_point() + scale_color_manual(values = pop_colors) +
  geom_smooth(method="lm", formula = y ~ x ) + 
  stat_poly_eq() + 
  theme_classic() + theme(legend.position="none")


dev.off()


pdf("20220607_corr_av_pi_vs_av_nonsyn.pdf", height = 3, width = 3)

my.formula = y ~ x
ggplot(average_pi_mutations_per_pop, aes(x = av_pi, y = av_nonsyn, color = pop, group = 1)) +
  xlab("average pi")+ylab("average # nonsynonymous mutations")+
  geom_point() + scale_color_manual(values = pop_colors) +
  geom_smooth(method="lm", formula = y ~ x ) + 
  stat_poly_eq() + 
  theme_classic() + theme(legend.position="none")


dev.off()



pdf("20220607_corr_av_pi_vs_av_noncoding.pdf", height = 3, width = 3)

my.formula = y ~ x
ggplot(average_pi_mutations_per_pop, aes(x = av_pi, y = av_noncoding, color = pop, group = 1)) +
  xlab("average pi")+ylab("average # non-coding mutations")+
  geom_point() + scale_color_manual(values = pop_colors) +
  geom_smooth(method="lm", formula = y ~ x ) + 
  stat_poly_eq() + 
  theme_classic() + theme(legend.position="none")


dev.off()


#num neutral nonsynonymous

pdf("20220607_corr_av_pi_vs_av_num_neutr_nonsyn.pdf", height = 3, width = 3)

my.formula = y ~ x
ggplot(average_pi_mutations_per_pop, aes(x = av_pi, y = av_neutr_nonsyn, color = pop, group = 1)) +
  xlab("average pi")+ylab("average # neutral nonsynonymous mutations")+
  geom_point() + scale_color_manual(values = pop_colors) +
  geom_smooth(method="lm", formula = y ~ x ) + 
  stat_poly_eq() + 
  theme_classic() + theme(legend.position="none")


dev.off()


#jpeg("20220224_corr_stom_dens_gr_seas_precip.jpeg", res=300, height = 6, width = 6, units = "in")

#ggscatter(dataset_quant, x = "mean_growing_season_precipitation_mm", y = "average_stomatal_density_per_cm2",
#          add = "reg.line", conf.int = TRUE,
#          cor.coef = TRUE, cor.method = "pearson",
#          xlab = "median growing season precipitation (mm)", ylab = "stomatal density per cm2",
#          color = "purple")

#dev.off()
