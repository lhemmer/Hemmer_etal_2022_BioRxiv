#########
## Tajima's D and Pi for all the DGRP samples
#########

## load libraries

library(ggplot2)


#### load data for distributions of large chromosomes

#pwd <- "/Users/lucashemmer/Documents/mel_centromere/DGRP3/GATK/output/all/"

## Going with 500 bp windows

## 2L_1
l2.1 <- read.csv(paste(pwd,"2L_1/DGRP3.all.2L_1.500.csv",sep=""), header = T, stringsAsFactors=FALSE)
#l2.1 <- read.csv(paste(pwd,"2L_1/DGRP3.all.2L_1.1000.csv",sep=""), header = T, stringsAsFactors=FALSE)
#l2.1 <- read.csv(paste(pwd,"2L_1/DGRP3.all.2L_1.10000.csv",sep=""), header = T, stringsAsFactors=FALSE)

## 2R_21
r2.21 <- read.csv(paste(pwd,"2R_21/DGRP3.all.2R_21.500.csv",sep=""), header = T, stringsAsFactors=FALSE)
#r2.21 <- read.csv(paste(pwd,"2R_21/DGRP3.all.2R_21.1000.csv",sep=""), header = T, stringsAsFactors=FALSE)
#r2.21 <- read.csv(paste(pwd,"2R_21/DGRP3.all.2R_21.10000.csv",sep=""), header = T, stringsAsFactors=FALSE)

## 3L_1
l3.1 <- read.csv(paste(pwd,"3L_1/DGRP3.all.3L_1.500.csv",sep=""), header = T, stringsAsFactors=FALSE)
#l3.1 <- read.csv(paste(pwd,"3L_1/DGRP3.all.3L_1.1000.csv",sep=""), header = T, stringsAsFactors=FALSE)
#l3.1 <- read.csv(paste(pwd,"3L_1/DGRP3.all.3L_1.10000.csv",sep=""), header = T, stringsAsFactors=FALSE)

## 3R_28
r3.28 <- read.csv(paste(pwd,"3R_28/DGRP3.all.3R_28.500.csv",sep=""), header = T, stringsAsFactors=FALSE)
#r3.28 <- read.csv(paste(pwd,"3R_28/DGRP3.all.3R_28.1000.csv",sep=""), header = T, stringsAsFactors=FALSE)
#r3.28 <- read.csv(paste(pwd,"3R_28/DGRP3.all.3R_28.10000.csv",sep=""), header = T, stringsAsFactors=FALSE)

## X_1
x.1 <- read.csv(paste(pwd,"X_1/DGRP3.all.X_1.500.csv",sep=""), header = T, stringsAsFactors=FALSE)
#x.1 <- read.csv(paste(pwd,"X_1/DGRP3.all.X_1.1000.csv",sep=""), header = T, stringsAsFactors=FALSE)
#x.1 <- read.csv(paste(pwd,"X_1/DGRP3.all.X_1.10000.csv",sep=""), header = T, stringsAsFactors=FALSE)

## other for 4_2, centromeres
other <- read.csv(paste(pwd,"other/DGRP3.all.other.500.csv",sep=""), header = T, stringsAsFactors=FALSE)
#other <- read.csv(paste(pwd,"other/DGRP3.all.other.1000.csv",sep=""), header = T, stringsAsFactors=FALSE)
#other <- read.csv(paste(pwd,"other/DGRP3.all.other.10000.csv",sep=""), header = T, stringsAsFactors=FALSE)


#### combine data together

all.chrom <- rbind.data.frame(l2.1, r2.21, l3.1, r3.28, x.1, other[other$chrom %in% c("4_2","tig00057289", "3R_5", "Contig119", "Contig79"),])
all.chrom$chrom <- factor(all.chrom$chrom, levels = c("2L_1","tig00057289","2R_21","3L_1","3R_5","3R_28", "Contig119","4_2","Contig79","X_1"), ordered = T)


#### MWU tests

pi.mwu <- pairwise.wilcox.test(all.chrom$pi_DGRP, all.chrom$chrom, p.adjust.method = "bonferroni")
wilcox.test(all.chrom$pi_DGRP[all.chrom$chrom=="2L_1"],all.chrom$pi_DGRP[all.chrom$chrom=="tig00057289"])
wilcox.test(all.chrom$pi_DGRP[all.chrom$chrom=="2R_21"],all.chrom$pi_DGRP[all.chrom$chrom=="tig00057289"])
wilcox.test(all.chrom$pi_DGRP[all.chrom$chrom=="3L_1"],all.chrom$pi_DGRP[all.chrom$chrom=="3R_5"])
wilcox.test(all.chrom$pi_DGRP[all.chrom$chrom=="3R_28"],all.chrom$pi_DGRP[all.chrom$chrom=="3R_5"])
wilcox.test(all.chrom$pi_DGRP[all.chrom$chrom=="X_1"],all.chrom$pi_DGRP[all.chrom$chrom=="Contig79"])
wilcox.test(all.chrom$pi_DGRP[all.chrom$chrom=="4_2"],all.chrom$pi_DGRP[all.chrom$chrom=="Contig119"])



tmp <- data.frame(chrom = c("2L_1","tig00057289","3L_1","3R_5","Contig119","Contig79"), y.height = c(0.0192, 0.0182, 0.0182, 0.0172, 0.0092, 0.0132))

#### Figures

#pdf("/Users/lucashemmer/Documents/Figure_S7_all.boxplot.pi.pdf", width=8,height=6)
all.chrom %>% ggplot(aes(x=chrom, y=pi_DGRP), ylim=c(0,0.02)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,0.02) + 
  scale_x_discrete(labels=c("2L","Cen2","2R","3L","Cen3","3R","Cen4","4","CenX","X")) +
  theme_bw() + ylab("Pi") + xlab("Chromosome") +
  geom_text(data = tmp, aes(x = chrom, y = y.height, label = "***"), nudge_x = 0.5) +
  # 2L
  geom_segment(x = "2L_1", xend = "2L_1", y = 0.019, yend = 0.0185,colour = "black") +
  geom_segment(x = "tig00057289", xend = "tig00057289", y = 0.019, yend = 0.0185,colour = "black") +
  geom_segment(x = "2L_1", xend = "tig00057289", y = 0.019, yend = 0.019,colour = "black") +
  # 2R
  geom_segment(x = "tig00057289", xend = "tig00057289", y = 0.018, yend = 0.0175,colour = "black") +
  geom_segment(x = "2R_21", xend = "2R_21", y = 0.018, yend = 0.0175,colour = "black") +
  geom_segment(x = "2R_21", xend = "tig00057289", y = 0.018, yend = 0.018,colour = "black") +
  # 3L
  geom_segment(x = "3L_1", xend = "3L_1", y = 0.018, yend = 0.0175,colour = "black") +
  geom_segment(x = "3R_5", xend = "3R_5", y = 0.018, yend = 0.0175,colour = "black") +
  geom_segment(x = "3L_1", xend = "3R_5", y = 0.018, yend = 0.018,colour = "black") +
  # 3R
  geom_segment(x = "3R_28", xend = "3R_28", y = 0.017, yend = 0.0165,colour = "black") +
  geom_segment(x = "3R_5", xend = "3R_5", y = 0.017, yend = 0.0165,colour = "black") +
  geom_segment(x = "3R_28", xend = "3R_5", y = 0.017, yend = 0.017,colour = "black") + 
  # 4
  geom_segment(x = "4_2", xend = "4_2", y = 0.009, yend = 0.0085,colour = "black") +
  geom_segment(x = "Contig119", xend = "Contig119", y = 0.009, yend = 0.0085,colour = "black") +
  geom_segment(x = "4_2", xend = "Contig119", y = 0.009, yend = 0.009,colour = "black") + 
  # X
  geom_segment(x = "X_1", xend = "X_1", y = 0.013, yend = 0.0125,colour = "black") +
  geom_segment(x = "Contig79", xend = "Contig79", y = 0.013, yend = 0.0125,colour = "black") +
  geom_segment(x = "X_1", xend = "Contig79", y = 0.013, yend = 0.013,colour = "black") 
  


#dev.off()



tajD.mwu <- pairwise.wilcox.test(all.chrom$TajD_DGRP, all.chrom$chrom, p.adjust.method = "bonferroni")
wilcox.test(all.chrom$TajD_DGRP[all.chrom$chrom=="2L_1"],all.chrom$TajD_DGRP[all.chrom$chrom=="tig00057289"])
wilcox.test(all.chrom$TajD_DGRP[all.chrom$chrom=="2R_21"],all.chrom$TajD_DGRP[all.chrom$chrom=="tig00057289"])
wilcox.test(all.chrom$TajD_DGRP[all.chrom$chrom=="3L_1"],all.chrom$TajD_DGRP[all.chrom$chrom=="3R_5"])
wilcox.test(all.chrom$TajD_DGRP[all.chrom$chrom=="3R_28"],all.chrom$TajD_DGRP[all.chrom$chrom=="3R_5"])
wilcox.test(all.chrom$TajD_DGRP[all.chrom$chrom=="X_1"],all.chrom$TajD_DGRP[all.chrom$chrom=="Contig79"])
wilcox.test(all.chrom$TajD_DGRP[all.chrom$chrom=="4_2"],all.chrom$TajD_DGRP[all.chrom$chrom=="Contig119"])

tmp2 <- data.frame(chrom = c("3L_1","3R_5"), y.height = c(2.55, 2.35))

#pdf("/Users/lucashemmer/Documents/Figure_S8_all.boxplot.TajD.pdf", width=6,height=6)
all.chrom %>% ggplot(aes(x=chrom, y=TajD_DGRP)) +
  geom_boxplot(outlier.shape = NA) + ylim(-2.5,3) + 
  scale_x_discrete(labels=c("2L","Cen2","2R","3L","Cen3","3R","Cen4","4","CenX","X")) +
  theme_bw() + xlab("Chromosome") + ylab("Tajima's D") +
  geom_text(data = tmp2, aes(x = chrom, y = y.height, label = "*"), nudge_x = 0.5) +
  # 3L
  geom_segment(x = "3L_1", xend = "3L_1", y = 2.5, yend = 2.4,colour = "black") +
  geom_segment(x = "3R_5", xend = "3R_5", y = 2.5, yend = 2.4,colour = "black") +
  geom_segment(x = "3L_1", xend = "3R_5", y = 2.5, yend = 2.5,colour = "black") +
  # 3R
  geom_segment(x = "3R_28", xend = "3R_28", y = 2.3, yend = 2.2,colour = "black") +
  geom_segment(x = "3R_5", xend = "3R_5", y = 2.3, yend = 2.2,colour = "black") +
  geom_segment(x = "3R_28", xend = "3R_5", y = 2.3, yend = 2.3,colour = "black") 
  
dev.off()



sum.tab <- all.chrom %>% group_by(chrom) %>% dplyr::summarise(mean(pi_DGRP),median(pi_DGRP), mean(TajD_DGRP, na.rm=T), median(TajD_DGRP, na.rm=T))

write.table(sum.tab, file="/Users/lucashemmer/Documents/DGRP_pi_TajD_summary.txt", quote = F,sep="\t",row.names = F)
