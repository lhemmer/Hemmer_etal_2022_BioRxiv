########
#### Pairwise analysis of G2/Jockey-3 in melanogaster and figures
########

#### load libraries

library(ape)
library(phytools)
library(geiger)
library(evobiR)
library(ggnewscale)
library(ggtree)
library(dplyr)
library(tidyr)
library(ggpubr)
options(scipen=999)


#### load data 

## pairwise

pwd <- "/Users/lucashemmer/Documents/mel_centromere/ReAnnotation/phylogeny/Dmel_Ref_Denovo_aligntoRef_Y_chrom/"
g2.mat1 <- read.csv("/Users/lucashemmer/Documents/mel_centromere/ReAnnotation/GFF/melanogaster/082020/dmel_scaffold2_plus0310_2.fasta.out.G2.pairwise1.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)
g2.mat2 <- read.csv("/Users/lucashemmer/Documents/mel_centromere/ReAnnotation/GFF/melanogaster/082020/dmel_scaffold2_plus0310_2.fasta.out.G2.pairwise2.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)

## info
ref.calls <- read.table("/Users/lucashemmer/Documents/mel_centromere/ReAnnotation/GFF/melanogaster/082020/dmel_scaffold2_plus0310_2.fasta.out.G2.info.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)


#### process

ref.calls <- ref.calls[!is.na(ref.calls$ref.chromatin),]
ref.calls <- ref.calls[!is.na(ref.calls$uniq.subs),]
ref.calls <- ref.calls %>% mutate(centromere2=ifelse(ref.region=="centromere" & ref.chrom!="Y","Non-Y Cen", ifelse(ref.region=="centromere" & ref.chrom=="Y","CenY", ifelse(ref.region!="centromere" & ref.chromatin=="het","Heterochromatin","Euchromatin"))))

ref.calls$centromere2 <- factor(ref.calls$centromere2, levels=c("Euchromatin","Heterochromatin","Non-Y Cen", "CenY"))

ref.calls <- ref.calls %>% mutate(contig2=ifelse(centromere2=="CenY", grouping, contig))

ref.calls <- ref.calls %>% mutate(centromere3=ifelse(ref.region=="centromere", "Centromere", ifelse(ref.region!="centromere" & ref.chromatin=="het","Heterochromatin","Euchromatin")))


#### formatting data for use

g2.mat <- g2.mat1

for (i in 1:nrow(g2.mat)) {
  for (j in 1:ncol(g2.mat)) {
    if (i ==j) {
      g2.mat[i,j] <- 0.0
    }
    else if (is.na(g2.mat[i,j])) {
      g2.mat[i,j] <- g2.mat2[i,j]
    } 
  }
}



ncol(g2.mat)
#g2.df <- g2.mat
#g2.mat <- g2.df

g2.mat <- g2.mat[colnames(g2.mat) != "3R_5|68215-68413|+|", row.names(g2.mat) != "3R_5|68215-68413|+|"]
g2.mat <- g2.mat[colnames(g2.mat) != "3R_5|68414-68544|-|", row.names(g2.mat) != "3R_5|68414-68544|-|"]

g2.mat1 <- g2.mat1[colnames(g2.mat1) != "3R_5|68215-68413|+|", row.names(g2.mat1) != "3R_5|68215-68413|+|"]
g2.mat1 <- g2.mat1[colnames(g2.mat1) != "3R_5|68414-68544|-|", row.names(g2.mat1) != "3R_5|68414-68544|-|"]


g2.df <- (g2.mat1)
g2.df$te.name <- row.names(g2.df)

g2.df1 <- left_join(g2.df, ref.calls %>% filter(!is.na(uniq.subs )))


#clades <- unique(ref.calls$centromere2)
#clades <- unique(ref.calls$centromere3)

clades <- unique(ref.calls$grouping)
clade.list <- list()

for (i in 1:length(clades)) {
  clade1 <- ref.calls %>% filter(grouping==clades[i])
  g2.names1 <- clade1$te.name
  g2.df.plot <- g2.df[row.names(g2.df) %in% clade1$te.name, colnames(g2.df) %in% clade1$te.name]
  #MyHeatmap <- LDheatmap(g2.df.plot, genetic.distances = g2.names1$mid.var1, color=cols, add.key=TRUE, flip=TRUE)
  g2.pair <- reshape2::melt(base::as.matrix(g2.df.plot))
  g2.pair <- g2.pair[!is.na(g2.pair$value),]
  g2.pair <- g2.pair[!duplicated(g2.pair), ]
  clade.list[[i]] <- g2.pair
  clade.list[[i]]$clade <- clades[i]
}



clade.dat <- do.call(rbind, clade.list)
clade.dat$clade <- factor(clade.dat$clade , levels=c("Euchromatin","Heterochromatin","Non-Y Cent", "Clade1", "Clade2", "Clade3", "Clade4", "Clade4A"))


clade.dat1 <- clade.dat %>% filter(value < 0.4) ## none are comparing within the same clade

clade.dat1$clade <- factor(clade.dat1$clade , levels=c("Euchromatin","Heterochromatin","Non-Y Cent", "Clade1", "Clade2", "Clade3", "Clade4", "Clade4A"))


clade.pvalue <- compare_means(value ~ clade,  data = clade.dat1)

write.csv(clade.pvalue, file="/Users/lucashemmer/Documents/DGRP3_pairwise_clade_MWU.csv", quote = F, row.names = F)

set.seed(2023)
p5 <- clade.dat1 %>% ggplot(aes(x=clade, y = value, color=clade)) + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_violin(outlier.shape = NA) +
  geom_jitter(height = 0) + 
  #ylim(0,0.4) + 
  theme_bw() + ylab("Genetic Divegence") + xlab("Grouping") +
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#FFEA46FF", "#7C7B78FF", "#00204DFF","#00204DFF", "#00204DFF","#00204DFF", "#00204DFF","#00204DFF")) +
  scale_x_discrete(labels = c("Euchromatin","Heterochromatin","Non-Y Cent", "Class 1", "Class 2", "Class 3", "Class 4", "Class 4A")) #+
  #stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                   ref.group = "NULL")


div.anova <- aov(value ~ clade, data = clade.dat)
TukeyHSD(div.anova)

tajD.mwu <- pairwise.wilcox.test(clade.dat$value, clade.dat$clade, p.adjust.method = "bonferroni")
tajD.mwu <- pairwise.wilcox.test(clade.dat$value, clade.dat$clade)

#### breaking it down by contig now


clades <- c("2L_1", "3L_1","Contig135","Contig5","tig00022795","Y_Contig10","Y_Contig20","Y_Contig70","Y_scaffold4","Y_scaffold5","Y_scaffold7","Contig79","3R_5","Contig119","Clade1","Clade2","Clade3","Clade4","Clade4A")
clade.list <- list()

for (i in 1:length(clades)) {
  clade1 <- ref.calls %>% filter(contig2==clades[i])
  g2.names1 <- clade1$te.name
  g2.df.plot <- g2.df[row.names(g2.df) %in% clade1$te.name, colnames(g2.df) %in% clade1$te.name]
  #MyHeatmap <- LDheatmap(g2.df.plot, genetic.distances = g2.names1$mid.var1, color=cols, add.key=TRUE, flip=TRUE)
  g2.pair <- reshape2::melt(base::as.matrix(g2.df.plot))
  g2.pair <- g2.pair[!is.na(g2.pair$value),]
  g2.pair <- g2.pair[!duplicated(g2.pair), ]
  clade.list[[i]] <- g2.pair
  clade.list[[i]]$clade <- clades[i]
}


clade.dat <- do.call(rbind, clade.list)
clade.dat1 <- clade.dat %>% filter(value < 0.4) ## none are comparing within the same clade

clade.dat1$clade <- factor(clade.dat1$clade , levels=c("2L_1", "3L_1","Contig135","Contig5","tig00022795","Y_Contig10","Y_Contig20","Y_Contig70","Y_scaffold4","Y_scaffold5","Y_scaffold7","Contig79","3R_5","Contig119","Clade1","Clade2","Clade3","Clade4","Clade4A"))


clade.pvalue <- compare_means(value ~ clade,  data = clade.dat1)
write.csv(clade.pvalue, file="/Users/lucashemmer/Documents/DGRP3_pairwise_contig_MWU.csv", quote = F, row.names = F)

set.seed(2023)
p6 <- clade.dat1 %>% ggplot(aes(x=clade, y = value, color=clade)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(height = 0) + 
  #ylim(0,0.4) + 
  theme_bw() + ylab("Genetic Divegence") + xlab("Grouping") +
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#FFEA46FF", "#FFEA46FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF",  "#7C7B78FF", "#00204DFF","#00204DFF", "#00204DFF","#00204DFF", "#00204DFF","#00204DFF", "#00204DFF", "#00204DFF")) +
  scale_x_discrete(labels = c("2L_1", "3L_1","Contig135","Contig5","tig00022795","Y_Contig10","Y_Contig20","Y_Contig70","Y_scaffold4","Y_scaffold5","Y_scaffold7","Contig79","3R_5","Contig119","Class1","Class2","Class3","Class4","Class4A")) #+
 #+
#stat_compare_means(label = "p.signif", method = "wilcox.test",
#                   ref.group = "NULL")



#pdf("/Users/lucashemmer/Documents/Figure_S8_Jockey-3_DmelRef_intraRegion_intraContig_pairwise.pdf",width=16,height=8)
ggarrange(p5,p6,nrow=2,ncol=1,labels=c("A","B"), heights = c(1,1))
dev.off()










### break glass in case of emergency

for (i in 1:length(clades)) {
  clade1 <- ref.calls %>% filter(centromere2==clades[i])
  g2.names1 <- clade1$te.name
  g2.df.plot <- g2.df[row.names(g2.df) %in% clade1$te.name, colnames(g2.df) %in% clade1$te.name]
  #MyHeatmap <- LDheatmap(g2.df.plot, genetic.distances = g2.names1$mid.var1, color=cols, add.key=TRUE, flip=TRUE)
  g2.pair <- reshape2::melt(base::as.matrix(g2.df.plot))
  g2.pair <- g2.pair[!is.na(g2.pair$value),]
  g2.pair <- g2.pair[!duplicated(g2.pair), ]
  clade.list[[i]] <- g2.pair
  clade.list[[i]]$clade <- clades[i]
}

for (i in 1:length(clades)) {
  clade1 <- ref.calls %>% filter(centromere3==clades[i])
  g2.names1 <- clade1$te.name
  g2.df.plot <- g2.df[row.names(g2.df) %in% clade1$te.name, colnames(g2.df) %in% clade1$te.name]
  #MyHeatmap <- LDheatmap(g2.df.plot, genetic.distances = g2.names1$mid.var1, color=cols, add.key=TRUE, flip=TRUE)
  g2.pair <- reshape2::melt(base::as.matrix(g2.df.plot))
  g2.pair <- g2.pair[!is.na(g2.pair$value),]
  g2.pair <- g2.pair[!duplicated(g2.pair), ]
  clade.list[[i]] <- g2.pair
  clade.list[[i]]$clade <- clades[i]
}

for (i in 1:length(clades)) {
  clade1 <- ref.calls %>% filter(contig2==clades[i])
  g2.names1 <- clade1$te.name
  g2.df.plot <- g2.df[row.names(g2.df) %in% clade1$te.name, colnames(g2.df) %in% clade1$te.name]
  #MyHeatmap <- LDheatmap(g2.df.plot, genetic.distances = g2.names1$mid.var1, color=cols, add.key=TRUE, flip=TRUE)
  g2.pair <- reshape2::melt(base::as.matrix(g2.df.plot))
  g2.pair <- g2.pair[!is.na(g2.pair$value),]
  g2.pair <- g2.pair[!duplicated(g2.pair), ]
  clade.list[[i]] <- g2.pair
  clade.list[[i]]$clade <- clades[i]
}
