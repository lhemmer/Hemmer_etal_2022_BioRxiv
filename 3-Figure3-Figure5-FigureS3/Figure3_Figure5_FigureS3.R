######## 
#### Analysis of Distribution of Transposable Elements
######## 

#### Load libraries

library(DescTools)
library(viridis)
library(ggplot2)
library(tidyr)
library(dplyr)
options(scipen=999)

#### Load data

DGRP.all <- read.table("File_S2_McClintock_output_DGRP.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DGRP.all$ID[DGRP.all$ID=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"


## higher coverage sample calls
DGRP.cent.cov <- read.table("Table_S3_DGRP_total_cent_cov_500bpWin.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
cent.hcov.pop <- DGRP.cent.cov %>% filter(centromere=="cent" & med_cov > 20 & population!=340 & population!=373) %>% pull(population) 




######################################################
#### De novo calls need to either be in more than 1 pop or confirmed by multiple methods
######################################################



DGRP.all.filt <- DGRP.all %>% mutate(centromere=ifelse(region=="centromere","Centromere",
  ifelse(region!="centromere" & chromatin.state=="het","Heterochromatin","Euchromatin"))) %>% ## filter out piRNA clusters for popgen analysis
  filter(!(contig=="2L_1" & start < 20270000 & start > 20140000)) %>% ## 38C1/C2
  filter(!(contig=="3L_1" & start < 23310000 & start > 23260000)) %>% ## 80F
  filter(!(contig=="2R_21" & start < 2780000 & start > 2520000)) %>% ## 42AB
  filter(!(contig=="X_1" & start < 21880000 & start > 21631000)) %>% ## flamenco
  filter(!(contig=="X_1" & start < 21560000 & start > 21500000)) %>% ## 20A
  filter(!(chromosome=="Y")) %>% #%>% filter(ref=="non-reference") ## Y chromosome calls unrealiable in DGRP
  filter(!(ID %in% c("Gypsy6_I_Dmoj","Gypsy1_I_Dmoj"))) %>% ## repeats misannotated as Gypsy elements from D. mojavensis
  filter(!(population %in% c(138, 340))) %>% ## much higher TE calls than other samples, could lead to too many false positives
  filter(population %in% cent.hcov.pop) %>% 
  filter(identifier!="Jockey_3_Dmel_08212020.2R_19.0")

#### reconfigure pop counts to account for those taken out by removing 138 and 340

## re-adjust population counts for our sample size, just takes a few seconds
DGRP.all.filt$pop.count <- NA
DGRP.all.filt <- split(DGRP.all.filt,f=DGRP.all.filt$identifier)
for (k in 1:length(DGRP.all.filt)) {
  tmp <- length(unique(DGRP.all.filt[[k]]$population))
  DGRP.all.filt[[k]]$pop.count <- tmp
}
DGRP.all.filt <- do.call(rbind, DGRP.all.filt)
rownames(DGRP.all.filt) <- 1:nrow(DGRP.all.filt)

## Number of population samples
pop.num <- length(unique(DGRP.all.filt$population))


## Remove repeat calls for same TE insertions if confirmed by multiple methods
nrow(DGRP.all.filt)
DGRP.all.filt1 <- DGRP.all.filt[!duplicated(DGRP.all.filt[,c("mcc.identifier")]),]
nrow(DGRP.all.filt1)


## IMPORTANT filter out non-reference calls if not found in more than one population or confirmed by multiple methods
DGRP.all.filt2 <- DGRP.all.filt1 %>% filter(!(ref=="non-reference" & call.count==1 & pop.count==1))
nrow(DGRP.all.filt2)


## Keeping all unique insertions in entire population 
DGRP.pop <- DGRP.all.filt2[!duplicated(DGRP.all.filt2[,c("identifier")]),]
nrow(DGRP.pop)
count.by.region <- as.data.frame(DGRP.pop %>% mutate(centromere=ifelse(region=="centromere","Centromere",
  ifelse(region!="centromere" & chromatin.state=="het","Heterochromatin","Euchromatin"))) %>% 
  count(ref,centromere) %>% ungroup())

## Categorize insertions based on frequency, rare < 10%, Intermediate 10-50%, Common > 50%, also add easy region
DGRP.pop.fig <- as.data.frame(
  DGRP.pop %>% mutate(centromere=ifelse(region=="centromere","Centromere",
  ifelse(region!="centromere" & chromatin.state=="het","Heterochromatin","Euchromatin"))) %>% #filter(ID %in% xxx$ID) %>%
    mutate(freq.class = ifelse(pop.count <= (pop.num* 0.1), "Rare",
    ifelse(pop.count > (pop.num* 0.1) & pop.count < (pop.num* 0.5), "Intermediate","Common")))
)



#### distributon of TEs by region

dmel.tes <- read.table("uniq.DGRP.tes.mel.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DGRP.pop.fig.alt <- DGRP.pop.fig %>% filter(ID %in% dmel.tes$TE) 


count.tots <- as.data.frame(DGRP.pop.fig.alt %>% filter(freq.class=="Singleton") %>% count(ID, name = "num") %>% arrange(desc(num)) %>% ungroup()) 
count.tots <- count.tots %>% mutate(ID = str_replace(ID, "_I", ""), ID = str_replace(ID, "-I", ""), ID = str_replace(ID, "_LTR", ""), ID = str_replace(ID, "-LTR", "")) %>% distinct(ID, .keep_all=T)



DGRP.pop.fig.alt$centromere <- factor(DGRP.pop.fig.alt$centromere, levels=c("Centromere","Heterochromatin","Euchromatin"))
DGRP.pop.fig.alt$freq.class <- factor(DGRP.pop.fig.alt$freq.class, levels=c("Rare","Intermediate","Common"))
DGRP.pop.fig.alt$ref <- factor(DGRP.pop.fig.alt$ref, levels=c("reference","non-reference"))


#pdf("/Users/lucashemmer/Documents/S15_DGRP3_totalTE_count_dist_byRegion.pdf",width=6,height=10)
ggplot(DGRP.pop.fig.alt, aes(x=freq.class)) +
  geom_bar(aes(fill=freq.class),position="dodge", stat="count") +
  ggtitle("Transposable Element Count per Region") + 
  facet_grid(rows = vars(centromere),scales = "free_y") + 
  theme_bw() + ylab("Count") + xlab("Frequency Class") + theme(legend.position = "none") +
  #theme(strip.background =element_rect(fill="white")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(discrete = T, option = "D", direction = 1)
#dev.off()


#### Distribution of centromere TEs

## select group of TEs
#original #cent.te.to.keep <- c("Copia_I","DMCR1A","DOC","DOC2_DM","DOC6_DM","FW_DM","G5_DM","I_DM","Jockey_3_DM","LINEJ1_DM","NOMAD_I","PROTOP_A","R1_DM","ROO_I")
cent.te.to.keep <- c("BS","DOC","DOC2_DM","DOC6_DM","FW_DM","G5_DM","Jockey_3_DM","LINEJ1_DM","NOMAD_I","PROTOP_A","R1_DM","ROO_I","Jockey_3_DM")
#Super strict $cent.te.to.keep <- c("DOC","Jockey_3_DM")


te.table.all <- as.data.frame(DGRP.pop.fig.alt %>% filter(ID %in% cent.te.to.keep) %>% dplyr::count(ID,freq.class,centromere,.drop = FALSE) %>% ungroup()) 
te.table.all <- as.data.frame(DGRP.pop.fig.alt %>% filter(ID %in% cent.te.to.keep, ref=="non-reference") %>% dplyr::count(ID,freq.class,centromere,.drop = FALSE) %>% ungroup()) 
te.table.sum <- as.data.frame(te.table.all %>% group_by(centromere) %>% summarise(sum.g = sum(n)))



pwd0 <- "/Users/lucashemmer/Documents/"
chrom.spans <- read.table(paste(pwd,"dmel_scaffold2_plus0310_chromatin.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)

chrom.spans <- as.data.frame(chrom.spans %>% mutate(centromere=ifelse(region=="centromere","centromere",
  ifelse(region!="centromere" & chromatin=="het","heterochromatin","euchromatin"))))

chrom.span.sum <- chrom.spans %>% filter(chromosome!="Y") %>% group_by(centromere) %>% summarise(sum(length))


## without Y chromosome
te.rare <- te.table.all %>% filter(freq.class=="Rare")

te.rare$exp.n <- NA
for (i in 1:nrow(te.rare)) {
  tmp.id <- te.rare$ID[i]
  tmp.sum <- sum(te.rare$n[te.rare$ID==tmp.id])
  if (te.rare$centromere[i] == "Centromere") {
    te.rare$exp.n[i] <- round(tmp.sum * 0.002068555)
  } else if (te.rare$centromere[i] == "Heterochromatin") {
    te.rare$exp.n[i] <- round(tmp.sum * 0.1660067)
  } else if (te.rare$centromere[i] == "Euchromatin") {
    te.rare$exp.n[i] <- round(tmp.sum * 0.8319248)
  }
}

te.rare$fisher.p <- NA
#te.rare$p.sig <- NA
for (i in seq(1,nrow(te.rare),3)) {
  tmp.het.test <- fisher.test(t(matrix(c(te.rare$n[i],te.rare$exp.n[i],
    te.rare$n[i+1],te.rare$exp.n[i+1]),ncol=2,nrow=2)))
  tmp.euc.test <- fisher.test(t(matrix(c(te.rare$n[i],te.rare$exp.n[i],
    te.rare$n[i+2],te.rare$exp.n[i+2]),ncol=2,nrow=2)))
  te.rare$fisher.p[i+1] <- tmp.het.test$p.value
  te.rare$fisher.p[i+2] <- tmp.euc.test$p.value
}

te.rare <- te.rare %>% mutate(p.sig = ifelse(fisher.p < 0.05, "sig", "not")) %>%
  mutate(p.adj = p.adjust(fisher.p)) %>% mutate(p.adj.sig = ifelse(p.adj < 0.05, "sig", "not"))


te.rare <- te.rare %>% mutate(n.mb = ifelse(centromere=="Centromere", n/0.292483, ifelse(centromere=="Heterochromatin",n/23.472482,n/117.629863)))


#### plot

te.rare.plot <- te.rare %>% filter(ID %in% c("Jockey_3_DM", "DOC","FW_DM","R1_DM","LINEJ1_DM"))
te.rare.plot$ID <- factor(te.rare.plot$ID, levels=c("Jockey_3_DM", "DOC","FW_DM","R1_DM","LINEJ1_DM"))
te.rare.plot$centromere <- factor(te.rare.plot$centromere, levels=c("Centromere","Heterochromatin","Euchromatin"))

#pdf("/Users/lucashemmer/Documents/Fig3_DGRP3_singleton_byRegion_Yremoved_presentable.pdf",width=6,height=3)
ggplot(te.rare.plot, aes(x=ID,y=n.mb,fill=centromere)) +
  geom_bar(position="dodge", stat="identity") +
  scale_x_discrete(labels=c("Jockey_3_DM" = "G2/Jockey-3","DOC" = "Doc","FW_DM" = "F","R1_DM" = "R1","LINEJ1_DM" = "Jockey-1")) +
  theme_bw() + theme(legend.position=c(0.7,0.7), legend.title = element_blank()) + 
  ylab("Insertions / Mb") + xlab("Genomic Region") + ylim(0,50) +
  scale_fill_viridis(discrete = T, option = "E", direction = 1) +
  ggtitle("Distribution of Singleton Insertions Per Region") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  # Jockey-3 cent vs euc
  geom_text(x = 1,  y = 46.3, label = "***", colour = "black") +
  geom_segment(x = 0.66, xend = 0.66, y = 45, yend = 46,colour = "black") +
  geom_segment(x = 1.33, xend = 1.33, y = 45, yend = 46,colour = "black") +
  geom_segment(x = 0.66, xend = 1.33, y = 46, yend = 46,colour = "black") +
  # Jockey-3 cent vs het
  geom_text(x = 0.833,  y = 43.3, label = "***", colour = "black") +
  geom_segment(x = 0.66, xend = 0.66, y = 42, yend = 43,colour = "black") +
  geom_segment(x = 1, xend = 1, y = 42, yend = 43,colour = "black") +
  geom_segment(x = 0.66, xend = 1, y = 43, yend = 43,colour = "black") +
  # R1 cent vs euc
  geom_text(x = 4,  y = 16.3, label = "*", colour = "black") +
  geom_segment(x = 3.66, xend = 3.66, y = 15, yend = 16,colour = "black") +
  geom_segment(x = 4.33, xend = 4.33, y = 15, yend = 16,colour = "black") +
  geom_segment(x = 3.66, xend = 4.33, y = 16, yend = 16,colour = "black") #+
#dev.off() 



## also add in all insertions

te.table.all <- as.data.frame(DGRP.pop.fig.alt %>% filter(ID %in% cent.te.to.keep) %>% dplyr::count(ID,freq.class,centromere,.drop = FALSE) %>% ungroup()) 
te.table.all <- te.table.all[order(te.table.all$centromere,te.table.all$ID,te.table.all$freq.class),]

## get percentage in each region
te.table.all$percent <- NA
for (i in 1:nrow(te.table.all)) {
  tmp.id <- te.table.all$ID[i]
  tmp.region <- te.table.all$centromere[i]
  tmp.sum <- sum(te.table.all$n[te.table.all$ID==tmp.id & te.table.all$centromere==tmp.region])
  te.table.all$percent[i] <- round((te.table.all$n[i] / tmp.sum) * 100)
}

## fisher exact test of differences between region

te.table.all <- te.table.all[order(te.table.all$ID,te.table.all$centromere,te.table.all$freq.class),]
te.table.all$percent[is.na(te.table.all$percent)] <- 0
te.table.all$fisher.p <- NA
for (i in seq(1,nrow(te.table.all),9)) {
  #tmp.het.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+3)],te.table.all$percent[(i+4):(i+7)]),ncol=2,nrow=4)))
  tmp.het.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+2)],te.table.all$percent[(i+3):(i+5)]),ncol=2,nrow=3)))
  tmp.euc.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+2)],te.table.all$percent[(i+6):(i+8)]),ncol=2,nrow=3)))
  te.table.all$fisher.p[i+3] <- tmp.het.test$p.value
  te.table.all$fisher.p[i+6] <- tmp.euc.test$p.value
}


## add significance and adjusted p values

te.table.all <- te.table.all %>% mutate(p.sig = ifelse(fisher.p < 0.05, "sig", "not")) %>%
  mutate(p.adj = p.adjust(fisher.p)) %>% mutate(p.adj.sig = ifelse(p.adj < 0.05, "sig", "not"))

## output fisher table 

write.table(te.table.all, paste(pwd0,"DGRP_distribution_fisher_results.txt",sep = ""), quote = F, row.names = F, sep="\t")



#pdf("/Users/lucashemmer/Documents/Figure_5_DGRP3_distribution_byRegion_twoTEs.pdf",width=6,height=8)
te.table.all %>% filter(ID %in% c("DOC","Jockey_3_DM")) %>%
  ggplot() + 
  geom_bar(aes(fill=freq.class,y=percent, x=centromere),position="dodge", stat="identity") + 
  scale_x_discrete(labels=c("centromere" = "cen", "heterochromatin" = "het","euchromatin" = "eu")) +
  facet_wrap(~ID, nrow = 2, ncol = 1) + 
  ggtitle("Percentage of Insertions by Region") + 
  theme_bw() + ylab("Percentage") + xlab("Genomic Region") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.title.align=0.5) + 
  guides(fill=guide_legend(title = "Frequency Class", nrow=1,byrow=TRUE)) + 
  scale_fill_viridis(discrete = T, option = "D") +
  geom_text(data = te.table.sum %>% filter(ID %in% c("DOC","Jockey_3_DM")), 
            aes(x = centromere, y = 90, label = paste0("N = ", sum.g))) +
  # DOC significant centromere vs euchromatin
  geom_text(data = subset(te.table.sum, ID=="DOC"), x = 2,  y = 79, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 1, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 3, xend = 3, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 3, y = 78, yend = 78,colour = "black") +
  # DOC significant centromere vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="DOC"), x = 1.5,  y = 69, label = "*", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 1, y = 66, yend = 68,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 2, xend = 2, y = 66, yend = 68,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 2, y = 68, yend = 68,colour = "black") +
  # DOC significant euchromatin vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="DOC"), x = 2.3,  y = 59, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 2, xend = 2, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 3, xend = 3, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 2, xend = 3, y = 58, yend = 58,colour = "black") +
  # Jockey-3 significant centromere vs euchromatin
  geom_text(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2,  y = 79, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 1, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 3, xend = 3, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 3, y = 78, yend = 78,colour = "black") +
  # Jockey-3 not significant centromere vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1.5,  y = 72, label = "ns", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 1, y = 66, yend = 68,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2, xend = 2, y = 66, yend = 68,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 2, y = 68, yend = 68,colour = "black") +
  # Jockey-3 significant euchromatin vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2.3,  y = 59, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2, xend = 2, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 3, xend = 3, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2, xend = 3, y = 58, yend = 58,colour = "black") 
#dev.off()




#pdf("/Users/lucashemmer/Documents/Figure_S3_DGRP3_distribution_byRegion_allTEs.pdf",width=12,height=6)
te.table.all %>%
  ggplot() + 
  geom_bar(aes(fill=freq.class,y=percent, x=centromere),position="dodge", stat="identity") + 
  scale_x_discrete(labels=c("centromere" = "cen", "heterochromatin" = "het","euchromatin" = "eu")) +
  facet_wrap(~ID, nrow = 4, ncol = 4) + 
  ggtitle("Percentage of Insertions by Region") + 
  theme_bw() + ylab("Percentage") + xlab("Genomic Region") + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.title.align=0.5) + 
  guides(fill=guide_legend(title = "Frequency Class", nrow=1,byrow=TRUE)) + 
  scale_fill_viridis(discrete = T, option = "D") + 
  # DOC significant centromere vs euchromatin
  geom_text(data = subset(te.table.sum, ID=="DOC"), x = 2,  y = 79, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 1, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 3, xend = 3, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 3, y = 78, yend = 78,colour = "black") +
  # DOC significant centromere vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="DOC"), x = 1.5,  y = 69, label = "*", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 1, y = 66, yend = 68,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 2, xend = 2, y = 66, yend = 68,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 1, xend = 2, y = 68, yend = 68,colour = "black") +
  # DOC significant euchromatin vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="DOC"), x = 2.3,  y = 59, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 2, xend = 2, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 3, xend = 3, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="DOC"), x = 2, xend = 3, y = 58, yend = 58,colour = "black") +
  # Jockey-3 significant centromere vs euchromatin
  geom_text(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2,  y = 79, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 1, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 3, xend = 3, y = 76, yend = 78,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 3, y = 78, yend = 78,colour = "black") +
  # Jockey-3 not significant centromere vs heterochromatin
  #geom_text(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1.5,  y = 72, label = "ns", colour = "black") +
  #geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 1, y = 66, yend = 68,colour = "black") +
  #geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2, xend = 2, y = 66, yend = 68,colour = "black") +
  #geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 1, xend = 2, y = 68, yend = 68,colour = "black") +
  # Jockey-3 significant euchromatin vs heterochromatin
  geom_text(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2.3,  y = 59, label = "***", colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2, xend = 2, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 3, xend = 3, y = 56, yend = 58,colour = "black") +
  geom_segment(data = subset(te.table.sum, ID=="Jockey_3_DM"), x = 2, xend = 3, y = 58, yend = 58,colour = "black") 
#dev.off()
















###################
## differences in the distribution of frequency classes of different elements
###################

## reorder frequency class summary table from earlier

te.table.all <- te.table.all[order(te.table.all$centromere,te.table.all$ID,te.table.all$freq.class),]

## get percentage in each region
te.table.all$percent <- NA
for (i in 1:nrow(te.table.all)) {
  tmp.id <- te.table.all$ID[i]
  tmp.region <- te.table.all$centromere[i]
  tmp.sum <- sum(te.table.all$n[te.table.all$ID==tmp.id & te.table.all$centromere==tmp.region])
  te.table.all$percent[i] <- round((te.table.all$n[i] / tmp.sum) * 100)
}

## fisher exact test of differences between region

te.table.all <- te.table.all[order(te.table.all$ID,te.table.all$centromere,te.table.all$freq.class),]
te.table.all$percent[is.na(te.table.all$percent)] <- 0
te.table.all$fisher.p <- NA
for (i in seq(1,nrow(te.table.all),12)) {
  tmp.het.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+3)],te.table.all$percent[(i+4):(i+7)]),ncol=2,nrow=4)))
  tmp.euc.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+3)],te.table.all$percent[(i+8):(i+11)]),ncol=2,nrow=4)))
  te.table.all$fisher.p[i+4] <- tmp.het.test$p.value
  te.table.all$fisher.p[i+8] <- tmp.euc.test$p.value
}

## add significance and adjusted p values

te.table.all <- te.table.all %>% mutate(p.sig = ifelse(fisher.p < 0.05, "sig", "not")) %>%
  mutate(p.adj = p.adjust(fisher.p)) %>% mutate(p.adj.sig = ifelse(p.adj < 0.05, "sig", "not"))

## output table 

#write.table(te.table.all, paste(pwdO,"DGRP3_distribution_fisher_results.txt",sep = ""), quote = F, row.names = F, sep="\t")

te.table.sum <- te.table.all %>% group_by(ID,centromere) %>% summarise(sum.g = sum(n))
te.table.sum <- left_join(te.table.sum, te.table.all %>% drop_na(p.adj) %>% select(ID, centromere, fisher.p, p.adj,p.adj.sig))


#pdf(paste(pwdO,"DGRP3_distribution_byRegion_allTEs.pdf",sep = ""),width=16,height=6)
ggplot(te.table.all) + 
  geom_bar(aes(fill=freq.class,y=percent, x=centromere),position="dodge", stat="identity") + 
  scale_x_discrete(labels=c("centromere" = "cen", "heterochromatin" = "het","euchromatin" = "eu")) +
  facet_wrap(~ID, nrow = 4, ncol = 4) + 
  ggtitle("Percentage of Insertions by Region") + 
  #stat_n_text() +
  theme_bw() + ylab("Percentage") + xlab("Genomic Region") + 
  #theme(legend.position=c(0.8,0.1), legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position=c(0.8,0.1), plot.title = element_text(hjust = 0.5), legend.title.align=0.5) + 
  guides(fill=guide_legend(title = "Frequency Class", nrow=1,byrow=TRUE)) + 
  scale_fill_viridis(discrete = T, option = "D") +
  geom_text(data = te.table.sum, aes(x = centromere, y = 90, label = paste0("N = ", sum.g)))
#dev.off()




























## without Y chromosome
te.singles2 <- te.table.all %>% filter(freq.class=="Singleton")

te.table.all <- as.data.frame(DGRP.pop.fig.alt %>% filter(ID %in% cent.te.to.keep, ref=="non-reference") %>% dplyr::count(ID,freq.class,centromere,.drop = FALSE) %>% ungroup()) 
te.singles2 <- te.table.all %>% filter(freq.class=="Singleton")


te.singles2$exp.n <- NA
for (i in 1:nrow(te.singles2)) {
  tmp.id <- te.singles2$ID[i]
  tmp.sum <- sum(te.singles2$n[te.singles2$ID==tmp.id])
  if (te.singles2$centromere[i] == "Centromere") {
    te.singles2$exp.n[i] <- round(tmp.sum * 0.002068555)
  } else if (te.singles2$centromere[i] == "Heterochromatin") {
    te.singles2$exp.n[i] <- round(tmp.sum * 0.1660067)
  } else if (te.singles2$centromere[i] == "Euchromatin") {
    te.singles2$exp.n[i] <- round(tmp.sum * 0.8319248)
  }
}

te.singles2$fisher.p <- NA
#te.singles2$p.sig <- NA
for (i in seq(1,nrow(te.singles2),3)) {
  tmp.het.test <- fisher.test(t(matrix(c(te.singles2$n[i],te.singles2$exp.n[i],
                                         te.singles2$n[i+1],te.singles2$exp.n[i+1]),ncol=2,nrow=2)))
  tmp.euc.test <- fisher.test(t(matrix(c(te.singles2$n[i],te.singles2$exp.n[i],
                                         te.singles2$n[i+2],te.singles2$exp.n[i+2]),ncol=2,nrow=2)))
  te.singles2$fisher.p[i+1] <- tmp.het.test$p.value
  te.singles2$fisher.p[i+2] <- tmp.euc.test$p.value
}

te.singles2 <- te.singles2 %>% mutate(p.sig = ifelse(fisher.p < 0.05, "sig", "not")) %>%
  mutate(p.adj = p.adjust(fisher.p)) %>% mutate(p.adj.sig = ifelse(p.adj < 0.05, "sig", "not"))
#write.table(te.singles2, file=paste(pwdO,"DGRP3_singleton_fisherResults_Yremoved.txt",sep = ""), quote = F, row.names = F, sep="\t")


te.singles2 <- te.singles2 %>% mutate(n.mb = ifelse(centromere=="Centromere", n/0.292483, ifelse(centromere=="Heterochromatin",n/23.472482,n/117.629863)))

#pdf(paste(pwdO,"DGRP3_singleton_byRegion_Yremoved.pdf",sep = ""),width=16,height=6)
ggplot(te.singles2, aes(x=centromere,y=n.mb, fill=centromere)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels=c("Centromere" = "cen", "Heterochromatin" = "het","Euchromatin" = "eu")) +
  #scale_y_continuous(breaks=c(0, 2, 4, 8, 10))+
  facet_wrap(~ID, nrow = 2, ncol = 7) +
  ggtitle("Singleton Insertions Per Region")+
  theme_bw() + ylab("Insertions / Mb") + xlab("Genomic Region") + 
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
  scale_fill_viridis(discrete = T, option = "E") 
#dev.off()


## reorder frequency class summary table from earlier

te.table.all <- te.table.all[order(te.table.all$centromere,te.table.all$ID,te.table.all$freq.class),]

## get percentage in each region
te.table.all$percent <- NA
for (i in 1:nrow(te.table.all)) {
  tmp.id <- te.table.all$ID[i]
  tmp.region <- te.table.all$centromere[i]
  tmp.sum <- sum(te.table.all$n[te.table.all$ID==tmp.id & te.table.all$centromere==tmp.region])
  te.table.all$percent[i] <- round((te.table.all$n[i] / tmp.sum) * 100)
}

## fisher exact test of differences between region

te.table.all <- te.table.all[order(te.table.all$ID,te.table.all$centromere,te.table.all$freq.class),]
te.table.all$percent[is.na(te.table.all$percent)] <- 0
te.table.all$fisher.p <- NA
for (i in seq(1,nrow(te.table.all),12)) {
  tmp.het.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+3)],te.table.all$percent[(i+4):(i+7)]),ncol=2,nrow=4)))
  tmp.euc.test <- fisher.test(t(matrix(c(te.table.all$percent[i:(i+3)],te.table.all$percent[(i+8):(i+11)]),ncol=2,nrow=4)))
  te.table.all$fisher.p[i+4] <- tmp.het.test$p.value
  te.table.all$fisher.p[i+8] <- tmp.euc.test$p.value
}

## add significance and adjusted p values

te.table.all <- te.table.all %>% mutate(p.sig = ifelse(fisher.p < 0.05, "sig", "not")) %>%
  mutate(p.adj = p.adjust(fisher.p)) %>% mutate(p.adj.sig = ifelse(p.adj < 0.05, "sig", "not"))

## output table 

#write.table(te.table.all, paste(pwdO,"DGRP3_distribution_fisher_results.txt",sep = ""), quote = F, row.names = F, sep="\t")

te.table.sum <- te.table.all %>% group_by(ID,centromere) %>% summarise(sum.g = sum(n))
te.table.sum <- left_join(te.table.sum, te.table.all %>% drop_na(p.adj) %>% select(ID, centromere, fisher.p, p.adj,p.adj.sig))


#pdf(paste(pwdO,"DGRP3_distribution_byRegion_allTEs.pdf",sep = ""),width=16,height=6)
ggplot(te.table.all) + 
  geom_bar(aes(fill=freq.class,y=percent, x=centromere),position="dodge", stat="identity") + 
  scale_x_discrete(labels=c("centromere" = "cen", "heterochromatin" = "het","euchromatin" = "eu")) +
  facet_wrap(~ID, nrow = 4, ncol = 4) + 
  ggtitle("Percentage of Insertions by Region") + 
  #stat_n_text() +
  theme_bw() + ylab("Percentage") + xlab("Genomic Region") + 
  #theme(legend.position=c(0.8,0.1), legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), legend.title.align=0.5) + 
  guides(fill=guide_legend(title = "Frequency Class", nrow=1,byrow=TRUE)) + 
  scale_fill_viridis(discrete = T, option = "D") +
  geom_text(data = te.table.sum, aes(x = centromere, y = 90, label = paste0("N = ", sum.g)))
#dev.off()






