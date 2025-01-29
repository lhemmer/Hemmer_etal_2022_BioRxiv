######## 
#### Figure showing the distribution of G2/Jockey-3 in the DGRP
######## 

#### Load libraries

library(DescTools)
library(viridis)
library(ggplot2)
library(tidyr)
library(dplyr)
options(scipen=999)


#### Load data

pwd <- "/Users/lucashemmer/Documents/mel_centromere/Paper/"
DGRP.all <- read.table(paste(pwd,"File_S2_McClintock_output_DGRP.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DGRP.all$ID[DGRP.all$ID=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"


## using jockey-3 data only

DGRP.g2 <- DGRP.all %>% filter(ID=="Jockey_3_DM") %>% distinct(identifier, .keep_all = TRUE)


#### changing positions for illustration

## documenting reverse contigs, which contigs belong to which chromosome

reverse_contigs=c('2L_2','2R_contig_19','3L_4','3R_5','3_scaffold1','Contig5','Y_scaffold7','Y_Contig10','Y_Contig104')
chrom.x <- c('X_1','Contig5','X_2','X_3','X_7','X_9','Contig135','X_10','X3X4_6_D1712_2','Contig95','Contig79')
chrom.2 <- c('2L_1','2L_2','tig00057289','Contig142','Rsp_2','2R_16','2R_18','2R_19','2R_20','2R_21')
chrom.3 <- c('3L_1','3L_3','3L_4','3L_5','3L_8','3L_10','3R_1','3R_2','3R_3','3R_4','3R_5','tig00022795','id=67822_0','3R_6','3R_8','3R_9','3R_10','Contig11','Contig145','3_scaffold1','3_scaffold2','3R_28')
chrom.4 <- c('Contig119','4_2')
#chrom.Y <- c('Y_scaffold6','Y_scaffold7','Y_scaffold4','Y_Contig140','Y_Contig143','Y_Contig10','Y_Contig104','Y_Contig6','Y_Contig2','Y_Contig26','Y_scaffold5','Y_scaffold3')

chroms <- list(chrom.x,chrom.2,chrom.3,chrom.4)

## separate out jockey-3 calls by chromosome

te.by.chrom <- list(chromX=NA,chrom2=NA,chrom3=NA,chrom4=NA)
for (i in 1:length(chroms)) {
  te.by.chrom[[i]] <- DGRP.g2[DGRP.g2$contig %in% chroms[[i]],]
  te.by.chrom[[i]]$start.new <- NA
  te.by.chrom[[i]]$end.new <- NA
}

## load sizes of each contig
pwd <- "/Users/lucashemmer/Documents/mel_centromere/reference/"
sizes <- read.csv(paste(pwd,"dmel_scaffold2_plus0310_sizes.csv",sep=""), header=FALSE)
colnames(sizes) <- c("contig","len")


#### run loop to change position based on contig, forward or reversed, add 100,000 bp between contig segments

for (k in 1:length(chroms)) { 
  for (i in 1:length(te.by.chrom[[k]]$contig)) {
    ## intial position for each chromosome
    pos <- 0
    for (j in 1:length(chroms[[k]])) {
      if (te.by.chrom[[k]]$contig[i]==chroms[[k]][j]) {
        # for reverse contigs
        if (chroms[[k]][j] %in% reverse_contigs) {
          te.by.chrom[[k]]$start.new[i] <- sizes$len[sizes$contig==chroms[[k]][j]] - te.by.chrom[[k]]$start[i] + pos
          te.by.chrom[[k]]$end.new[i] <- sizes$len[sizes$contig==chroms[[k]][j]] - te.by.chrom[[k]]$end[i] + pos
        } else { # for other contigs
          te.by.chrom[[k]]$start.new[i] <- te.by.chrom[[k]]$start[i] + pos
          te.by.chrom[[k]]$end.new[i] <- te.by.chrom[[k]]$end[i] + pos
        }
      }
      ## add 100000 bp for breaks between contigs
      pos <- pos+sizes$len[sizes$contig==chroms[[k]][j]]+100000
    }
  }
}

## transform list into dataframe again, reorder

te.by.chrom <- do.call(rbind, te.by.chrom)
rownames(te.by.chrom) <- 1:nrow(te.by.chrom)
te.by.chrom <- te.by.chrom[order(te.by.chrom$chromosome,te.by.chrom$start.new),]


#### Plotting different chromosomes

## X chromosome

x.arms <- te.by.chrom %>% filter(chromosome == "X") %>%
  ggplot(aes(x=(start.new/1000000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start.new/1000000,xmax=start.new/1000000+(0.015*25.5),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3"))) +
  theme(panel.border = element_blank(),axis.line = element_line(color='black'),panel.background=element_rect(fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(5.3,15.3,25.3),labels=c(-20,-10,0),limits=c(0,26)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Ch X position (Mb)") +
  ylab("") 


## second chromosome

two.arms <- te.by.chrom %>% filter(chromosome == "2", contig!="2R_19") %>%
  ggplot(aes(x=(start.new/1000000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start.new/1000000,xmax=start.new/1000000+(0.015*25.5),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3"))) +
  theme(panel.border = element_blank(),axis.line = element_line(color='black'),panel.background=element_rect(fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(3.75,13.75,23.75,33.75,43.75),labels=c(-20,-10,0,10,20),limits=c(0,50.5)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Ch 2 position (Mb)") +
  ylab("") 


## third chromosome

three.arms <- te.by.chrom %>% filter(chromosome == "3") %>%
  ggplot(aes(x=(start.new/1000000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start.new/1000000,xmax=start.new/1000000+(0.015*25.5),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3"))) +
  theme(panel.border = element_blank(),axis.line = element_line(color='black'),panel.background=element_rect(fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(9.5,19.5,29.5,39.5,49.5,59.5),labels=c(-20,-10,0,10,20,30),limits=c(0,65)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Ch 3 position (Mb)") +
  ylab("")


## fourth chromosome

four.arms <-  te.by.chrom %>% filter(chromosome == "4") %>%
  ggplot(aes(x=(start.new/1000000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start.new/1000000,xmax=start.new/1000000+(0.025*1.6),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3"))) +
  theme(panel.border = element_blank(),axis.line = element_line(color='black'),panel.background=element_rect(fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(-0.05,0.95),labels=c(0,1),limits=c(0,1.6)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Ch 4 position (Mb)") +
  ylab("") 


## x centromere

x.cent <-  te.by.chrom %>% filter(contig == "Contig79") %>% 
  ggplot(aes(x=(start/10000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start/10000,xmax=start/10000+(0.01*(sizes$len[sizes$contig=="Contig79"])/10000),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3"))) +
  theme(panel.background=element_rect(colour = "black", fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12))+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(1.75,3.5,5.25),labels=c(-20,0,20),limits=c(0,(sizes$len[sizes$contig=="Contig79"])/10000)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("CenX position (kb)") +
  ylab("") 


## 2nd cent

two.cent <- te.by.chrom %>% filter(contig == "tig00057289") %>% 
  ggplot(aes(x=(start/10000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start/10000,xmax=start/10000+(0.02*(sizes$len[sizes$contig=="tig00057289"])/10000),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("maroon", "maroon"),name = expression(italic("G2/Jockey-3")))+
  theme(panel.background=element_rect(colour = "black", fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12))+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c((sizes$len[sizes$contig=="tig00057289"])/40000,(sizes$len[sizes$contig=="tig00057289"])/20000,(sizes$len[sizes$contig=="tig00057289"])*3/40000),labels=c(-10,0,10),limits=c(0,(sizes$len[sizes$contig=="tig00057289"])/10000)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Cen2 position (kb)") +
  ylab("")


## 3rd centromere

three.cent <- te.by.chrom %>% filter(contig == "3R_5") %>% 
  ggplot(aes(x=(start/10000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start/10000,xmax=start/10000+(0.015*(sizes$len[sizes$contig=="3R_5"])/10000),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3")))+
  theme(panel.background=element_rect(colour = "black", fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12))+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(1.19135,3.19135,(sizes$len[sizes$contig=="3R_5"])/20000,7.19135,9.19135),labels=c(-40,-20,0,20,40),limits=c(0,(sizes$len[sizes$contig=="3R_5"])/10000)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Cen3 position (kb)") +
  ylab("")


## 4th centromere

four.cent <- te.by.chrom %>% filter(contig == "Contig119") %>% 
  ggplot(aes(x=(start/10000), y=pop.count)) +
  geom_rect(mapping=aes(xmin=start/10000,xmax=start/10000+(0.015*(sizes$len[sizes$contig=="Contig119"])/10000),ymin=1,ymax=pop.count+1,fill=ref),alpha=0.7) +
  scale_fill_manual(values=c("royalblue1", "maroon"),name = expression(italic("G2/Jockey-3"))) +
  theme(panel.background=element_rect(colour = "black", fill = "white")) +
  theme(panel.grid.minor = element_line(linetype = 'blank'),panel.grid.major = element_line(linetype = 'blank'),axis.text=element_text(size=12))+
  theme(legend.position="none")+
  scale_x_continuous(breaks=c(0.6957,2.6957,(sizes$len[sizes$contig=="Contig119"])/20000,6.6957,8.6957),labels=c(-40,-20,0,20,40),limits=c(0,(sizes$len[sizes$contig=="Contig119"])/10000)) +
  scale_y_continuous(breaks=c(1,27,54,81,108),labels=c(0,0.2,0.4,0.6,0.8),limits=c(-1.5,135)) +
  xlab("Cen4 position (kb)") +
  ylab("") 


#### export plots

## empty plot for spacing
df <- data.frame()
pg <- ggplot(df) + geom_blank(mapping=NULL) + theme_void()

pdf("Figure_3_DGRP3_insertion_freq_chromosome_all_mod.pdf",width=10,height=12)
x <- ggarrange(two.arms,ggarrange(two.cent,pg,three.cent,ncol=3,widths=c(2,1,4)),
          three.arms,ggarrange(four.arms,x.arms,pg,ncol=3,widths=c(2,3,1),labels=c("C","D","")),
          ggarrange(four.cent,pg,x.cent,ncol=3,widths=c(3,1,3)),nrow=5,heights = c(1,0.7,1,1,0.7),
          labels=c("A","","B","","",""))

annotate_figure(x, left = text_grob("Frequency in DGRP", rot = 90))

dev.off()


