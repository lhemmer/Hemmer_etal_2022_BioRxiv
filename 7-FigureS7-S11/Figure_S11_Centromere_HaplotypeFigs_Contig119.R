#### ####  CONTIG119 AKA 4TH CHROMOSOME CENTROMERE

############## starting over again from the top



#### load libraries

library(ape)
library(vcfR)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)

#### #### 
#### dendrogram based on snps
#### #### 

library(gdsfmt)
library(SNPRelate)
snpgdsVCF2GDS("All_filtered_snps.variant.select.Contig119.nosat.vcf","Contig119.nosat.gds",method ="biallelic.only")

genofile1<-snpgdsOpen("Contig119.nosat.gds")
set.seed(100)
ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile1,num.thread=2, autosome.only=FALSE))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram,main="Dendrogram based on IBS")


samples.ordered <- rv$sample.id[rv$samp.order]


#### figure from vcf


library(ape)
library(vcfR)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)

#### load data

vcf <- read.vcfR("All_filtered_snps.variant.select.Contig119.2.vcf", verbose = FALSE)
dna <- ape::read.dna("Contig119.fasta", format = "fasta")
gff <- read.table("Contig119.gff", sep="\t", quote="")


# Create a chromR object.
chrom <- create.chromR(name="Contig119", vcf=vcf, seq=dna, ann=gff, verbose=TRUE)

x <- vcfR2genlight(vcf)

## working with genotypes
gt <- extract.gt(vcf, element = "GT")

## melt the dataframe
gt.melt <- melt(gt)
colnames(gt.melt) <- c("contig_pos","sample","genotype")

## recoding the genotypes for plotting

## first is changing het genotypes to R/A

gt.proc <- gt.melt %>% separate(contig_pos, into=c("contig","pos"),sep="_",remove=F) %>%
  mutate(geno.recode = if_else(is.na(genotype), "NA", 
                               if_else(genotype %in% c("0/0","0|0"), "R", 
                                       if_else(genotype %in% c("0/1", "0|1", "0/2", "0|2", "0/3", "0|3", "0/4", "0|4"), "R/A", "A"))))

## second is changing het genotypes to NA

gt.proc2 <- gt.melt %>% separate(contig_pos, into=c("contig","pos"),sep="_",remove=F) %>%
  mutate(geno.recode = if_else(is.na(genotype), "NA", 
                               if_else(genotype %in% c("0/0","0|0"), "R", 
                                       if_else(genotype %in% c("0/1", "0|1", "0/2", "0|2", "0/3", "0|3", "0/4", "0|4"), "NA", "A"))))

gt.proc <- gt.proc %>% mutate(pos = as.numeric(pos), sample = as.character(sample))
gt.proc2 <- gt.proc2 %>% mutate(pos = as.numeric(pos), sample = as.character(sample))


#### loading gff file

# process
gff1 <- gff %>% separate(V9, into=c("Target", "repeats", "te.start","te.end"), sep="\\s", remove = T) %>% 
  mutate(repeats = str_remove_all(repeats, "Motif:")) %>% mutate(repeats = str_remove_all(repeats, '"')) %>%
  select(-c( V2, V3, V7, V8, Target)) %>% rename(contig="V1", contig.start = "V4", contig.end = "V5", div = "V6") %>%
  mutate(te.start = as.numeric(te.start), te.end = as.numeric(te.end))

## isolate large satellite expanses and jockey-3
gff.sat <- gff1 %>% filter(str_detect(repeats, "\\)n"), te.end - te.start > 500)
gff.g2 <- gff1 %>% filter(str_detect(repeats, "Jockey-3_Dmel_08212020"))


gt.proc$geno.recode <- factor(gt.proc$geno.recode, levels = c("R", "R/A", "A", "NA"))

gt.proc2$geno.recode <- factor(gt.proc2$geno.recode, levels = c("R", "A", "NA"))


### get matches for satellite and G2/Jockey-3 dataframes

## subset
gt.samp <- gt.proc %>% filter(sample == sort(unique(gt.proc$sample[1])))
gt.samp$contig_pos <- as.character(gt.samp$contig_pos)

## getting positions for satellite dna
gff.sat$contig_pos1 <- NA
gff.sat$contig_pos2 <- NA
for (i in 1:nrow(gff.sat)) {
  tmp.start <- which(abs(gt.samp$pos-gff.sat$contig.start[i])==min(abs(gt.samp$pos - gff.sat$contig.start[i])))
  tmp.end <- which(abs(gt.samp$pos-gff.sat$contig.end[i])==min(abs(gt.samp$pos - gff.sat$contig.end[i])))
  gff.sat$contig_pos1[i] <- gt.samp$contig_pos[tmp.start]
  gff.sat$contig_pos2[i] <- gt.samp$contig_pos[tmp.end]
}
gff.sat$ymin <- -2
gff.sat$ymax <- -1

## getting information for jockey-3
gff.g2$contig_pos1 <- NA
gff.g2$contig_pos2 <- NA
for (i in 1:nrow(gff.g2)) {
  tmp.start <- which(abs(gt.samp$pos-gff.g2$contig.start[i])==min(abs(gt.samp$pos - gff.g2$contig.start[i])))
  tmp.end <- which(abs(gt.samp$pos-gff.g2$contig.end[i])==min(abs(gt.samp$pos - gff.g2$contig.end[i])))
  gff.g2$contig_pos1[i] <- gt.samp$contig_pos[tmp.start]
  gff.g2$contig_pos2[i] <- gt.samp$contig_pos[tmp.end]
}
gff.g2$ymin <- -1
gff.g2$ymax <- -0


## eliminate satellite DNA, order based on dendrograme

gt.proc3 <- gt.proc %>% mutate(pos = as.numeric(pos), sample = as.character(sample)) %>% 
  arrange(sample,pos) %>% filter(pos >= 29253, pos <=69147)

gt.proc3$sample <- factor(gt.proc3$sample, levels=samples.ordered)


gt.proc4 <- gt.proc2 %>% mutate(pos = as.numeric(pos), sample = as.character(sample)) %>% 
  arrange(sample,pos) %>% filter(pos >= 29253, pos <=69147)

gt.proc4$sample <- factor(gt.proc4$sample, levels=samples.ordered)


#####  plot

## R/A replaced with NA 

#pdf("SNP_poly_Contig119_ordered_11.pdf",width=10,height=8)
gt.proc4 %>% ggplot() +
  geom_raster(aes(x = contig_pos, y = sample, fill = geno.recode)) +
  scale_y_discrete(expand = expansion(add=c(1, 0))) +
  #geom_rect(aes(ymin=-2, xmin=1, ymax=-1, xmax=1000), fill="steelblue") +
  #geom_rect(data = gff.sat, aes(ymin=ymin, xmin=contig_pos1, ymax=ymax, xmax=contig_pos2), fill="steelblue") +
  geom_rect(data = gff.g2, aes(ymin=ymin, xmin=contig_pos1, ymax=ymax, xmax=contig_pos2), fill="firebrick") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank()) +
  #scale_fill_manual(values=c("#3C7DC4", "#FF8F00","#7C7B78FF"), name = "Genotype")
  #scale_fill_manual(values=c("#243743", "#28B78D","#7C7B78FF"), name = "Genotype")
  scale_fill_manual(values=c("#2F4E6F", "#E15119","#7C7B78FF"), name = "Genotype")
#dev.off()



#### #### #### 
#### adding DGRP

DGRP.all <- read.table("File_S2_McClintock_output_DGRP.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DGRP.all$ID[DGRP.all$ID=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"

pops <- as.numeric(as.vector(unique(gt.proc$sample)))

DGRP.select <- DGRP.all %>% filter(contig %in% c("Contig79", "Contig119", "tig00057289", "3R_5"), population %in% pops, ID!="Gypsy6_I_Dmoj", ID!="Gypsy1_I_Dmoj")

gt.te <- DGRP.select %>% select(identifier, population, ref)
colnames(gt.te) <- c("identifier", "sample", "ref")
gt.te <- gt.te %>% separate(identifier, into=c("genotype","contig", "pos"), sep="\\.")

gt.te <- gt.te %>% unite(col=contig_pos, contig, pos, sep="_", remove=F) %>% 
  mutate(geno.recode=if_else(genotype=="Jockey_3_Dmel_08212020", "Jockey_3_DM", "Other"))

gt.te <- gt.te %>% select(c(contig_pos, contig, pos, sample, genotype, geno.recode, ref))

gt.te.4 <- gt.te %>% filter(contig=="Contig119")



#### only jockey-3

te.combos.j1 <- gt.te.4 %>% filter(genotype == "Jockey_3_Dmel_08212020") %>% tidyr::expand(contig_pos, sample) #, geno.recode)
te.combos.j2 <- te.combos.j1 %>% left_join(gt.te.4 %>% select(c(contig_pos, sample, genotype, geno.recode ))) 
te.combos.j2 <- te.combos.j2 %>% separate(contig_pos, into=c("contig", "pos"), sep="_", remove=F)
te.combos.j2$geno.recode[is.na(te.combos.j2$geno.recode)] <- "No insertion"
te.combos.j2 <- te.combos.j2 %>% select(c(contig_pos, contig, pos, sample, genotype, geno.recode))

te.combos.j2$genotype[te.combos.j2$genotype=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"

te.combos.j2 <- te.combos.j2 %>% mutate(pos = as.numeric(pos), sample = as.character(sample)) %>% 
  arrange(sample,pos) %>% filter(pos >= 29253, pos <=69147)

te.combos.j3 <- te.combos.j2 %>% distinct() %>% filter()
te.combos.cast.j3 <- dcast(te.combos.j3, sample ~ contig_pos , fun.aggregate = NULL)

te.combos.cast.j3[te.combos.cast.j3 == "No insertion"] <- 0
te.combos.cast.j3[te.combos.cast.j3 == "Jockey_3_DM"] <- 1

te.combos.cast.j3[, c(2:ncol(te.combos.cast.j3))] <- sapply(te.combos.cast.j3[, c(2:ncol(te.combos.cast.j3))], as.numeric)

ncol(te.combos.cast.j3)



te.combos.cast.jh <- as.data.frame(scale(te.combos.cast.j3[,2:ncol(te.combos.cast.j3)]))
rownames(te.combos.cast.jh) <- te.combos.cast.j3[,1]
te.combos.dist.jc <- dist(te.combos.cast.jh,method = 'euclidean')

hclust_j3 <- hclust(te.combos.dist.jc, method = 'average')
plot(hclust_j3)


te.combos.j3$geno.recode <- factor(te.combos.j3$geno.recode, levels = c("Jockey_3_DM", "No insertion"))
te.combos.j3$sample <- factor(te.combos.j3$sample, levels=samples.ordered)


#pdf("/Users/lucashemmer/Documents/TE_poly_allTE_Contig119_ordered_1.pdf",width=10,height=8)
te.combos.j3 %>% ggplot() +
  geom_raster(aes(x = contig_pos, y = sample, fill = geno.recode)) +
  #scale_y_discrete(expand = expansion(add=c(2, 0))) +
  scale_y_discrete(position = "right") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "left") +
  scale_fill_manual(values=c("#E15119", "#2F4E6F"), name = "Genotype")
#dev.off()



### snps vs all TEs
den <- as.dendrogram(ibs.hc$dendrogram)
denTE <- as.dendrogram(hclust_allTE)
denJ3 <- as.dendrogram(hclust_j3)


### j3 vs snps

d2 <- dendextend::dendlist(
  den %>% 
    dendextend::set("branches_lty", 1), #%>%
  denJ3 %>% 
    dendextend::set("branches_lty", 1) #%>%
)

dendextend::tanglegram(d2, 
                       common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
                       margin_inner=5, lwd=2
)


#pdf("/Users/lucashemmer/Documents/tanglegrame_Contig119_SNPs_J3_1.pdf",width=10,height=8)
dendextend::dendlist(
  den %>% 
    dendextend::set("branches_lty", 1), #%>%
  denJ3 %>% 
    dendextend::set("branches_lty", 1) ) %>% 
  dendextend::untangle(method = "step1side") %>%  
  dendextend::tanglegram(common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
                         margin_inner=2, lwd=2)
#dev.off()



### put it all together now

pdf("Contig119_SNPraster_p1.pdf",width=8,height=6)
p1 <- gt.proc4 %>% ggplot() +
  geom_raster(aes(x = contig_pos, y = sample, fill = geno.recode)) +
  scale_y_discrete(expand = expansion(add=c(1, 0))) +
  #geom_rect(aes(ymin=-2, xmin=1, ymax=-1, xmax=1000), fill="steelblue") +
  #geom_rect(data = gff.sat, aes(ymin=ymin, xmin=contig_pos1, ymax=ymax, xmax=contig_pos2), fill="steelblue") +
  geom_rect(data = gff.g2, aes(ymin=ymin, xmin=contig_pos1, ymax=ymax, xmax=contig_pos2), fill="firebrick") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank()) +
  #scale_fill_manual(values=c("#3C7DC4", "#FF8F00","#7C7B78FF"), name = "Genotype")
  #scale_fill_manual(values=c("#243743", "#28B78D","#7C7B78FF"), name = "Genotype")
  scale_fill_manual(values=c("#2F4E6F", "#E15119","#7C7B78FF"), name = "Genotype")
p1
dev.off()


pdf("Contig119_J3raster_p2.pdf",width=8,height=6)
p2 <- te.combos.j3 %>% ggplot() +
  geom_raster(aes(x = contig_pos, y = sample, fill = geno.recode)) +
  #scale_y_discrete(expand = expansion(add=c(2, 0))) +
  scale_y_discrete(position = "right", limits=samples.ordered, expand = expansion(add=c(2, 0))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "left") +
  scale_fill_manual(values=c("#D95F02", "#1B9E77"), name = "Genotype")
p2
dev.off()


pdf("Contig119_SNPs_J3_tanglegrame_p3.pdf",width=14,height=6)
dendextend::dendlist(
  den %>% 
    dendextend::set("branches_lty", 1), #%>%
  denJ3 %>% 
    dendextend::set("branches_lty", 1) ) %>% 
  dendextend::untangle(method = "step1side") %>%  
  dendextend::tanglegram(common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
                         margin_inner=3, lwd=2)
dev.off()





