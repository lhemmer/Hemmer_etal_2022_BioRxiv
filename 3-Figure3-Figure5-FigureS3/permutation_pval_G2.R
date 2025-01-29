######## 
#### Analysis of Distribution of Transposable Elements, Permutation test
######## 

#### Load libraries

library(DescTools)
library(viridis)
library(ggplot2)
library(tidyr)
library(dplyr)
library(regioneR)
library(ape)
library(stringr)
options(scipen=999)

#### Load data

DGRP.all <- read.table("File_S2_McClintock_output_DGRP.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DGRP.all$ID[DGRP.all$ID=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"


## higher coverage sample calls
DGRP.cent.cov <- read.table("Table_S3_DGRP_total_cent_cov_500bpWin.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
cent.hcov.pop <- DGRP.cent.cov %>% filter(centromere=="cent" & med_cov > 20 & population!=340 & population!=373) %>% pull(population) 


## long reads
jock.info <- read.table("File_S5_Jockey-3_Phylogeny_Info.txt",header=T,sep="\t",stringsAsFactors = F, fill=T)


## spans
chrom.spans <- read.table("dmel_scaffold2_plus0310_chromatin.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)


## chromosome sizes
chrom.bed <- read.table("dmel_scaffold2_plus0310_sizes.bed", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)


########################
#### Long read permutation test
########################


#### adding info 

## presence or absense of complete ORF (not necessarily coding)

jock.info <- jock.info %>% mutate(orf1 = if_else(te.start <= 186 & te.end >= 1469, "Complete", "Incomplete")) %>% 
  mutate(orf2 = if_else(te.start <= 1457 & te.end >= 4156, "Complete", "Incomplete")) %>% 
  mutate(orfs.present = if_else(orf1 == "Complete" | orf2 == "Complete", "ORF1 or ORF2", "Incomplete"))


#### filter out de novo

jock.bed.denovo <-  jock.info %>% filter(denovo == "yes", keep.delete=="keep", ref.contig!="unknown") %>% 
  select(ref.contig, ref.start, ref.chromatin, ref.region) %>% 
  mutate(ref.end = ref.start + 1) %>% 
  rename(chr = ref.contig) %>% 
  select(chr, ref.start, ref.end, ref.chromatin, ref.region) 


## divide into regions
jock.cent.denovo <- jock.bed.denovo %>% filter(ref.region == "centromere")
jock.euc.denovo <- jock.bed.denovo %>% filter(ref.chromatin == "eu")
jock.het.denovo <- jock.bed.denovo %>% filter(ref.chromatin == "het", ref.region != "centromere")


#### split up genome into 100Kb tiles

cent.spans <- chrom.spans %>% filter(region=="centromere") %>% select(contig_name, start.range, end.range) %>% mutate(chromatin = "Centromere")
euc.spans <- chrom.spans %>% filter(chromatin=="eu") %>% select(contig_name, start.range, end.range)
het.spans <- chrom.spans %>% filter(chromatin=="het", region!="centromere") %>% select(contig_name, start.range, end.range)


## convert to genomic ranges
cent.range <- toGRanges(cent.spans)
euc.range <- toGRanges(euc.spans)
het.range <- toGRanges(het.spans)


## tile up het and euchromatin

euc.tiles <- tile(euc.range, width = 90000)
het.tiles <- tile(het.range, width = 90000)

tmp.euc.tiles <- as.data.frame(euc.tiles)
tmp.het.tiles <- as.data.frame(het.tiles)

cent.wind <- cent.spans %>% rename(chr = contig_name, start = start.range, end = end.range)
euc.wind <- tmp.euc.tiles %>% select(seqnames, start, end) %>% rename(chr = seqnames) %>% mutate(chromatin = "Euchromatin")
het.wind <- tmp.het.tiles %>% select(seqnames, start, end) %>% rename(chr = seqnames) %>% mutate(chromatin = "Heterochromatin")


#### overlaps

## centromere

x <- toGRanges(cent.wind)
y <- toGRanges(jock.cent.denovo)

overlaps <- countOverlaps(x, y)
cent.wind$G2.count <- as.vector(overlaps)

## euchromatin

x <- toGRanges(euc.wind)
y <- toGRanges(jock.euc.denovo)

overlaps <- countOverlaps(x, y)
euc.wind$G2.count <- as.vector(overlaps)

## heterochromatin

x <- toGRanges(het.wind)
y <- toGRanges(jock.het.denovo)

overlaps <- countOverlaps(x, y)
het.wind$G2.count <- as.vector(overlaps)


#### combine altogether

genome.wind <- rbind.data.frame(cent.wind, euc.wind, het.wind)

genome.wind <- genome.wind %>% mutate(Ychrom = if_else(str_detect(chr, "Y"), "Y", "other"))


## jockey-3 counts 


cent.count <- genome.wind %>% filter(chromatin == "Centromere") %>% summarise(total = sum(G2.count)) %>% .$total
euc.count <- genome.wind %>% filter(chromatin == "Euchromatin") %>% summarise(total = sum(G2.count)) %>% .$total
het.count <- genome.wind %>% filter(chromatin == "Heterochromatin") %>% summarise(total = sum(G2.count)) %>% .$total



#### scramble / permutate


cent.list <- vector()
euc.list <- vector()
het.list <- vector()

set.seed(2021)
for (i in 1:10000) {
  tmp.df <- genome.wind
  tmp.df <- tmp.df %>% mutate(chromatin = sample(chromatin))
  
  cent.list[i] <- tmp.df %>% filter(chromatin == "Centromere") %>% summarise(total = sum(G2.count)) %>% .$total
  euc.list[i] <- tmp.df %>% filter(chromatin == "Euchromatin") %>% summarise(total = sum(G2.count)) %>% .$total
  het.list[i] <- tmp.df %>% filter(chromatin == "Heterochromatin") %>% summarise(total = sum(G2.count)) %>% .$total
}



permutation.test <- function(df, n, cent.count, euc.count, het.count){
  # set up variables
  tmp.df <- df
  cen.dist=c()
  euc.dist=c()
  het.dist=c()
  cen.result=0
  euc.result=0
  het.result=0
  # scramble loop
  for(i in 1:n){
    tmp.df <- tmp.df %>% mutate(chromatin = sample(chromatin))
    cen.dist[i] <- tmp.df %>% filter(chromatin == "Centromere") %>% summarise(total = sum(G2.count)) %>% .$total
    euc.dist[i] <- tmp.df %>% filter(chromatin == "Euchromatin") %>% summarise(total = sum(G2.count)) %>% .$total
    het.dist[i] <- tmp.df %>% filter(chromatin == "Heterochromatin") %>% summarise(total = sum(G2.count)) %>% .$total
  }
  # p vals
  cen.result <- sum(abs(cen.dist) >= abs(cent.count))/(n)
  euc.result <- sum(abs(euc.dist) >= abs(euc.count))/(n)
  het.result <- sum(abs(het.dist) >= abs(het.count))/(n)
  
  return(list(cen.result, cen.dist, euc.result, euc.dist, het.result, het.dist))
}



permutation.test2 <- function(df, n, cent.count, het.count){
  # set up variables
  tmp.df <- df
  cen.dist=c()
  euc.dist=c()
  het.dist=c()
  cen.result=0
  #euc.result=0
  het.result=0
  # scramble loop
  for(i in 1:n){
    tmp.df <- tmp.df %>% mutate(chromatin = sample(chromatin))
    cen.dist[i] <- tmp.df %>% filter(chromatin == "Centromere") %>% summarise(total = sum(G2.count)) %>% .$total
    #euc.dist[i] <- tmp.df %>% filter(chromatin == "Euchromatin") %>% summarise(total = sum(G2.count)) %>% .$total
    het.dist[i] <- tmp.df %>% filter(chromatin == "Heterochromatin") %>% summarise(total = sum(G2.count)) %>% .$total
  }
  # p vals
  cen.result <- sum(abs(cen.dist) >= abs(cent.count))/(n)
  #euc.result <- sum(abs(euc.dist) >= abs(euc.count))/(n)
  het.result <- sum(abs(het.dist) >= abs(het.count))/(n)
  
  return(list(cen.result, cen.dist, het.result, het.dist))
}



permutation.test3 <- function(df, n, cent.count, euc.count){
  # set up variables
  tmp.df <- df
  cen.dist=c()
  euc.dist=c()
  #het.dist=c()
  cen.result=0
  euc.result=0
  #het.result=0
  # scramble loop
  for(i in 1:n){
    tmp.df <- tmp.df %>% mutate(chromatin = sample(chromatin))
    cen.dist[i] <- tmp.df %>% filter(chromatin == "Centromere") %>% summarise(total = sum(G2.count)) %>% .$total
    euc.dist[i] <- tmp.df %>% filter(chromatin == "Euchromatin") %>% summarise(total = sum(G2.count)) %>% .$total
    #het.dist[i] <- tmp.df %>% filter(chromatin == "Heterochromatin") %>% summarise(total = sum(G2.count)) %>% .$total
  }
  # p vals
  cen.result <- sum(abs(cen.dist) >= abs(cent.count))/(n)
  euc.result <- sum(abs(euc.dist) >= abs(euc.count))/(n)
  #het.result <- sum(abs(het.dist) >= abs(het.count))/(n)
  
  return(list(cen.result, cen.dist, euc.result, euc.dist))
}



#### everything, including y chromosome

test1 <- permutation.test(genome.wind %>% filter(chr!="tig00057289"), 10000, cent.count, euc.count, het.count)


#### no Y chromosome, including y chromosome

cent.count.y <- genome.wind %>% filter(chromatin == "Centromere", Ychrom=="other") %>% summarise(total = sum(G2.count)) %>% .$total
euc.count.y <- genome.wind %>% filter(chromatin == "Euchromatin", Ychrom=="other") %>% summarise(total = sum(G2.count)) %>% .$total
het.count.y <- genome.wind %>% filter(chromatin == "Heterochromatin", Ychrom=="other") %>% summarise(total = sum(G2.count)) %>% .$total


test2 <- permutation.test(genome.wind %>% filter(Ychrom=="other", chr!="tig00057289"), 10000, cent.count.y, euc.count.y, het.count.y)



########################
#### Short read permutation test
########################



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


## centromere TEs to keep
#cent.te.to.keep <- c("BS","DOC","DOC2_DM","DOC6_DM","FW_DM","G5_DM","Jockey_3_DM","LINEJ1_DM","NOMAD_I","PROTOP_A","R1_DM","ROO_I","Jockey_3_DM")

## rare TEs are "de novo" or young
te.rare <- DGRP.pop.fig %>% filter(freq.class=="Rare", ID %in% cent.te.to.keep) 



#### filter out de novo


jock.bed.denovo.all <-  te.rare %>% #filter(denovo == "yes", keep.delete=="keep", ref.contig!="unknown") %>% 
  select(contig, start, chromatin.state, region, ID) %>% 
  mutate(ref.end = start + 1) %>% 
  rename(chr = contig) %>% rename(ref.start = start, ref.chromatin = chromatin.state, ref.region = region) %>% 
  select(chr, ref.start, ref.end, ref.chromatin, ref.region, ID) 

  
jock.bed.denovo <- jock.bed.denovo.all %>% filter(ID=="Jockey_3_DM")
#jock.bed.denovo <- jock.bed.denovo.all %>% filter(ID=="LINEJ1_DM")
#jock.bed.denovo <- jock.bed.denovo.all %>% filter(ID=="DOC")



## divide into regions
jock.cent.denovo <- jock.bed.denovo %>% filter(ref.region == "centromere")
jock.euc.denovo <- jock.bed.denovo %>% filter(ref.chromatin == "eu")
jock.het.denovo <- jock.bed.denovo %>% filter(ref.chromatin == "het", ref.region != "centromere")


#### split up genome into 100Kb tiles

cent.spans <- chrom.spans %>% filter(region=="centromere") %>% select(contig_name, start.range, end.range) %>% mutate(chromatin = "Centromere")
euc.spans <- chrom.spans %>% filter(chromatin=="eu") %>% select(contig_name, start.range, end.range)
het.spans <- chrom.spans %>% filter(chromatin=="het", region!="centromere") %>% select(contig_name, start.range, end.range)


## convert to genomic ranges
cent.range <- toGRanges(cent.spans)
euc.range <- toGRanges(euc.spans)
het.range <- toGRanges(het.spans)


## tile up het and euchromatin

euc.tiles <- tile(euc.range, width = 90000)
het.tiles <- tile(het.range, width = 90000)

tmp.euc.tiles <- as.data.frame(euc.tiles)
tmp.het.tiles <- as.data.frame(het.tiles)

cent.wind <- cent.spans %>% rename(chr = contig_name, start = start.range, end = end.range)
euc.wind <- tmp.euc.tiles %>% select(seqnames, start, end) %>% rename(chr = seqnames) %>% mutate(chromatin = "Euchromatin")
het.wind <- tmp.het.tiles %>% select(seqnames, start, end) %>% rename(chr = seqnames) %>% mutate(chromatin = "Heterochromatin")


#### overlaps

## centromere

x <- toGRanges(cent.wind)
y <- toGRanges(jock.cent.denovo)

overlaps <- countOverlaps(x, y)
cent.wind$G2.count <- as.vector(overlaps)

## euchromatin

x <- toGRanges(euc.wind)
y <- toGRanges(jock.euc.denovo)

overlaps <- countOverlaps(x, y)
euc.wind$G2.count <- as.vector(overlaps)

## heterochromatin

x <- toGRanges(het.wind)
y <- toGRanges(jock.het.denovo)

overlaps <- countOverlaps(x, y)
het.wind$G2.count <- as.vector(overlaps)


#### combine altogether

genome.wind <- rbind.data.frame(cent.wind, euc.wind, het.wind)

genome.wind <- genome.wind %>% mutate(Ychrom = if_else(str_detect(chr, "Y"), "Y", "other"))


## jockey-3 counts


cent.count <- genome.wind %>% filter(chromatin == "Centromere") %>% summarise(total = sum(G2.count)) %>% .$total
euc.count <- genome.wind %>% filter(chromatin == "Euchromatin") %>% summarise(total = sum(G2.count)) %>% .$total
het.count <- genome.wind %>% filter(chromatin == "Heterochromatin") %>% summarise(total = sum(G2.count)) %>% .$total


#### test 

test3 <- permutation.test(genome.wind %>% filter(Ychrom=="other", chr!="tig00057289"), 10000, cent.count, euc.count, het.count)


test4 <- permutation.test2(genome.wind %>% filter(Ychrom=="other", chr!="tig00057289", chromatin != "Euchromatin"), 10000, cent.count, het.count)


test5 <- permutation.test3(genome.wind %>% filter(Ychrom=="other", chr!="tig00057289", chromatin != "Heterochromatin"), 10000, cent.count, euc.count)


