## load libraries

library(readr)
library(dplyr)
library(tidyr)

## load data from entire directory

pwd = "/Users/lucashemmer/Documents/mel_centromere/DGRP2/coverage/"
bed_files <- list.files(path = pwd, pattern = "*.bed.gz")

tmp.bed <- read.table(paste(pwd,bed_files[1],sep=""), header=FALSE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)

depth.bed <- data.frame()
for (i in 1:length(bed_files)) {
  tmp.bed <- read.table(paste(pwd,bed_files[i],sep=""), header=FALSE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
  strain.num <- gsub(".bam_cov.regions.bed.gz","",bed_files[i])
  tmp.bed$V5 <- strain.num
  depth.bed <- rbind.data.frame(depth.bed,tmp.bed)
}

colnames(depth.bed) <- c("contig","win.start","win.end","coverage","population")
depth.bed$win.mid <- (depth.bed$win.end + depth.bed$win.start) / 2

### coverage stats

## coverage per sample overall 

mean.cov.pop <- group_by(depth.bed,population) %>%
  summarize(mean_cov = mean(coverage), med_cov=median(coverage))

#write.table(mean.cov.pop, file="DGRP_total_cov.txt", quote = F, row.names = F, sep="\t")

## coverage in euchromatin and centromeres

mean.cov.pop.cent <- mutate(depth.bed,centromere = ifelse(contig=="X_1" | contig=="2L_1" | contig=="2R_21" | contig=="3L_1" | contig=="3R_28" | contig=="4_2", "eu",
    ifelse(contig=="Contig79" & win.end > 5655 & win.start < 48432 |contig=="tig00057289" & win.end > 5561 & win.start < 7864 |
      contig=="3R_5" & win.end > 19152 & win.start < 69541 | contig=="Contig79" & win.end > 5655 & win.start < 48432 |
      contig=="Contig119" & win.end > 29252 & win.start < 69146 | contig=="Y_contig26", "cent", NA))) %>%
  group_by(population, centromere) %>%
  summarize(mean_cov = mean(coverage), med_cov=median(coverage)) %>% drop_na()

#write.table(mean.cov.pop.cent, file="Table_S3_DGRP_total_cent_cov_500bpWin.txt", quote = F, row.names = F, sep="\t")

## coverage for individual contigs

mean.cov.contig <- filter(depth.bed,
                   contig=="X_1" | contig=="2L_1" | contig=="2R_21" | contig=="3L_1" | contig=="3R_28" | contig=="4_2" |
                   contig=="Contig79" & win.end > 5655 & win.start < 48432 |
                   contig=="tig00057289" & win.end > 5561 & win.start < 7864 |
                   contig=="3R_5" & win.end > 19152 & win.start < 69541 |
                   contig=="Contig79" & win.end > 5655 & win.start < 48432 |
                   contig=="Contig119" & win.end > 29252 & win.start < 69146 |
                   contig=="Y_contig26") %>%
  group_by(population, contig) %>%
  summarize(mean_cov = mean(coverage), med_cov=median(coverage))

#write.table(mean.cov.contig, file="DGRP2_total_contig_cov.txt", quote = F, row.names = F, sep="\t")
