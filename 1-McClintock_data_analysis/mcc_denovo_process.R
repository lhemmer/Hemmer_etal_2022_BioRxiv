#!/usr/bin/env Rscript

library(DescTools)

### load dgrp data

#setwd("/Users/lucashemmer/Desktop/test")
dat <- read.table("DGRP_calls.txt", header=TRUE, sep="\t",stringsAsFactors = F)
dat <- dat[!(dat$called.by=="temp" & dat$ref=="reference"),]

dat1 <- dat[dat$ref=="non-reference",]
dat2 <- dat[dat$ref=="reference",]
remove(dat)
remove(dat2)

## Remove Y CHromosome calls

te.by.chrom <- dat1[-grep("Y",dat1$chromosome),]

## Order and set up ID field

te.by.chrom <- te.by.chrom[order(te.by.chrom$ID,te.by.chrom$contig,te.by.chrom$start),]
te.by.chrom$identifier <- NA

## Sliding window size set at 500 bp

wind <- 500

for (i in 1:(nrow(te.by.chrom)-1)) {
  if (is.na(te.by.chrom$identifier[i])) {
    te.by.chrom$identifier[i] <- paste(te.by.chrom$ID[i],te.by.chrom$contig[i],te.by.chrom$start[i],sep=".")
    for (j in 1:500 ) {
      if (i+j == nrow(te.by.chrom)) {break}
      else if (te.by.chrom$ID[i] != te.by.chrom$ID[i+j] | te.by.chrom$contig[i] != te.by.chrom$contig[i+j]) {break}
      else if (te.by.chrom$start[i+j] - te.by.chrom$start[i] < wind | te.by.chrom$end[i+j] - te.by.chrom$end[i] < wind ) {
        te.by.chrom$identifier[i+j] <- te.by.chrom$identifier[i] 
        } else if (te.by.chrom$start[i+j] - te.by.chrom$start[i] > wind & te.by.chrom$end[i+j] - te.by.chrom$end[i] > wind ) {break}
      }
    } 
  print(i)
}


write.table(te.by.chrom,"DGRP_denovo_calls_processed.txt",quote = F,sep="\t",row.names = F)

