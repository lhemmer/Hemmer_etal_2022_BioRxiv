## load DGRP data

pwd <- "/Users/lucashemmer/Desktop/test"

## load reference calls 

# this one has NA removed already
rc.by.chrom1 <- read.table(paste(pwd,"DGRP_ref_calls_processed1.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F)
nrow(rc.by.chrom1)
rc.by.chrom1 <- rc.by.chrom1[!is.na(rc.by.chrom1$identifier),]
nrow(rc.by.chrom1)

# this one has NA removed already
rc.by.chrom2 <- read.table(paste(pwd,"DGRP_ref_calls_NAprocessed1.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F)
nrow(rc.by.chrom2)
rc.by.chrom2 <- rc.by.chrom2[!is.na(rc.by.chrom2$identifier),]
nrow(rc.by.chrom2)

# this one has NAs manually corrected
rc.by.chrom3 <- read.table(paste(pwd,"DGRP_ref_calls_NAprocessed2.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F)
nrow(rc.by.chrom3)

## Combine

rc.by.chrom <- rbind.data.frame(rc.by.chrom1,rc.by.chrom2,rc.by.chrom4)
remove(rc.by.chrom1,rc.by.chrom2,rc.by.chrom3)


## now importing  the de novo calls

te.by.chrom <- read.table(paste(pwd,"DGRP_denovo_calls_processed.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F)

## Combine

all.te <- rbind.data.frame(rc.by.chrom, te.by.chrom, rc.by.chrom.low)

## Set up an identifier for each unique McClintock call
all.te$mcc.identifier <- paste(all.te$identifier,all.te$new.pop,sep=".")

## Load library

library(dplyr)

## Number of times a TE is called
all.te2 <- all.te %>% 
  add_count(mcc.identifier, name = "call.count")

all.te2 <- as.data.frame(all.te2)

## How many times does a specific TE appear in the populaiton
all.te2$pop.count <- NA

all.te3 <- split(all.te2,f=all.te2$identifier)

for (k in 1:length(all.te3)) {
  tmp <- length(unique(all.te3[[k]]$population))
  all.te3[[k]]$pop.count <- tmp
}

## Combine
rc.to.export1 <- do.call(rbind, all.te3)

write.table(rc.to.export1,"DGRP_calls_withFreq.txt",quote = F,sep="\t",row.names = F)


