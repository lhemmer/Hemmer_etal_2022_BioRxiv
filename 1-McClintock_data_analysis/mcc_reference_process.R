#!/usr/bin/env Rscript

## Load Libraries

library(DescTools)

## Load data

#setwd("/Users/lucashemmer/Desktop/test")
dat.ref <- read.table("dmel.chromosomes.fa.TE.mcClintock.2020.mod.gff", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
sizes <- read.csv(paste(pwd,"dmel_scaffold2_plus0310_sizes.csv",sep=""), header=FALSE)
colnames(sizes) <- c("contig","len")


### load dgrp data

dat <- read.table("DGRP_calls.txt", header=TRUE, sep="\t",stringsAsFactors = F)
dat <- dat[!(dat$called.by=="temp" & dat$ref=="reference"),]
print(nrow(dat))

### remove non-reference calls

dat1 <- dat[dat$ref=="non-reference",]
dat2 <- dat[dat$ref=="reference",]
remove(dat)
remove(dat1)

### remove Y chromosome calls

rc.by.chrom <- dat2[-grep("Y",dat2$chromosome),]
remove(dat2)

### identify the repeats and contigs which were detected to reduce run time

rc.contig <- unique(rc.by.chrom$contig)
rc.repeats <- unique(rc.by.chrom$ID)
ref.by.chrom <- ref.by.chrom[ref.by.chrom$contig %in% rc.contig,]
ref.by.chrom <- ref.by.chrom[ref.by.chrom$repeats %in% rc.repeats,]

### order and add identifier column

rc.by.chrom <- rc.by.chrom[order(rc.by.chrom$contig,rc.by.chrom$start,rc.by.chrom$ID),] #rc.by.chrom$population,
rc.by.chrom$identifier <- NA

### split the dataframe into a list of dataframes for faster looping

rc.by.chrom <- split(rc.by.chrom,f=rc.by.chrom$contig)
ref.by.chrom <- ref.by.chrom[order(ref.by.chrom$contig),]
ref.by.chrom <- split(ref.by.chrom,f=ref.by.chrom$contig)


### loop to find overlaps between reference TE position and calls

for (k in 1:length(rc.by.chrom)) {
  print(paste("TE # ",k))
  for (i in 1:nrow(ref.by.chrom[[k]])) {
    for (j in which(rc.by.chrom[[k]]$ID == ref.by.chrom[[k]]$repeats[i] & is.na(rc.by.chrom[[k]]$identifier))) {
      if ( c(ref.by.chrom[[k]]$start[i],ref.by.chrom[[k]]$end[i])  %overlaps% c(rc.by.chrom[[k]]$start[j],rc.by.chrom[[k]]$end[j]) ) {
        rc.by.chrom[[k]]$identifier[j] <- paste(ref.by.chrom[[k]]$repeats[i],ref.by.chrom[[k]]$contig[i],ref.by.chrom[[k]]$start[i],sep=".")
      }
    } 
    print(i)
  }
}

### combine back into dataframe

rc.by.chrom <- do.call(rbind, rc.by.chrom)

### remove those that did not work and export

rc.to.export <- rc.by.chrom[!is.na(rc.by.chrom$identifier),]

write.table(rc.to.export,"DGRP_ref_calls_processed1.txt",quote = F,sep="\t",row.names = F)
write.table(rc.by.chrom,"DGRP_ref_calls_processed_NAincluded.txt",quote = F,sep="\t",row.names = F)

### taking care of first round of missing values, keep those missing

rc.by.chrom <- rc.by.chrom[is.na(rc.by.chrom$identifier),]

### sort

rc.by.chrom <- rc.by.chrom[order(rc.by.chrom$contig,rc.by.chrom$start,rc.by.chrom$ID),]
rc.contig <- unique(rc.by.chrom$contig)
ref.by.chrom <- ref.by.chrom[ref.by.chrom$contig %in% rc.contig,]
rc.by.chrom <- split(rc.by.chrom,f=rc.by.chrom$contig)
ref.by.chrom <- ref.by.chrom[order(ref.by.chrom$contig),]
ref.by.chrom <- split(ref.by.chrom,f=ref.by.chrom$contig)

### loop 

for (k in 1:length(rc.by.chrom)) {
  print(paste("TE # ",k))
  for (i in 1:nrow(ref.by.chrom[[k]])) {
    for (j in which(is.na(rc.by.chrom[[k]]$identifier) & (ref.by.chrom[[k]]$start[i] > rc.by.chrom[[k]]$start - 1000) ) ) {  
      print(j)
      if ( c(ref.by.chrom[[k]]$start[i],ref.by.chrom[[k]]$end[i])  %overlaps% c(rc.by.chrom[[k]]$start[j],rc.by.chrom[[k]]$end[j]) 
           && Overlap(c(ref.by.chrom[[k]]$start[i],ref.by.chrom[[k]]$end[i]), c(rc.by.chrom[[k]]$start[j],rc.by.chrom[[k]]$end[j])) > 30) {
        rc.by.chrom[[k]]$identifier[j] <- paste(ref.by.chrom[[k]]$repeats[i],ref.by.chrom[[k]]$contig[i],ref.by.chrom[[k]]$start[i],sep=".")
      }
    } 
    print(i)
  }
}

### combine back into dataframe
rc.by.chrom <- do.call(rbind, rc.by.chrom)

### export the new identified ones

rc.to.export <- rc.by.chrom[!is.na(rc.by.chrom$identifier),]

write.table(rc.to.export,"DGRP_ref_calls_NAprocessed1.txt",quote = F,sep="\t",row.names = F)


### dealing with missing one last time

rc.by.chrom <- rc.by.chrom[is.na(rc.by.chrom$identifier),]

### sort

rc.by.chrom <- rc.by.chrom[is.na(rc.by.chrom$identifier),]
rc.contig <- unique(rc.by.chrom$contig)
ref.by.chrom <- ref.by.chrom[ref.by.chrom$contig %in% rc.contig,]
rc.by.chrom <- split(rc.by.chrom,f=rc.by.chrom$contig)
ref.by.chrom <- ref.by.chrom[order(ref.by.chrom$contig),]
ref.by.chrom <- split(ref.by.chrom,f=ref.by.chrom$contig)

### loop
for (k in 1:length(rc.by.chrom)) {
  print(paste("TE # ",k))
  for (i in 1:nrow(ref.by.chrom[[k]])) {
    for (j in which(rc.by.chrom[[k]]$ID == ref.by.chrom[[k]]$repeats[i] & is.na(rc.by.chrom[[k]]$identifier))) {
      if ( c(ref.by.chrom[[k]]$start[i]-500,ref.by.chrom[[k]]$end[i]+500)  %overlaps% c(rc.by.chrom[[k]]$start[j],rc.by.chrom[[k]]$end[j]) ) {
        rc.by.chrom[[k]]$identifier[j] <- paste(ref.by.chrom[[k]]$repeats[i],ref.by.chrom[[k]]$contig[i],ref.by.chrom[[k]]$start[i],sep=".")
      }
    } 
    print(i)
  }
}

pat <- "(.*)_.*"

for (k in 1:length(rc.by.chrom)) {
  print(paste("TE # ",k))
  for (i in 1:nrow(ref.by.chrom[[k]])) {
    for (j in which(is.na(rc.by.chrom[[k]]$identifier) ) ) {
      if ( sub(pat,"\\1",ref.by.chrom[[k]]$repeats[i]) == sub(pat,"\\1",rc.by.chrom[[k]]$ID[j]) ) {
        if (c(ref.by.chrom[[k]]$start[i]-250,ref.by.chrom[[k]]$end[i]+250)  %overlaps% c(rc.by.chrom[[k]]$start[j],rc.by.chrom[[k]]$end[j])) {
          print("yes")
          rc.by.chrom[[k]]$identifier[j] <- paste(ref.by.chrom[[k]]$repeats[i],ref.by.chrom[[k]]$contig[i],ref.by.chrom[[k]]$start[i],sep=".")
        }
      }
    } 
  }
}

### export
rc.to.export <- do.call(rbind, rc.by.chrom)

write.table(rc.to.export,"DGRP_ref_calls_NAprocessed2.txt",quote = F,sep="\t",row.names = F)


