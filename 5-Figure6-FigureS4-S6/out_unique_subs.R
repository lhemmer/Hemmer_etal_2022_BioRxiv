######## 
#### Getting unique subs for each sequence in an alignment
######## 

#### Import libraries

library(BiocManager)
#BiocManager::install("Biostrings")
library(Biostrings)
library(tidyr)
library(dplyr)

#### import alignment data

args = commandArgs(trailingOnly=TRUE)

inFile <- args[1]

outFile <- args[2]

setwd("/Users/lucashemmer/Documents/mel_centromere/age_test/DGRP3/")

teAlign <- readDNAStringSet(filepath = inFile,format="fasta")

#teAlign <- readDNAStringSet(filepath = "/Users/lucashemmer/Downloads/alignment_iso1_Jockey-3_DM_all_noDsim.fasta",format="fasta")
#jockeyAlign <- readDNAStringSet(filepath = "/Users/lucashemmer/Documents/mel_centromere/DGRP2/alignment_all_jockey-3_DM_DmelRef_DGRPcalls_wRef.fasta",format="fasta")
names(teAlign)

freq.list <- list()
for (i in 1:width(teAlign[1])) {
  freq.list[[i]] <- as.data.frame(table(substr(teAlign, i, i)))
}

## does not count indels
for (i in 1:length(freq.list)) {
  if ("-" %in% freq.list[[i]]$Var1 == F) {
    print(i)
  }
}
pos <- vector()
nuc <- vector()
for (i in 1:length(freq.list)) {
  if (length(freq.list[[i]]$Freq) > 2) {
    for (j in 1:length(freq.list[[i]]$Freq)) {
      if (freq.list[[i]]$Freq[j] == 1 & freq.list[[i]]$Var1[j] != '-') {
        pos <- c(pos, i)
        nuc <- c(nuc, as.character(freq.list[[i]]$Var1[j]))
      }
    }
  }
}


numb.vect <- vector(length = length(teAlign))
numb.vect <- rep(0, length(teAlign))

for (i in 1:length(pos)) {
  tmp.pos <- substr(teAlign, pos[i], pos[i])
  for (j in 1:length(teAlign)) {
    if (tmp.pos[j] == nuc[i]) {
      numb.vect[j] <- numb.vect[j] + 1
    }
  }
}


numb.vect <- as.data.frame(numb.vect)
numb.vect$TE <- names(teAlign)

numb.vect1 <- numb.vect %>% separate(col = TE, into = c("contig","span","direction."), sep = "\\|", remove=F) %>%
  separate(col = span, into=c("start","end"), sep = "-",remove = T) #mutate(mid.var1=end.var1-start.var1) %>%

write.table(numb.vect1, file=outFile, quote = F, row.names = F, sep="\t")

#write.table(g2.denovo, file=paste(pwdO,"DGRP3_singleton_fisherResults_Yremoved_20x.txt",sep = ""), quote = F, row.names = F, sep="\t")

