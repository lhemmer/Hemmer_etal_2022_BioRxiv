######## 
#### Filter out TEs based on number of times they are called within McClintock
######## 

#### Load libraries

library(DescTools)
library(viridis)
library(ggplot2)
library(tidyr)
library(dplyr)
options(scipen=999)

#### Load data

DGRP.all <- read.table("DGRP_calls_withFreq.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DGRP.all$ID[DGRP.all$ID=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"


## higher coverage sample calls

DGRP.cent.cov <- read.table(paste(pwd,"DGRP_total_cent_cov_500bpWin.txt",sep=""), header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
cent.hcov.pop <- DGRP.cent.cov %>% filter(centromere=="cent" & med_cov > 20 & population!=340 & population!=373 & population!=138) %>% pull(population) 

## note we had problems with 340, 373, and 138 because they had lots of duplicated calls which seemed abnormal given their coverage


######################################################
#### De novo calls need to either be in more than 1 pop or confirmed by multiple methods
######################################################

######## 
#### 
######## 



DGRP.all.filt <- DGRP.all %>% mutate(centromere=ifelse(region=="centromere","Centromere",
  ifelse(region!="centromere" & chromatin.state=="het","Heterochromatin","Euchromatin"))) %>% ## filter out piRNA clusters for popgen analysis
  filter(!(chromosome=="Y")) %>% #%>% filter(ref=="non-reference") ## Y chromosome calls unrealiable in DGRP
  filter(!(ID %in% c("Gypsy6_I_Dmoj","Gypsy1_I_Dmoj"))) %>% ## repeats misannotated as Gypsy elements from D. mojavensis
  filter(population %in% cent.hcov.pop) %>% 
  filter(identifier!="Jockey_3_Dmel_08212020.2R_19.0") ## independently confirmed these reads were misplaced

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


write.table("File_S2_McClintock_output_DGRP.txt",quote = F,sep="\t",row.names = F)



