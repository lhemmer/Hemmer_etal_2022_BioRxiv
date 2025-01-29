#####################
#### Coverage analysis of DGRP, GDL, and DSPR
#####################



## load libraries

library(DescTools)
library(viridis)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
options(scipen=999)


## load data
DGRP.cent.cov <- read.table("Table_S3_DGRP_total_cent_cov_500bpWin.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DSPR.cent.cov <- read.table("DSPR_total_cent_cov.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
GDL.cent.cov <- read.table("GDL_total_cent_cov.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)

###$ add source
DGRP.cent.cov$sources <- "DGRP"
DSPR.cent.cov$sources <- "DSPR"
GDL.cent.cov$sources <- "GDL"

cencov <- rbind.data.frame(DGRP.cent.cov, DSPR.cent.cov, GDL.cent.cov)

#pdf(paste(pwdO,"Figure_ST1_DGRP3_GDL_DSPR_coverage_EuCen.pdf",sep = ""),width=10,height=6)
ggplot(cencov, aes(x=centromere, y=med_cov)) + geom_boxplot() + theme_bw() +
  ylab("Median Genomic Coverage per Sample") + xlab("Region") + scale_x_discrete(labels= c("Centromere", "Euchromatin")) +
  facet_grid(cols=vars(sources))
#dev.off()



#### now per contig

DGRP.contig.cov <- read.table("DGRP2_total_contig_cov_500bpWin.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
DSPR.contig.cov <- read.table("DSPR_total_contig_cov.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
GDL.contig.cov <- read.table("GDL_total_contig_cov.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)

###$ add source
DGRP.contig.cov$sources <- "DGRP"
DSPR.contig.cov$sources <- "DSPR"
GDL.contig.cov$sources <- "GDL"

contigcov <- rbind.data.frame(DGRP.contig.cov, DSPR.contig.cov, GDL.contig.cov)
contigcov$contig <- factor(contigcov$contig, levels=c("X_1","Contig79","2L_1","tig00057289","2R_21","3L_1","3R_5","3R_28","4_2","Contig119"))

#pdf(paste(pwdO,"Figure_ST2_DGRP3_GDL_DSPR_coverage_perContig.pdf",sep = ""),width=10,height=10)
ggplot(contigcov, aes(x=contig, y=med_cov)) + geom_boxplot() + theme_bw() +
  ylab("Median Genomic Coverage per Sample") + xlab("Chromosome / Region") + 
  scale_x_discrete(labels= c("X", "CenX", "2L", "Cen2", "2R", "3L","Cen3","3R","4","Cen4")) +
  facet_grid(rows=vars(sources))
#dev.off()
