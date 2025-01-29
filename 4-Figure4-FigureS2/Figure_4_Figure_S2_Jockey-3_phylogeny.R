########
#### Phylogenetic analysis of G2/Jockey-3 in melanogaster and figures
########

#### load libraries

library(ape)
library(phytools)
library(geiger)
library(evobiR)
library(ggnewscale)
library(ggtree)
library(dplyr)
library(tidyr)
library(tidytree)
library(ggnewscale)
library(ggpubr)
library(ggtree)
options(scipen=999)


#### load data

## original tree
jock.tree <- read.tree(file="File_S6_RAxML_bipartitions.alignment_Dmel_Ref_Denovo_DsimOut_unknownRemoved_Mod_Ordered.automre")


## useful infromation for plotting
jock.info <- read.table("File_S5_Jockey-3_Phylogeny_Info.txt",header=T,sep="\t",stringsAsFactors = F, fill=T)

## check to make sure jockey-3 tree carries same tips as the information 
name.check(jock.tree,jock.info)

## sort information to match tip labels of tree
jock.info <- jock.info[match(jock.tree$tip.label, jock.info$te.name),]


#### sorting data and adding some columns


## presence or absense of complete ORF (not necessarily coding)

jock.info <- jock.info %>% mutate(orf1 = if_else(te.start <= 186 & te.end >= 1469, "Complete", "Incomplete")) %>% 
  mutate(orf2 = if_else(te.start <= 1457 & te.end >= 4156, "Complete", "Incomplete")) %>% 
  mutate(orfs.present = if_else(orf1 == "Complete" | orf2 == "Complete", "ORF1 or ORF2", "Incomplete"))
#jock.info <- jock.info %>% mutate(orfs.present = if_else(orf1 == "complete" & orf2 == "complete", "ORF1 & ORF2", if_else(orf1 == "incomplete" & orf2 == "complete", "ORF2", "Incomplete")))


## prescence of TAGTTT on either 5' or 3' end

jock.info$tag <- "no"
jock.info$tag[jock.info$tsd.5=="yes" & jock.info$tsd.adjacent == "yes" | 
                jock.info$tsd.3.forward=="yes" & jock.info$tsd.adjacent == "yes"|
                jock.info$tsd.3.reverse=="yes" & jock.info$tsd.adjacent == "yes"] <- "yes"

## presence of TAGTTT and whether the corresponding region has rDNA or not
jock.info$tag.nts <- "no"
jock.info$tag.nts[jock.info$tsd.5=="yes" & jock.info$tsd.adjacent == "yes" | 
                    jock.info$tsd.3.forward=="yes" & jock.info$tsd.adjacent == "yes"|
                    jock.info$tsd.3.reverse=="yes" & jock.info$tsd.adjacent == "yes"] <- "yes"
jock.info$tag.nts[ jock.info$tag.nts == "yes" & jock.info$nts=="yes" ] <- "nts"
jock.info$tag.nts[ jock.info$tag.nts == "no" & jock.info$nts=="yes" ] <- "nts.notag"


#### #### #### 
#### plotting for paper
#### #### #### 


## make sure information label and tree label are same for plotting

jock.info.gg <- as_tibble(jock.info)
jock.info.gg <- jock.info.gg %>% mutate(y.chrom = if_else(ref.chrom == "Y", "Y", "other"))


## test to see if NTS / TAGTTT is more associated with young elements over older elements

jock.info.gg %>% dplyr::count(denovo, tag.nts)
## denovo no, nts 17, nts notag 2, tag%nts 7
## denovo yes, nts 14, nts notag 3, tag%nts 9
fisher.test(matrix(c(121, 26, 95, 26), 2, 2)) #p = 0.4425
## not significant, association with rDNA does not differ between young and old



#### dataframes for tree annotation

## centromere
df.cent <- as.data.frame(jock.info.gg %>% select(centromere))
rownames(df.cent) <- jock.tree$tip.label

## young insertions
df.denovo <- as.data.frame(jock.info.gg %>% select(denovo))
rownames(df.denovo) <- jock.tree$tip.label

## presence of tagttt or variants
df.tag <- as.data.frame(jock.info.gg %>% select(tag.nts))
rownames(df.tag) <- jock.tree$tip.label
df.tag[df.tag != "no",] <- "yes"
df.tag[rownames(df.tag) == "Jockey-3_Dsim_dmel-like",] <- "outgroup"

## complete ORF2
df.orf <- as.data.frame(jock.info.gg %>% select(orf2))
rownames(df.orf) <- jock.tree$tip.label
df.orf$orf2[is.na(df.orf$orf2)] <- "outgroup"

## y chromosome
df.y <- as.data.frame(jock.info.gg %>% select(y.chrom))
rownames(df.y) <- jock.tree$tip.label
df.y$y.chrom[is.na(df.y$y.chrom)] <- "other"


### figures

pdf("Figure_4_Jockey-3_DmelRef_Denovo_DsimOut_unknownRemoved_Mod_Fig2_Ychrom_Mod_noRDNA.pdf",width=10,height=8)

circ <- ggtree(jock.tree,layout="rectangular") %<+% jock.info.gg + 
  #geom_point(aes(color=orf2)) + 
  geom_tippoint(aes(color=orf2)) +
  geom_treescale(offset=2, linesize = 0.5, x=0.01, y=20) + #geom_tiplab( align=T, linetype = 'dotted', offset=0.1)
  #scale_color_manual(breaks=c("ORF1 & ORF2","ORF2", "Incomplete"),
  scale_color_manual(breaks=c("Complete","Incomplete", NA),
                     values=c("black", "grey", ""), name="Complete ORF2") #+ geom_tiplab(aes(label=""), align=T, linetype = 'dotted', offset=0)

p1 <- gheatmap(circ, df.denovo, offset=0.0, width=0.125,colnames=F) +
  scale_fill_manual(breaks=c("yes", "no", "outgroup"), labels=c("Young", "Old", ""),
                    values=c("#9C179EFF", "#ED7953FF", "darkgreen"), name="Relative Age") 

p2 <- p1 + new_scale_fill()

p2 <- gheatmap(p2, df.y, offset=0.015, width=0.125, colnames=F) +
  scale_fill_manual(breaks=c("Y", "other"), labels = c("Y", "Other"),
                    values=c("#2C5F2DFF", "#FFE77AFF", ""), name="Y Chrom") 

p3 <- p2 + new_scale_fill()

p3 <- gheatmap(p3, df.cent, offset=0.03, width=0.125,colnames=F) +
  scale_fill_manual(breaks=c("centromere","non-centromere","outgroup"), labels = c("Centromere", "Non-Centromere",""),
                    values=c("#21908CFF","#440154FF", "darkgreen"), name="Region") 

p3

dev.off()



### figure

pdf("Figure_S2_Jockey-3_DmelRef_Denovo_DsimOut_unknownRemoved_Mod_Fig2_Ychrom_Mod.pdf",width=10,height=8)

circ <- ggtree(jock.tree,layout="rectangular") %<+% jock.info.gg + 
  #geom_point(aes(color=orf2)) + 
  geom_tippoint(aes(color=orf2)) +
  geom_treescale(offset=2, linesize = 0.5, x=0.01, y=20) + #geom_tiplab( align=T, linetype = 'dotted', offset=0.1)
  #scale_color_manual(breaks=c("ORF1 & ORF2","ORF2", "Incomplete"),
  scale_color_manual(breaks=c("Complete","Incomplete", NA),
                     values=c("black", "grey", ""), name="Complete ORF2") #+ geom_tiplab(aes(label=""), align=T, linetype = 'dotted', offset=0)

p1 <- gheatmap(circ, df.denovo, offset=0.0, width=0.125,colnames=F) +
  scale_fill_manual(breaks=c("yes", "no", "outgroup"), labels=c("Young", "Old", ""),
                    values=c("#9C179EFF", "#ED7953FF", "darkgreen"), name="Relative Age") 

p2 <- p1 + new_scale_fill()

p2 <- gheatmap(p2, df.y, offset=0.015, width=0.125, colnames=F) +
  scale_fill_manual(breaks=c("Y", "other"), labels = c("Y", "Other"),
                    values=c("#2C5F2DFF", "#FFE77AFF", ""), name="Y Chrom") 

p3 <- p2 + new_scale_fill()

p3 <- gheatmap(p3, df.cent, offset=0.03, width=0.125,colnames=F) +
  scale_fill_manual(breaks=c("centromere","non-centromere","outgroup"), labels = c("Centromere", "Non-Centromere",""),
                    values=c("#21908CFF","#440154FF", "darkgreen"), name="Region") 

p4 <- p3 + new_scale_fill()

p4 <- gheatmap(p4, df.tag, offset=0.045, width=0.125, colnames=F) +
  scale_fill_manual(breaks=c("yes", "no", "outgroup"), labels = c("Present", "Absent",""),
                    values=c("steelblue", "firebrick", "darkgreen"), name="rDNA") 

p4

dev.off()



