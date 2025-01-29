################
#### Phylogenetic tree of the separate clades
################


#### load libraries


library(ape)
library(phytools)
library(ggtree)

#### load data

jock.tree <- read.tree(file="RAxML_bipartitions.alignment_Jockey-3_Ychrom_all.automre")


#pdf(Figure_S12_RAxML_bipartitions.alignment_Jockey-3_Ychrom_all.automre_phylogeny_highlight.pdf",width=6,height=6)
ggtree(jock.tree) + 
  geom_highlight(node=364, fill="#0D0887FF", alpha=0.5) +
  geom_highlight(node=201, fill="#FBD424FF", alpha=0.5) +
  geom_highlight(node=326, fill="#FCA338FF", alpha=0.5) +
  geom_highlight(node=234, fill="#DE5F65FF", alpha=0.5) +
  geom_highlight(node=258, fill="#A92395FF", alpha=0.5) +
  geom_cladelabel(node=364, label="Class 1\nYoung Insertions", align=T, color='black', offset.text = 0.005) +
  geom_cladelabel(node=201, label="Class 2\nShort 3' Fragments", align=T, color='black', offset.text = 0.005) +
  geom_cladelabel(node=326, label="Class 3\nOlder Internal Fragments", align=T, color='black', offset.text = 0.005) +
  geom_cladelabel(node=234, label="Class 4\nLong Old 3' Fragment", align=T, color='black', offset.text = 0.005) +
  geom_cladelabel(node=258, label="Class 4A\nLong Old 3' Fragments\nfrom Region3 on CenY", align=T, color='black', offset.text = 0.004) + 
  xlim(0,0.18)
#dev.off()


