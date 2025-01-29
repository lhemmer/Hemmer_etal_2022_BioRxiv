########
#### Age analysis of G2/Jockey-3 in melanogaster and figures
########

#### load libraries

library(ggnewscale)
library(ggtree)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)

options(scipen=999)


#### load data


ref.calls <- read.table("dmel_scaffold2_plus0310_2.fasta.out.G2.info.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)


ref.calls <- ref.calls[!is.na(ref.calls$ref.chromatin),]
ref.calls <- ref.calls[!is.na(ref.calls$uniq.subs),]
ref.calls <- ref.calls %>% mutate(centromere2=ifelse(ref.region=="centromere" & ref.chrom!="Y","Non-Y Cen", ifelse(ref.region=="centromere" & ref.chrom=="Y","CenY", ifelse(ref.region!="centromere" & ref.chromatin=="het","Heterochromatin","Euchromatin"))))
ref.calls <- ref.calls %>% mutate(centromere=ifelse(ref.region=="centromere", "Centromere", ifelse(ref.region!="centromere" & ref.chromatin=="het","Heterochromatin","Euchromatin")))

ref.calls$centromere2 <- factor(ref.calls$centromere2, levels=c("Euchromatin","Heterochromatin","Non-Y Cen", "CenY"))
ref.calls$centromere <- factor(ref.calls$centromere, levels=c("Euchromatin","Heterochromatin","Centromere"))



#### Divergence from the consensus

my_comparisons <- list( c("CenY", "Heterochromatin"), c("Heterochromatin", "Non-Y Cen"))
symnum.args <- list(cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))


p1 <- ref.calls %>% 
  ggplot(aes(y=div, x=centromere2, color=centromere2)) +
  geom_boxplot() + geom_jitter() + 
  theme_bw() + 
  #xlab("") +
  ylab("Divergence (%) from Consensus") +
  theme(legend.position="bottom", legend.title = element_blank(),plot.title = element_text(hjust = 0.5),         
        axis.title.x=element_blank()
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
  ) + scale_color_manual(values=c("#FFEA46FF", "#7C7B78FF", "#00204DFF","#00204DFF")) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, label.y = c(13, 11), method = "wilcox.test", symnum.args = symnum.args) #+

compare_means(div~centromere2, ref.calls, method = "wilcox.test", paired = FALSE)
div.anova <- aov(div ~ centromere2, data = ref.calls)
TukeyHSD(div.anova)

#$centromere2
#diff         lwr       upr     p adj
#Heterochromatin-Euchromatin -0.7338317 -2.47191556 1.0042522 0.6955074
#Non-Y Cen-Euchromatin        0.4537287 -1.38749266 2.2949500 0.9201171
#CenY-Euchromatin             0.1844273 -1.48067365 1.8495282 0.9918223
#Non-Y Cen-Heterochromatin    1.1875603 -0.05941527 2.4345359 0.0683761
#CenY-Heterochromatin         0.9182590 -0.04994212 1.8864600 0.0701156
#CenY-Non-Y Cen              -0.2693014 -1.41235592 0.8737532 0.9292965

wilcox.test(ref.calls$div[ref.calls$centromere2=="CenY"], ref.calls$div[ref.calls$centromere2=="Heterochromatin"])            #p-value = 0.00000001644; CenY mean = 3.989189 > Het mean 3.07093
wilcox.test(ref.calls$div[ref.calls$centromere2=="Non-Y Cen"], ref.calls$div[ref.calls$centromere2=="Heterochromatin"])       #p-value = 0.0001579; Non-Y Cen mean = 4.258491 > Het mean 3.07093
wilcox.test(ref.calls$div[ref.calls$centromere2=="CenY"], ref.calls$div[ref.calls$centromere2=="Euchromatin"])                #p-value = 0.5449; CenY mean = 3.989189 ~= Euch mean 3.804762
wilcox.test(ref.calls$div[ref.calls$centromere2=="Non-Y Cen"], ref.calls$div[ref.calls$centromere2=="Euchromatin"])           #p-value = 0.7325;  Non-Y Cen mean = 4.258491 ~= Euch mean 3.804762
wilcox.test(ref.calls$div[ref.calls$centromere2=="Heterochromatin"], ref.calls$div[ref.calls$centromere2=="Euchromatin"])     #p-value = 0.1076; Het mean 3.07093  ~= Euch mean 3.804762
wilcox.test(ref.calls$div[ref.calls$centromere2=="CenY" | ref.calls$centromere2=="Non-Y Cen"], ref.calls$div[ref.calls$centromere2=="Euchromatin"])     #p-value = 0.5774; Cen mean 4.060199  ~= Euch mean 3.804762
wilcox.test(ref.calls$div[ref.calls$centromere2=="CenY" | ref.calls$centromere2=="Non-Y Cen"], ref.calls$div[ref.calls$centromere2=="Heterochromatin"])     #p-value = 0.00000001008; Cen mean 4.060199  ~= Het mean 3.07093


ref.calls <- ref.calls %>% mutate(te.length = te.end - te.start) %>% mutate(uniq.per.bp = uniq.subs / te.length)
cor.test(ref.calls$div, ref.calls$te.length, method = 'spearman')


my_comparisons <- list( c("Non-Y Cen", "Heterochromatin"), c("CenY", "Euchromatin"), c("Heterochromatin", "Euchromatin"))

p2 <- ref.calls %>% 
  ggplot(aes(y=uniq.per.bp, x=centromere2, color=centromere2)) +
  geom_boxplot() + geom_jitter() + 
  theme_bw() + #ylab("Percentage") +
  ylab("Unique Substitutions / bp") +
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank()
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
  ) +
  scale_color_manual(values=c("#FFEA46FF", "#7C7B78FF", "#00204DFF","#00204DFF")) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, label.y = c(0.03, 0.06, 0.04), method = "wilcox.test", symnum.args = symnum.args)


wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="CenY"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Heterochromatin"])            #p-value = 0.2419; CenY mean = 0.004928424 ~= Het mean 0.007924991
wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="Non-Y Cen"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Heterochromatin"])       #p-value = 0.000006741; Non-Y Cen mean = 0.01223649 > Het mean 0.007924991
wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="CenY"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Euchromatin"])                #p-value = 0.0169; CenY mean = 0.004928424 > Euch mean 0.006092076
wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="Non-Y Cen"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Euchromatin"])           #p-value = 0.07698;  Non-Y Cen mean = 0.01223649 ~= Euch mean 0.006092076
wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="Heterochromatin"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Euchromatin"])     #p-value = 0.00206; Het mean 0.007924991 < Euch mean 0.006092076
wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="CenY" | ref.calls$centromere2=="Non-Y Cen"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Euchromatin"])     #p-value = 0.213; Cen mean 0.006855426  ~= Euch mean 0.006092076
wilcox.test(ref.calls$uniq.per.bp[ref.calls$centromere2=="CenY" | ref.calls$centromere2=="Non-Y Cen"], ref.calls$uniq.per.bp[ref.calls$centromere2=="Heterochromatin"])     #p-value = 0.0116; Cen mean 0.006855426  < Het mean 0.007924991




#pdf("Figure_S5_Jockey-3_DmelRef_divergence_uniquesubs_boxplot.pdf",width=6,height=8)

ggarrange(p1,p2,nrow=2,ncol=1,labels=c("A","B"), heights = c(1,1.2), common.legend = F)

#dev.off()




#### just centromere, no Y chromosome separation

ref.calls <- ref.calls  %>% mutate(cent=ifelse(ref.region=="centromere","Centromere", ifelse(ref.region!="centromere" & ref.chromatin=="het","Heterochromatin","Euchromatin")))
ref.calls$cent <- factor(ref.calls$cent, levels=c("Euchromatin","Heterochromatin","Centromere"))


my_comparisons <- list( c("Heterochromatin", "Centromere"), c("Euchromatin", "Heterochromatin"), c("Euchromatin", "Centromere"))

p6a <- ref.calls %>% 
  ggplot(aes(y=div, x=cent, color=cent)) +
  geom_boxplot() + geom_jitter() +
  theme_bw() + 
  ylab("Divergence (%) from Consensus") +
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  scale_color_manual(values=c("#FFEA46FF", "#7C7B78FF", "#00204DFF")) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, label.y = c(11, 15, 12.5), method = "wilcox.test", symnum.args = symnum.args)


wilcox.test(ref.calls$div[ref.calls$cent=="Heterochromatin"], ref.calls$div[ref.calls$cent=="Euchromatin"])     #p-value = 0.1076; Het mean 3.07093 ~= Euch mean 3.804762
wilcox.test(ref.calls$div[ref.calls$cent=="Centromere"], ref.calls$div[ref.calls$cent=="Euchromatin"])     #p-value = 0.5774; Cen mean 4.060199  ~= Euch mean 3.804762
wilcox.test(ref.calls$div[ref.calls$cent=="Centromere"], ref.calls$div[ref.calls$cent=="Heterochromatin"])     #p-value = 0.00000001008; Cen mean 4.060199  > Het mean 3.07093


#pdf("/Users/lucashemmer/Documents/Figure_6A_Jockey-3_DmelRef_divergence_uniquesubs_centcombined.pdf",width=6,height=3)
p6a
#dev.off()


#### Other TE divergence

## just go ahead and load the entire GFF file
ref.calls <- read.table("dmel_scaffold2_plus0310_2.fasta.out.gff", header=FALSE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)

## process
ref.calls1 <- ref.calls %>% separate(V9, into=c("Target", "TE", "te.start","te.end"), sep="\\s", remove = T) %>%
  mutate(TE = str_replace(TE, "Motif:", "")) %>% filter(!(str_detect(TE, '\\)n')))  %>% 
  dplyr::select(-c(V2, V3, V7, V8)) %>% mutate(strain = "ISO1") %>%
  dplyr::rename(contig="V1", contig.start = "V4", contig.end = "V5", div = "V6") %>%
  mutate(te.start = as.numeric(te.start)) %>% mutate(te.end = as.numeric(te.end)) %>%
  mutate(TE = str_replace(TE, "-int", ""))

ref.calls2 <- ref.calls1


## select out TEs that we want that are found in the centromer

ref.calls.select <- ref.calls2 %>% dplyr::filter(TE %in% c("BS", "DMRT1B", "DOC", "DOC2_DM", "G5_DM", "NOMAD_I", "Transib5"))



#### preparing information for sorting chromatin status

## load chromatin ranges
chrom.spans <- read.table("dmel_scaffold2_plus0310_chromatin.txt", header=TRUE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)
chrom.spans <- as.data.frame(chrom.spans %>% mutate(centromere=ifelse(region=="centromere","Centromere",
  ifelse(region!="centromere" & chromatin=="het","Heterochromatin","Euchromatin"))))

##process
'%!in%' <- function(x,y)!('%in%'(x,y))
ref.calls.select$ref.chrom <- NA
ref.calls.select$ref.chromatin <- NA
ref.calls.select$ref.region <- NA
ref.calls.select$centromere <- NA
## adding information on chromosome, chromatin, and region
for (i in 1:nrow(ref.calls.select)) {
  print(i)
  if (ref.calls.select$contig[i] %!in% chrom.spans$contig_name) { next }
  tmp.spans <- chrom.spans[chrom.spans$contig_name == ref.calls.select$contig[i],]
  print(i)
  for (j in 1:nrow(tmp.spans)) {
    if (c(ref.calls.select$contig.start[i] >= tmp.spans$start.range[j] & ref.calls.select$contig.start[i] <= tmp.spans$end.range[j])) {
      ref.calls.select$ref.chrom[i] <- tmp.spans$chromosome[j]
      ref.calls.select$ref.chromatin[i] <- tmp.spans$chromatin[j]
      ref.calls.select$ref.region[i] <- tmp.spans$region[j]
      ref.calls.select$centromere[i] <- tmp.spans$centromere[j]
    }
  }
}

ref.calls.select <- ref.calls.select %>% mutate(centromere2=ifelse(ref.chromatin=="eu", "Euchromatin", ifelse(ref.chromatin=="het" & ref.chrom=="Y","Y Chromosome", "Heterochromatin")))


## split into two different dataframes for plotting

ref.calls.select.noncent <- ref.calls.select %>% dplyr::filter(centromere!="Centromere", div <= 20.0)
ref.calls.select.cent <- ref.calls.select %>% dplyr::filter(centromere=="Centromere", div <= 20.0)


## factor to match dataframes together

ref.calls.select.noncent$centromere2 <- factor(ref.calls.select.noncent$centromere2, levels=c("Euchromatin","Heterochromatin","Y Chromosome"))
ref.calls.select.cent$centromere2 <- factor(ref.calls.select.cent$centromere2, levels=c("Euchromatin","Heterochromatin","Y Chromosome"))
ref.calls.select.noncent$TE <- factor(ref.calls.select.noncent$TE, levels=c("BS","DMRT1B","DOC","DOC2_DM","G5_DM","NOMAD_I","Transib5"))
ref.calls.select.cent$TE <- factor(ref.calls.select.cent$TE, levels=c("BS","DMRT1B","DOC","DOC2_DM","G5_DM","NOMAD_I","Transib5"))



ref.calls.plot <- ref.calls.select %>% dplyr::filter(div <= 20.0) %>% filter(!is.na(centromere))
ref.calls.plot$centromere <- factor(ref.calls.plot$centromere, levels=c("Euchromatin","Heterochromatin","Centromere"))
ref.calls.plot$TE <- factor(ref.calls.plot$TE, levels=c("BS","DMRT1B","DOC","DOC2_DM","G5_DM","NOMAD_I","Transib5"))

#pdf("Figure_S4_OtherTEs_DmelRef_divergence_boxplot.pdf",width=5,height=10)
ref.calls.plot %>%
  ggplot(aes(y=div, color=centromere, x=centromere)) +
  geom_boxplot() + geom_jitter() +
  facet_grid(rows=vars(TE), scales = "free_y") +
  theme_bw() + #ylab("Count") + 
  ylab("Divergence (%)") +
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank()
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
  ) +
  scale_color_manual(values=c("#FFEA46FF", "#7C7B78FF", "#00204DFF"))
#dev.off()


## original with het / y chromosome separate
#pdf("Figure_S4_OtherTEs_DmelRef_divergence_boxplot.pdf",width=5,height=10)
ref.calls.select.noncent %>%
  ggplot(aes(y=div, color=centromere2, x=centromere2)) +
  geom_boxplot() + geom_jitter() +
  facet_grid(rows=vars(TE), scales = "free_y") +
  theme_bw() +
  ylab("Divergence (%)") +
  theme(legend.position="none", legend.title = element_blank(),plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank()
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
        ) +
  scale_color_manual(values=c("#FFEA46FF", "#7C7B78FF", "#7C7B78FF")) +
  geom_jitter(data=ref.calls.select.cent, aes(y=div, x=centromere2), colour="#00204DFF", shape=17, size=3)
#dev.off()










#### #### 
#### age of allele test
#### #### 


#### load data


age.all <- read.csv("DGRP3_all_math_output.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F)

row.names(age.all) <- 1:nrow(age.all)


#### calcaulate p value of observing as many copies or fewer in the population

start.index <- which(colnames(age.all) == "prob1")#+1
end.index <- ncol(age.all)

age.all.calc <- age.all

## calculating p values for observing as many or fewer in the population
age.all.calc$p.low <- NA
## calculating p values for observing as many or greater in the population
age.all.calc$p.high <- NA

## note, this is assuming the TE is observed in the original reference sequence

for (i in 1:nrow(age.all.calc)) {
  age.all.calc$p.low[i] <- sum(age.all.calc[i,start.index:(start.index+age.all.calc$NC_presence_highcov[i])]) # as many or fewer copies
  age.all.calc$p.high[i] <- sum(age.all.calc[i,(start.index+age.all.calc$NC_presence_highcov[i]):end.index]) # as many or more copies
}


#### calculating the expected copy number equal to the sum probability of observing the number 

age.all.calc$exp.num <- NA

for (i in 1:nrow(age.all.calc)) {
  tmp.sum <- 0
  for (j in start.index:end.index) {
    tmp.sum1 <- age.all.calc[i,j]
    
    obs.num <- colnames(age.all.calc[j])
    obs.num <- as.numeric(str_remove(obs.num,"prob"))
    
    tmp.sum1 <- tmp.sum1*(obs.num)
    tmp.sum <- sum(tmp.sum,tmp.sum1)
  }
  age.all.calc$exp.num[i] <- tmp.sum
}

## and the difference between observed and expected

age.all.calc$dif <- (age.all.calc$NC_presence_highcov+1) - age.all.calc$exp.num



#### doing some stats, removing the a lot of the probablility table and arranging the data for plotting

cols.keep <- colnames(age.all.calc[c(1:start.index,end.index:ncol(age.all.calc))])
non.class <- colnames(age.all.calc[c(1:start.index-1)])
classes <- colnames(age.all.calc[c(start.index,end.index:ncol(age.all.calc))])

age.all.fig <- dplyr::select(age.all.calc, cols.keep) %>% 
  gather(classes, key="prob.class", -non.class) %>%
  arrange(identifier) %>%
  mutate(obs = str_replace(prob.class, "prob", "")) %>%
  mutate(divergence = numPS/totalPos)


## add centromere, heterochromatin, and euchromatin as factors

age.all.calc.dif <- age.all.calc %>% 
  mutate(dif = (NC_presence_highcov+1) - exp.num) %>% 
  mutate(centromere=ifelse(region=="centromere","Centromere", ifelse(region!="centromere" & chromatin.state=="het","Heterochromatin","Euchromatin"))) %>%
  arrange(dif)


age.all.calc.dif$identifier <- factor(age.all.calc.dif$identifier, levels = age.all.calc.dif$identifier)
age.all.calc.dif$centromere <- factor(age.all.calc.dif$centromere, levels = c("Euchromatin", "Heterochromatin", "Centromere"))


#pdf(paste("/Users/lucashemmer/Documents/mel_centromere/age_test/DGRP3/figures/DGRP3_ageTest_barplot_expVSobs_Jockey3_labeled.pdf",sep = ""),width=6,height=3)

p6b <- ggplot(age.all.calc.dif %>% filter(class == "Jockey-3_DM", !(identifier %in% c("Jockey-3_DM.3R_5|67216-68648|+|","Jockey-3_DM.Contig79|5655-7419|+|","Jockey-3_DM.3R_5|19174-20286|-|")))) + ## These weren't detected in the DGRP
  geom_bar(aes(x=identifier, y = dif, fill = centromere), stat = "identity") +
  coord_flip()+
  scale_fill_manual(values=c("#FFEA46FF", "#7C7B78FF", "#00204DFF"), name="Region") +
  theme_bw() + 
  theme(panel.grid.major.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_rect(fill = NA, size = 0, color = "white", linetype = "blank"),
        panel.border = element_blank(),
        legend.position = "top") +
  ylim(c(-32,32)) + xlab("Unique Insertions") + ylab("Observed - Expected") + scale_x_discrete(expand=c(0, 1)) +
  geom_text(aes(label=ifelse(p.low < 0.10, "*", ""), x = identifier, y = dif - 1)) +
  geom_text(aes(label=ifelse(p.high < 0.10, "*", ""), x = identifier, y = dif + 1))

#dev.off()


### Combine 6A and 6B

#pdf("Figure_6_DGRP3_ageTest_boxDiverg_expVSobs_Jockey3_labeled.pdf",width=6,height=6, pointsize = 12)
ggarrange(p6a,p6b,nrow=2,ncol=1,labels=c("A","B"), heights = c(1,1))
#dev.off()



## doing this for the other four transposable elements

other.te <- c("BS","G5_DM","DOC","DOC2_DM")
dif.plot.list <- list()

for (i in 1:length(other.te)) {
  dif.plot.list[[i]] <- age.all.calc.dif %>% filter(class == other.te[i]) %>%
    arrange(desc(dif)) %>%
    ggplot() +
    geom_bar(aes(x=reorder(identifier, dif), y = dif, fill = centromere), stat = "identity") +
    coord_flip()+
    scale_fill_viridis(discrete = T, option = "E", direction = 1,
                       labels = c("Centromere","Heterochromatin","Euchromatin"), name="Region") +
    theme_bw() + 
    theme(panel.grid.major.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          strip.background = element_rect(fill = NA, size = 0, color = "white", linetype = "blank"),
          panel.border = element_blank(),
          legend.position = "top") +
    ylim(c(-28,28)) + xlab("Unique Insertions") + ylab("Observed - Expected") +
    geom_text(aes(label=ifelse(p.low < 0.10, "*", ""), x = identifier, y = dif - 1)) +
    geom_text(aes(label=ifelse(p.high < 0.10, "*", ""), x = identifier, y = dif + 1))
}

#pdf(paste("Figure_S4_DGRP3_ageTest_barplot_expVSobs_others_labelled.pdf",sep = ""),width=12,height=12)
ggarrange(dif.plot.list[[1]],dif.plot.list[[2]],dif.plot.list[[3]],dif.plot.list[[4]],
          ncol=2,nrow=2,common.legend = TRUE, labels=other.te)
#dev.off()















