############
## Doing some DSPR stats 
############

#### load libraries

library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)
library(ggpubr)

#### importing gff for calls 

ref.calls <- read.table("scaffolds.fa.out.gff", header=FALSE, sep="\t",stringsAsFactors = F, fill = F, check.names=F)



ref.calls1 <- ref.calls %>% mutate(V9 = str_remove( V9, "-int")) %>% 
  filter(!str_detect(V9, '\\)n')) %>% filter(!str_detect(V9, '-rich')) %>%
  separate(V9, into=c("Target", "TE", "te.start","te.end"), sep="\\s", remove = T) %>% 
  mutate(TE = str_remove( TE, "Motif:")) %>%
  dplyr::select(-c(V2, V3, V8, Target)) %>% dplyr::rename(contig=V1, contig.start= V4, contig.end = V5, div = V6, direction=V7) %>%
  mutate(te.start = as.numeric(te.start)) %>% mutate(te.end = as.numeric(te.end)) %>% filter(!str_detect(TE, 'SAT')) %>% 
  mutate(TE = str_replace(TE, "-", "_")) %>% filter(!str_detect(contig, "chr")) %>%
  separate(contig, into=c("population","contig"), sep="\\.", remove=T)



tes.to.keep <- c("ACCORD_I","ACCORD_LTR","ACCORD2_I","ACCORD2_LTR","Baggins1","BARI_DM","BARI1",
                 "BATUMI_I","BATUMI_LTR","BEL_I","BEL_LTR","Bica_I","Bica_LTR","BLASTOPIA_I",
                 "BLASTOPIA_LTR","BLOOD_I","BLOOD_LTR","BS","BS2","BS3_DM","BS4_DM","BURDOCK_I",
                 "BURDOCK_LTR","Chimpo_I","Chimpo_LTR","Chouto_I","Chouto_LTR","CIRCE",
                 "Copia_I","Copia_LTR", "Copia1_I_DM","Copia1_LTR_DM","Copia2_I","Copia2_LTR_DM",
                 "DIVER_I","DIVER_LTR","DIVER2_I","DIVER2_LTR","DM1731_I","DM1731_LTR","DM176_I",
                 "DM176_LTR","DM297_I","DM297_LTR","DM412","DM412B_LTR","DMCR1A","DMLTR5","DMRT1A",
                 "DMRT1B","DMRT1C","DNAREP1_DM","DNAREP1_DSim","DOC","DOC2_DM","DOC3_DM","DOC4_DM",
                 "DOC5_DM","DOC6_DM","FB4_DM","FROGGER_I","FROGGER_LTR","FW_DM","FW2_DM","FW3_DM",
                 "G_DM","G3_DM","G4_DM","G5_DM","G5A_DM","G6_DM","G7_DM","Galileo_DB","Galileo_DM",
                 "GTWIN_I","GTWIN_LTR","Gypsy_I","Gypsy_LTR","Gypsy1_I_DM","Gypsy1_LTR_DM","Gypsy11_I",
                 "Gypsy11_LTR","Gypsy12_LTR","Gypsy12_I","Gypsy2_I","Gypsy2_I_DM","Gypsy2_LTR",
                 "Gypsy2_LTR_DM","Gypsy3_I","Gypsy3_LTR","Gypsy4_I","Gypsy4_LTR","Gypsy5_I",
                 "Gypsy5_LTR","Gypsy6_I","Gypsy6_LTR","Gypsy6A_LTR","Gypsy7_I","Gypsy7_LTR",
                 "Gypsy8_I","Gypsy9_I","Gypsy8_LTR","Gypsy9_I","Gypsy9_LTR","HETA","HOBO","I_DM",
                 "IDEFIX_I","IDEFIX_LTR","Jockey_3_Dmel_08212020","Jockey2","LINEJ1_DM","M4DM",
                 "MAX_I","MAX_LTR","MDG1_I","MDG1_LTR","MICROPIA_I","MICROPIA_LTR","NOMAD_I","NOMAD_LTR",
                 "POGO","PROTOP","PROTOP_A","PROTOP_B","QUASIMODO_I","QUASIMODO_LTR","QUASIMODO2_I_DM",
                 "QUASIMODO2_LTR_DM","R1_DM","R2_DM","ROO_I","ROO_LTR","ROVER_I_DM","ROVER_LTR_DM",
                 "S_DM","S2_DM","TABOR_I","TABOR_LTR","TAHRE","TART_A","TART_B1","TC1_2_DM","TC1_DM",
                 "TIRANT_I","TIRANT_LTR","TRANSIB1","TRANSIB2","TRANSIB3","TRANSIB4","TRANSPAC_I","TRANSPAC_LTR","ZAM_I","ZAM_LTR")

pop.to.keep <- c("A1","A2","A3","A4","A5","A6","A7","AB8","B1","B2","B3","B4","B6")  

ref.calls2 <- ref.calls1 %>% filter(TE %in% tes.to.keep) %>% mutate(population = str_to_upper(population)) %>% filter(population %in% unique(all.te.pop$population))
all.te.pop2 <- all.te.pop %>% filter(ID %in% tes.to.keep)  

ref.calls.sum <- ref.calls2 %>% dplyr::count(TE, population, name="long.calls.n") 
all.te.sum <- all.te.pop2 %>% dplyr::count(ID, population, name="mcc.calls.n") %>% dplyr::rename(TE=ID)

te.sum1 <- full_join(ref.calls.sum, all.te.sum)
te.sum1$population <- factor(te.sum1$population, levels=c("A1","A2","A3","A4","A5","A6","A7","AB8","B1","B2","B3","B4","B6"))

te.sum2 <- full_join(ref.calls2 %>% dplyr::count(TE, name="long.calls.n"), all.te.pop2 %>% dplyr::count(ID, name="mcc.calls.n") %>% rename(TE=ID))
te.sum2$TE[te.sum2$TE=="Jockey_3_Dmel_08212020"] <- "Jockey_3_DM"



pdf("Figure_ST3_DSPR_TE_count_longshortCompare_byTE.pdf",width=12,height=12)
p1 <- te.sum2 %>% ggplot(aes(x=long.calls.n, y=mcc.calls.n, label = TE)) +
  geom_point() + #geom_label_repel() +
  geom_smooth() +
  geom_text_repel(aes(label=TE)) +
  theme_bw() + ylab("# Detected by McClintock") + xlab("# Detected in Long-Read Assembly")
#xlim(0,1000) + ylim(0,500)

p2 <- te.sum2 %>% ggplot(aes(x=long.calls.n, y=mcc.calls.n, label = TE)) +
  geom_point() + #geom_label_repel() +
  geom_smooth() +
  geom_text_repel(aes(label=TE)) +
  theme_bw() + ylab("# Detected by McClintock") + xlab("# Detected in Long-Read Assembly") +
  xlim(0,1000) + ylim(0,500)

ggarrange(p1,p2,nrow=2,labels=c("A","B"))
dev.off()
