
################
#### Phylogenetic tree of the separate clades
################

#### Load libraries


library(ape)
library(phytools)
library(geiger)
library(evobiR)
library(reshape2)


#### load data


## load pairwise distance matrix outputted from MEGAX
g2.mat1 <- read.csv("Y_chrom_pairwise_dist1.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)
g2.mat2 <- read.csv("Y_chrom_pairwise_dist2.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)

#### formatting data for use

## add the two distance matrices together into one 

## function
combine.matrix <- function(mat1, mat2) {
  out.mat <- mat1
  for (i in 1:nrow(out.mat)) {
    for (j in 1:ncol(out.mat)) {
      if (i ==j) {
        out.mat[i,j] <- 0.0
      }
      else if (is.na(out.mat[i,j])) {
        out.mat[i,j] <- mat2[i,j]
      } 
    }
  }
  return(out.mat)
}

g2.mat <- g2.mat1

for (i in 1:nrow(g2.mat)) {
  for (j in 1:ncol(g2.mat)) {
    if (i ==j) {
      g2.mat[i,j] <- 0.0
    }
    else if (is.na(g2.mat[i,j])) {
      g2.mat[i,j] <- g2.mat2[i,j]
    } 
  }
}


ncol(g2.mat)

#### pairwise

g2.iso.pair <- g2.mat[grepl("Y_Contig26",row.names(g2.mat)),grepl("Y_Contig26",colnames(g2.mat))]
g2.iso.pair <- base::as.matrix(g2.iso.pair)
#g2.iso.melt <- melt(g2.iso.pair)

g2.names <- as.data.frame(base::colnames(g2.iso.pair))
colnames(g2.names) <- c("te.name")

g2.names <- g2.names %>% separate(col = te.name, into = c("species.var1","contig.var1","span","direction.var1"), sep = "\\|", remove=F) %>% 
  separate(col = span, into=c("start.var1","end.var1"), sep = "-",remove = T) %>%
  mutate(start.var1=as.numeric(start.var1)) %>%
  arrange(start.var1)


g2.iso.pair.plot <- g2.iso.pair[match(g2.names$te.name,row.names(g2.iso.pair)),match(g2.names$te.name,colnames(g2.iso.pair))]

#### TE annotations

## to make the ggplot process easier for levels in case we wanted to stagger the annotations
jock.info <- jock.info %>% mutate(level=if_else(te.classification == "other.te" | te.classification == "R2", 6, 
  if_else(te.classification == "other.rdna" | te.classification == "28S" | te.classification == "8S", 5,
  if_else(te.classification == "Clade4" | te.classification == "Clade4A", 4,
  if_else(te.classification == "Clade3", 3,
  if_else(te.classification == "Clade2", 2, 1))))))

## level color
jock.info <- jock.info %>% mutate(lev.col = if_else(level == 6, "grey", if_else(level == 5, "#000000",
  if_else(level == 1, "#0D0887FF", if_else(level == 2, "#FBD424FF",
  if_else(level == 3, "#FCA338FF", if_else(te.classification == "Clade4", "#DE5F65FF", "#A92395FF")))))))


## for arrow directions
jock.info <- jock.info %>% mutate(xstart = if_else(strand=="+", contig.start, contig.end)) %>%
  mutate(xend = if_else(strand=="+", contig.end, contig.start))


## for combining  
anno.plot <- jock.info %>% 
  ggplot(aes(x=xstart, xend=xend, y=level, yend=level)) +
  geom_segment(colour=jock.info$lev.col,arrow=arrow(length = unit(0.1, "cm"))) +
  theme_nothing() + xlim(c(0,139957)) +
  theme(legend.position = "none", axis.text=element_blank(), axis.ticks=element_blank(),
        axis.title=element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), #panel.border = element_blank())
        axis.ticks.length = unit(0, "mm")) + expand_limits(x = 0, y = 0)





#### heatmap again for presentation, mapping locations

cols <- viridis::viridis(120)
cols[121:1000] <- "#FFFFFF"
cols <- base::rev(cols)

g2.names$end.var1 <- as.numeric(g2.names$end.var1)
g2.names$start.var1 <- as.numeric(g2.names$start.var1)
g2.names$mid.var1 <- ((g2.names$end.var1+g2.names$start.var1) / 2)

MyHeatmap <- LDheatmap(g2.iso.pair.plot, genetic.distances = g2.names$mid.var1, color=cols, add.key=TRUE, flip=TRUE)


spacer1 <- rep(1, nrow(g2.iso.pair.plot))
g2.iso.pair.plot.space <- base::rbind(spacer1,g2.iso.pair.plot)
row.names(g2.iso.pair.plot.space)[1] <- "spacer1"
spacer1 <- rep(1, nrow(g2.iso.pair.plot.space))
g2.iso.pair.plot.space <- base::cbind(spacer1,g2.iso.pair.plot.space)
spacer2 <- rep(1, nrow(g2.iso.pair.plot.space))
g2.iso.pair.plot.space <- base::rbind(g2.iso.pair.plot.space, spacer2)
row.names(g2.iso.pair.plot.space)[nrow(g2.iso.pair.plot.space)] <- "spacer2"
spacer2 <- rep(1, nrow(g2.iso.pair.plot.space))
g2.iso.pair.plot.space <- base::cbind(g2.iso.pair.plot.space, spacer2)

g2.iso.pair.plot.space <- as.matrix(g2.iso.pair.plot.space)

MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names$mid.var1,139957), color=cols, add.key=TRUE, flip=TRUE)


#pdf(Figure_Ldheatmap_anno.pdf",width=10,height=10)
llplusgrob<-LDheatmap.addGrob(MyHeatmap, as.grob(anno.plot), height = 0.3)
#dev.off()




#### For the whole Y Chromosome


#### convert original g2 matrix for Y chromosome into actual matrix

g2.iso.pair <- base::as.matrix(g2.mat)

orig.col <- ncol(g2.iso.pair)

#### add sizes of each contig for conversion

sizes <- read.csv("dmel_scaffold2_plus0310_sizes.csv", header=FALSE)
colnames(sizes) <- c("contig","len")

#### convert contig position to chromosomal position

chrom.pos.convert <- function(df, sizes) { 
  
  ## make sure the variables are numeric
  df$start.var1 <- as.numeric(df$start.var1)
  df$end.var1 <- as.numeric(df$end.var1)
  
  ## organize contigs
  reverse_contigs=c('2L_2','2R_contig_19','3L_4','3R_5','3_scaffold1','Contig5','Y_scaffold7','Y_Contig10','Y_Contig104')
  chrX <- c('X_1','Contig5','X_2','X_3','X_7','X_9','Contig135','X_10','X3X4_6_D1712_2','Contig95','Contig79')
  chr2 <- c('2L_1','2L_2','tig00057289','Contig142','Rsp_2','2R_16','2R_18','2R_19','2R_20','2R_21')
  chr3 <- c('3L_1','3L_3','3L_4','3L_5','3L_8','3L_10','3R_1','3R_2','3R_3','3R_4','3R_5','tig00022795','id=67822_0','3R_6','3R_8','3R_9','3R_10','Contig11','Contig145','3_scaffold1','3_scaffold2','3R_28')
  chr4 <- c('Contig119','4_2')
  chrY <- c('Y_scaffold6','Y_scaffold7','Y_scaffold4','Y_Contig140','Y_Contig143','Y_Contig10','Y_Contig104','Y_Contig6','Y_Contig2','Y_Contig26','Y_scaffold5','Y_scaffold3')
  
  chroms <- list(chrX,chr2,chr3,chr4,chrY)
  names(chroms) <- c("X","2","3","4","Y")
  
  list.out <- list(chromX=NULL,chrom2=NULL,chrom3=NULL,chrom4=NULL,chromY=NULL)
  
  ## empty data frame for replacement if need be
  #emp.df <- data.frame("te.name" = NA, "species.var1" = NA,"contig.var1" = NA,"start.var1" = NA,"end.var1" = NA,"direction.var1" = NA)
  
  for (i in seq_along(chroms)) {
    if (nrow(df[df$contig.var1 %in% chroms[[i]],]) > 0) {
      list.out[[i]] <- df[df$contig.var1 %in% chroms[[i]],]
      list.out[[i]]$chrom <- names(chroms[i])
      list.out[[i]]$chrom.start.var1 <- NA
      list.out[[i]]$chrom.end.var1 <- NA
    } else {
      chroms[[i]] <- NA
    }
  }
  
  ## remove empty parts of that list
  list.out <- list.out[lapply(list.out, length) > 0]
  
  chroms <- chroms[lapply(chroms, length) > 1]
  
  #ref.by.chrom <- list.out
  for (k in 1:length(chroms)) { 
    for (i in 1:length(list.out[[k]]$contig.var1)) {
      pos <- 0
      for (j in 1:length(chroms[[k]])) {
        if (list.out[[k]]$contig.var1[i]==chroms[[k]][j]) {
          list.out[[k]]$chrom
          if (chroms[[k]][j] %in% reverse_contigs) {
            list.out[[k]]$chrom.start.var1[i] <- sizes$len[sizes$contig==chroms[[k]][j]] - list.out[[k]]$start.var1[i] + pos
            list.out[[k]]$chrom.end.var1[i] <- sizes$len[sizes$contig==chroms[[k]][j]] - list.out[[k]]$end.var1[i] + pos
          } else {
            list.out[[k]]$chrom.start.var1[i] <- list.out[[k]]$start.var1[i] + pos
            list.out[[k]]$chrom.end.var1[i] <- list.out[[k]]$end.var1[i] + pos
          }
        }
        pos <- pos+sizes$len[sizes$contig==chroms[[k]][j]]+100000
      }
    }
  }
  ## convert to dataframe to export
  
  df.out <- do.call(rbind, list.out)
  rownames(df.out) <- 1:nrow(df.out)
  
  return(df.out)
}

#### add spacer

add.spacer <- function(df) {
  
  ## adding space of NAs around
  spacer1 <- rep(NA, nrow(df))
  df.out <- base::rbind(spacer1,df)
  spacer1 <- rep(NA, nrow(df.out))
  df.out <- base::cbind(spacer1,df.out)
  spacer2 <- rep(NA, nrow(df.out))
  df.out <- base::rbind(df.out, spacer2)
  spacer2 <- rep(NA, nrow(df.out))
  df.out <- base::cbind(df.out, spacer2)
  
  return(as.matrix(df.out))
}


## info about the chromosomes for illustration

reverse_contigs=c('2L_2','2R_contig_19','3L_4','3R_5','3_scaffold1','Contig5','Y_scaffold7','Y_Contig10','Y_Contig104')
chrX <- c('X_1','Contig5','X_2','X_3','X_7','X_9','Contig135','X_10','X3X4_6_D1712_2','Contig95','Contig79')
chr2 <- c('2L_1','2L_2','tig00057289','Contig142','Rsp_2','2R_16','2R_18','2R_19','2R_20','2R_21')
chr3 <- c('3L_1','3L_3','3L_4','3L_5','3L_8','3L_10','3R_1','3R_2','3R_3','3R_4','3R_5','tig00022795','id=67822_0','3R_6','3R_8','3R_9','3R_10','Contig11','Contig145','3_scaffold1','3_scaffold2','3R_28')
chr4 <- c('Contig119','4_2')
chrY <- c('Y_scaffold6','Y_scaffold7','Y_scaffold4','Y_Contig140','Y_Contig143','Y_Contig10','Y_Contig104','Y_Contig6','Y_Contig2','Y_Contig26','Y_scaffold5','Y_scaffold3')

sizeX <- sum(sizes$len[sizes$contig %in% chrX]) + ((length(chrX)-1)*100000)
size2 <- sum(sizes$len[sizes$contig %in% chr2]) + ((length(chr2)-1)*100000)
size3 <- sum(sizes$len[sizes$contig %in% chr3]) + ((length(chr3)-1)*100000)
size4 <- sum(sizes$len[sizes$contig %in% chr4]) + ((length(chr4)-1)*100000)
sizeY <- sum(sizes$len[sizes$contig %in% chrY]) + ((length(chrY)-1)*100000)




## data table for whole Y chrom

g2.iso.pair <- base::as.matrix(g2.mat)

orig.col <- ncol(g2.iso.pair)

g2.names <- as.data.frame(base::colnames(g2.iso.pair))
colnames(g2.names) <- c("te.name")

g2.names <- g2.names %>% separate(col = te.name, into = c("species.var1","contig.var1","span","direction.var1"), sep = "\\|", remove=F) %>% 
  separate(col = span, into=c("start.var1","end.var1"), sep = "-",remove = T) %>%
  mutate(start.var1=as.numeric(start.var1)) 


g2.names.chrom <- chrom.pos.convert(g2.names, sizes)
g2.names.chrom <- g2.names.chrom %>% arrange(chrom.start.var1)
g2.iso.pair.plot <- g2.iso.pair[match(g2.names.chrom$te.name,row.names(g2.iso.pair)),match(g2.names.chrom$te.name,colnames(g2.iso.pair))]
g2.iso.pair.plot.space <- add.spacer(g2.iso.pair.plot)


g2.names.chrom$end.var1 <- as.numeric(g2.names.chrom$end.var1)
g2.names.chrom$start.var1 <- as.numeric(g2.names.chrom$start.var1)


## Figue for whole Y chromosome

cols <- viridis::viridis(120)
cols[121:1000] <- "#FFFFFF"
cols <- base::rev(cols)

## for highlighting centromere
cent.copy.index <- which(grepl("Y_Contig26", row.names(g2.iso.pair.plot.space)))

#pdf(Ldheatmap_wholeY.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names.chrom$chrom.start.var1,sizeY), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1], cent.copy.index[length(cent.copy.index)], fill="NA", lwd=2)
#dev.off()

#### lines for marking contig divisions

df <- data.frame()

sizesY <- sizes %>% filter(contig %in% chrY)
sizesY <- sizesY[match(chrY, sizesY$contig),]

#pdf(Ldheatmap_anno_ContigLinesY.pdf",width=10,height=2)
ggplot(df) + geom_point() +
  geom_rect(data=NULL, aes(xmin=0,xmax=100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sizesY$len[1],xmax=sizesY$len[1]+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:2])+(100000*1),xmax=sum(sizesY$len[1:2])+(100000*1)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:3])+(100000*2),xmax=sum(sizesY$len[1:3])+(100000*2)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:4])+(100000*3),xmax=sum(sizesY$len[1:4])+(100000*3)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:5])+(100000*4),xmax=sum(sizesY$len[1:5])+(100000*4)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:6])+(100000*5),xmax=sum(sizesY$len[1:6])+(100000*5)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:7])+(100000*6),xmax=sum(sizesY$len[1:7])+(100000*6)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:8])+(100000*7),xmax=sum(sizesY$len[1:8])+(100000*7)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:9])+(100000*8),xmax=sum(sizesY$len[1:9])+(100000*8)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:10])+(100000*9),xmax=sum(sizesY$len[1:10])+(100000*9)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:11])+(100000*10),xmax=sum(sizesY$len[1:11])+(100000*10)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesY$len[1:12])+(100000*10),xmax=sum(sizesY$len[1:12])+(100000*10)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  theme_bw() + #xlim(0, sizeY) +
  ylim(0,1) + 
  scale_x_continuous(breaks=c(0,50000000,sizeY), expand=c(0,0)) +
  #scale_y_continuous(breaks=c(0,1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#dev.off()




#### #### #### 
#### X Chromosome 

## import data
g2.mat1 <- read.csv("X_chrom_pairwise_dist1.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)
g2.mat2 <- read.csv("X_chrom_pairwise_dist2.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)


g2.iso.pair <- as.matrix(combine.matrix(g2.mat1, g2.mat2))

## dataframe for ordering matrix

g2.names <- as.data.frame(base::colnames(g2.iso.pair))
colnames(g2.names) <- c("te.name")

## adding names and position information
g2.names <- g2.names %>% separate(col = te.name, into = c("species.var1","contig.var1","span","direction.var1"), sep = "\\|", remove=F) %>% 
  separate(col = span, into=c("start.var1","end.var1"), sep = "-",remove = T) %>%
  mutate(start.var1=as.numeric(start.var1)) 

## convert coordinates and arrange by new coordinates
g2.names.chrom <- chrom.pos.convert(g2.names, sizes)
g2.names.chrom <- g2.names.chrom %>% arrange(chrom.start.var1)

## now match them up
g2.iso.pair.plot <- g2.iso.pair[match(g2.names.chrom$te.name,row.names(g2.iso.pair)),match(g2.names.chrom$te.name,colnames(g2.iso.pair))]

## add spacer to get full length of chromosome
g2.iso.pair.plot.space <- add.spacer(g2.iso.pair.plot)

## add colors and plot
cols <- viridis::viridis(120)
cols[121:1000] <- "#FFFFFF"
cols <- base::rev(cols)

cent.copy.index <- which(grepl("Contig79", row.names(g2.iso.pair.plot.space)))

#pdf("Ldheatmap_wholeX.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names.chrom$chrom.start.var1,sizeX), color=cols, add.key=TRUE, flip=TRUE, pop=F)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1], cent.copy.index[length(cent.copy.index)], fill="NA", lwd=2)
#dev.off()

## for X add one without spacer since most of the X chromosome doesn't have Jockey elements
#pdf("Ldheatmap_wholeX_zoomed.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot, genetic.distances = c(g2.names.chrom$chrom.start.var1), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1]-1, cent.copy.index[length(cent.copy.index)-1], fill="NA", lwd=2)
#dev.off()


#### lines for marking contig divisions

df <- data.frame()

sizesX <- sizes %>% filter(contig %in% chrX)
sizesX <- sizesX[match(chrX, sizesX$contig),]


## Whole X Chrom
#pdf("Ldheatmap_wholeX_ContigLines.pdf",width=10,height=2)
ggplot(df) + geom_point() +
  geom_rect(data=NULL, aes(xmin=0,xmax=100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sizesX$len[1],xmax=sizesX$len[1]+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:2])+(100000*1),xmax=sum(sizesX$len[1:2])+(100000*1)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:3])+(100000*2),xmax=sum(sizesX$len[1:3])+(100000*2)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:4])+(100000*3),xmax=sum(sizesX$len[1:4])+(100000*3)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:5])+(100000*4),xmax=sum(sizesX$len[1:5])+(100000*4)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:6])+(100000*5),xmax=sum(sizesX$len[1:6])+(100000*5)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:7])+(100000*6),xmax=sum(sizesX$len[1:7])+(100000*6)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:8])+(100000*7),xmax=sum(sizesX$len[1:8])+(100000*7)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:9])+(100000*8),xmax=sum(sizesX$len[1:9])+(100000*8)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:10])+(100000*9),xmax=sum(sizesX$len[1:10])+(100000*9)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:11])+(100000*9),xmax=sum(sizesX$len[1:11])+(100000*9)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  theme_bw() + #xlim(0, sizeX) + 
  ylim(0,1) + scale_x_continuous(breaks=c(0,50000000,sizeX), expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#dev.off()

## Subset of X Chrom
#pdf("Ldheatmap_wholeX_zoomed_ContigLines.pdf",width=10,height=2)
ggplot(df) + geom_point() +
  geom_rect(data=NULL, aes(xmin=sizesX$len[2],xmax=sizesX$len[2]+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  #geom_rect(data=NULL, aes(xmin=sizesX$len[2],xmax=sizesX$len[2]+10000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:3])+(100000*1),xmax=sum(sizesX$len[2:3])+(100000*1)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:4])+(100000*2),xmax=sum(sizesX$len[2:4])+(100000*2)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:5])+(100000*3),xmax=sum(sizesX$len[2:5])+(100000*3)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:6])+(100000*4),xmax=sum(sizesX$len[2:6])+(100000*4)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:7])+(100000*5),xmax=sum(sizesX$len[2:7])+(100000*5)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:8])+(100000*6),xmax=sum(sizesX$len[2:8])+(100000*6)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:9])+(100000*7),xmax=sum(sizesX$len[2:9])+(100000*7)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:10])+(100000*8),xmax=sum(sizesX$len[2:10])+(100000*8)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizesX$len[2:11])+(100000*8),xmax=sum(sizesX$len[2:11])+(100000*8)+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  #geom_rect(data=NULL, aes(xmin=sum(sizesX$len[1:11])+(100000*10),xmax=sum(sizesX$len[1:11])+(100000*10)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  theme_bw() + #xlim(55221, sizeX - (sizesX$len[1]+100000)) + 
  ylim(0,1) + scale_x_continuous(breaks=c(55221,sizeX - (sizesX$len[1]+100000)), expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#dev.off()





#### #### 
#### 2nd chrom

## import data
g2.mat1 <- read.csv("2_chrom_pairwise_dist1.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)
g2.mat2 <- read.csv("2_chrom_pairwise_dist2.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)


## combine matrices

g2.iso.pair <- as.matrix(combine.matrix(g2.mat1, g2.mat2))

## dataframe for ordering matrix

g2.names <- as.data.frame(base::colnames(g2.iso.pair))
colnames(g2.names) <- c("te.name")

## adding names and position information
g2.names <- g2.names %>% separate(col = te.name, into = c("species.var1","contig.var1","span","direction.var1"), sep = "\\|", remove=F) %>% 
  separate(col = span, into=c("start.var1","end.var1"), sep = "-",remove = T) %>%
  mutate(start.var1=as.numeric(start.var1)) 

## convert coordinates and arrange by new coordinates
g2.names.chrom <- chrom.pos.convert(g2.names, sizes)
g2.names.chrom <- g2.names.chrom %>% arrange(chrom.start.var1)

## now match them up
g2.iso.pair.plot <- g2.iso.pair[match(g2.names.chrom$te.name,row.names(g2.iso.pair)),match(g2.names.chrom$te.name,colnames(g2.iso.pair))]

## add spacer to get full length of chromosome
g2.iso.pair.plot.space <- add.spacer(g2.iso.pair.plot)

## add colors and plot
cols <- viridis::viridis(120)
cols[121:1000] <- "#FFFFFF"
cols <- base::rev(cols)

cent.copy.index <- which(grepl("tig00057289", row.names(g2.iso.pair.plot.space)))

#pdf("Ldheatmap_whole2.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names.chrom$chrom.start.var1,size2), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1], cent.copy.index[length(cent.copy.index)], fill="NA", lwd=2)
#dev.off()

#### lines for marking contig divisions

df <- data.frame()
sizes2 <- sizes %>% filter(contig %in% chr2)
sizes2 <- sizes2[match(chr2, sizes2$contig),]

#pdf("Ldheatmap_whole2_ContigLines.pdf",width=10,height=2)
ggplot(df) + geom_point() +
  geom_rect(data=NULL, aes(xmin=0,xmax=100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sizes2$len[1],xmax=sizes2$len[1]+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:2])+(100000*1),xmax=sum(sizes2$len[1:2])+(100000*1)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:3])+(100000*2),xmax=sum(sizes2$len[1:3])+(100000*2)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:4])+(100000*3),xmax=sum(sizes2$len[1:4])+(100000*3)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:5])+(100000*4),xmax=sum(sizes2$len[1:5])+(100000*4)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:6])+(100000*5),xmax=sum(sizes2$len[1:6])+(100000*5)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:7])+(100000*6),xmax=sum(sizes2$len[1:7])+(100000*6)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:8])+(100000*7),xmax=sum(sizes2$len[1:8])+(100000*7)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:9])+(100000*8),xmax=sum(sizes2$len[1:9])+(100000*8)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:10])+(100000*8),xmax=sum(sizes2$len[1:10])+(100000*8)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  #geom_rect(data=NULL, aes(xmin=sum(sizes2$len[1:11])+(100000*10),xmax=sum(sizes2$len[1:11])+(100000*10)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  theme_bw() + #xlim(0, size2) + 
  ylim(0,1) + scale_x_continuous(breaks=c(0,50000000,size2), expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#dev.off()


#### #### 
#### 3rd chrom

## import data
pwd <- "/Users/lucashemmer/Documents/mel_centromere/ReAnnotation/phylogeny/Dmel_Ref_Denovo_aligntoRef_Y_chrom/"
g2.mat1 <- read.csv("3_chrom_pairwise_dist1.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)
g2.mat2 <- read.csv("3_chrom_pairwise_dist2.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)


## combine matrices

g2.iso.pair <- as.matrix(combine.matrix(g2.mat1, g2.mat2))

## dataframe for ordering matrix

g2.names <- as.data.frame(base::colnames(g2.iso.pair))
colnames(g2.names) <- c("te.name")

## adding names and position information
g2.names <- g2.names %>% separate(col = te.name, into = c("species.var1","contig.var1","span","direction.var1"), sep = "\\|", remove=F) %>% 
  separate(col = span, into=c("start.var1","end.var1"), sep = "-",remove = T) %>%
  mutate(start.var1=as.numeric(start.var1)) 

## convert coordinates and arrange by new coordinates
g2.names.chrom <- chrom.pos.convert(g2.names, sizes)
g2.names.chrom <- g2.names.chrom %>% arrange(chrom.start.var1)

## now match them up
g2.iso.pair.plot <- g2.iso.pair[match(g2.names.chrom$te.name,row.names(g2.iso.pair)),match(g2.names.chrom$te.name,colnames(g2.iso.pair))]

## add spacer to get full length of chromosome
g2.iso.pair.plot.space <- add.spacer(g2.iso.pair.plot)

## add colors and plot
cols <- viridis::viridis(120)
cols[121:1000] <- "#FFFFFF"
cols <- base::rev(cols)

cent.copy.index <- which(grepl("3R_5", row.names(g2.iso.pair.plot.space)))

#pdf("Ldheatmap_whole3.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names.chrom$chrom.start.var1,size3), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1], cent.copy.index[length(cent.copy.index)], fill="NA", lwd=2)
#dev.off()

#### lines for marking contig divisions

df <- data.frame()

sizes3 <- sizes %>% filter(contig %in% chr3)
sizes3 <- sizes3[match(chr3, sizes3$contig),]

#pdf("Ldheatmap_whole3_ContigLines.pdf",width=10,height=2)
ggplot(df) + geom_point() +
  geom_rect(data=NULL, aes(xmin=0,xmax=100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sizes3$len[1],xmax=sizes3$len[1]+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:2])+(100000*1),xmax=sum(sizes3$len[1:2])+(100000*1)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:3])+(100000*2),xmax=sum(sizes3$len[1:3])+(100000*2)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:4])+(100000*3),xmax=sum(sizes3$len[1:4])+(100000*3)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:5])+(100000*4),xmax=sum(sizes3$len[1:5])+(100000*4)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:6])+(100000*5),xmax=sum(sizes3$len[1:6])+(100000*5)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:7])+(100000*6),xmax=sum(sizes3$len[1:7])+(100000*6)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:8])+(100000*7),xmax=sum(sizes3$len[1:8])+(100000*7)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:9])+(100000*8),xmax=sum(sizes3$len[1:9])+(100000*8)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:10])+(100000*9),xmax=sum(sizes3$len[1:10])+(100000*9)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:11])+(100000*10),xmax=sum(sizes3$len[1:11])+(100000*10)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:12])+(100000*11),xmax=sum(sizes3$len[1:12])+(100000*11)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:13])+(100000*12),xmax=sum(sizes3$len[1:13])+(100000*12)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:14])+(100000*13),xmax=sum(sizes3$len[1:14])+(100000*13)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:15])+(100000*14),xmax=sum(sizes3$len[1:15])+(100000*14)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:16])+(100000*15),xmax=sum(sizes3$len[1:16])+(100000*15)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:17])+(100000*16),xmax=sum(sizes3$len[1:17])+(100000*16)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:18])+(100000*17),xmax=sum(sizes3$len[1:18])+(100000*17)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:19])+(100000*18),xmax=sum(sizes3$len[1:19])+(100000*18)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:20])+(100000*19),xmax=sum(sizes3$len[1:20])+(100000*19)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:21])+(100000*20),xmax=sum(sizes3$len[1:21])+(100000*20)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes3$len[1:22])+(100000*20),xmax=sum(sizes3$len[1:22])+(100000*20)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  theme_bw() + #xlim(0, size3) + 
  ylim(0,1) + scale_x_continuous(breaks=c(0,size3), expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#dev.off()


#### #### 
#### 4th chrom

## import data
pwd <- "/Users/lucashemmer/Documents/mel_centromere/ReAnnotation/phylogeny/Dmel_Ref_Denovo_aligntoRef_Y_chrom/"
g2.mat1 <- read.csv("4_chrom_pairwise_dist1.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)
g2.mat2 <- read.csv("4_chrom_pairwise_dist2.csv", header=TRUE,stringsAsFactors = F, fill = F, check.names=F,row.names = 1)


## combine matrices

g2.iso.pair <- as.matrix(combine.matrix(g2.mat1, g2.mat2))

## dataframe for ordering matrix

g2.names <- as.data.frame(base::colnames(g2.iso.pair))
colnames(g2.names) <- c("te.name")

## adding names and position information
g2.names <- g2.names %>% separate(col = te.name, into = c("species.var1","contig.var1","span","direction.var1"), sep = "\\|", remove=F) %>% 
  separate(col = span, into=c("start.var1","end.var1"), sep = "-",remove = T) %>%
  mutate(start.var1=as.numeric(start.var1)) 

## convert coordinates and arrange by new coordinates
g2.names.chrom <- chrom.pos.convert(g2.names, sizes)
g2.names.chrom <- g2.names.chrom %>% arrange(chrom.start.var1)

## now match them up
g2.iso.pair.plot <- g2.iso.pair[match(g2.names.chrom$te.name,row.names(g2.iso.pair)),match(g2.names.chrom$te.name,colnames(g2.iso.pair))]

## add spacer to get full length of chromosome
g2.iso.pair.plot.space <- add.spacer(g2.iso.pair.plot)

## add colors and plot
cols <- viridis::viridis(120)
cols[121:1000] <- "#FFFFFF"
cols <- base::rev(cols)

cent.copy.index <- which(grepl("Contig119", row.names(g2.iso.pair.plot.space)))

#pdf("Ldheatmap_whole4.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names.chrom$chrom.start.var1,size4), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1], cent.copy.index[length(cent.copy.index)], fill="NA", lwd=2)
#dev.off()

#pdf("Ldheatmap_whole4_zoomed.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot.space, genetic.distances = c(0,g2.names.chrom$chrom.start.var1,93914), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1], cent.copy.index[length(cent.copy.index)], fill="NA", lwd=2)
#dev.off()

#pdf("Ldheatmap_whole4_zoomed1.pdf",width=10,height=10)
MyHeatmap <- LDheatmap(g2.iso.pair.plot, genetic.distances = c(g2.names.chrom$chrom.start.var1), color=cols, add.key=TRUE, flip=TRUE)
LDheatmap.highlight(MyHeatmap, i=cent.copy.index[1]-1, cent.copy.index[length(cent.copy.index)-1], fill="NA", lwd=2)
#dev.off()



df <- data.frame()

sizes4 <- sizes %>% filter(contig %in% chr4)
sizes4 <- sizes4[match(chr4, sizes4$contig),]

#pdf("Ldheatmap_whole4_ContigLines.pdf",width=10,height=2)
ggplot(df) + geom_point() +
  geom_rect(data=NULL, aes(xmin=0,xmax=50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sizes4$len[1],xmax=sizes4$len[1]+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  geom_rect(data=NULL, aes(xmin=sum(sizes4$len[1:2]),xmax=sum(sizes4$len[1:2])+50000,ymin=-Inf,ymax=Inf),fill="grey87") +
  #geom_rect(data=NULL, aes(xmin=sum(sizes4$len[1:2])+(100000*1),xmax=sum(sizes4$len[1:2])+(100000*1)+100000,ymin=-Inf,ymax=Inf),fill="grey87") +
  theme_bw() + #xlim(0, size4) + 
  ylim(0,1) + scale_x_continuous(breaks=c(0,size4), expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
#dev.off()




