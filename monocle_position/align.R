library(Seurat)
library(ggplot2)
library(cowplot)
library(monocle3)
library(dplyr)
library(readr)

datasets <- readRDS("../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")

tissues <- c("Meristem_2","Trasition_zone","Cortex","Endodermis")


phases <- c()
lengths <- c()
range <- data.frame(name=c(),meristem=c(),transition=c())
n=1
for(name in names(datasets@images))
    {
    name2 <- gsub("image_BIN40","",name)
    cor <- read.csv(paste0("../geneswitches/files/",name2,".csv"),header=T,row.names=1)
    datasets@images[[name]]@coordinates[,c('x','y')] <- cor
    temp <- subset(x=datasets,idents=tissues)
    tip <- temp@images[[name]]@coordinates[,'y'] %>% min()
    phase <- (datasets@images[[name]]@coordinates[,'y'] -tip)%/% 10
    length <- (datasets@images[[name]]@coordinates[,'y'] -tip)
    names(phase) <- rownames(datasets@images[[name]]@coordinates)
    names(length) <- rownames(datasets@images[[name]]@coordinates)
    phases <- c(phases,phase)
    lengths <-c(lengths,length)
    temp1 <- subset(x=datasets,idents=tissues[-1])
    transition_start <- temp1@images[[name]]@coordinates[,'y'] %>% min()
    temp2 <- subset(x=datasets,idents=tissues[-c(1,2)])
    cortex_start <- temp2@images[[name]]@coordinates[,'y'] %>% min()
    range[n,'name'] <- name2
    range[n,'meristem'] <- transition_start-tip
    range[n,'transitioin'] <- cortex_start-tip 
    n=n+1
    }

write.csv(file="range.csv",range,quote=F,row.names=F)
save.image(file="align.RData")
#datasets <- AddMetaData(object = datasets, metadata = phases,col.name = 'phase')
#datasets <- AddMetaData(object = datasets, metadata = lengths,col.name = 'length')
