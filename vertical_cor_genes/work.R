library(Seurat)
library(ggplot2)
library(cowplot)
library(monocle3)
library(dplyr)
library(readr)

datasets <- readRDS("../../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")

tissues <- c("Meristem_2","Endodermis","Cortex","Trasition_zone","Meristem_1","Stele_1","Stele_2","Stele_3")

load("../../monocle_position/align.RData")
datasets <- AddMetaData(object = datasets, metadata = phases,col.name = 'phase')
datasets <- AddMetaData(object = datasets, metadata = lengths,col.name = 'length')

## choose the interested clusters
datasetso <- datasets
datasets<-subset(x=datasets,idents=tissues)
datasets$type <- "cortex"
datasets@meta.data[datasets@meta.data$cell.type %in% c("Meristem_1","Stele_1","Stele_2","Stele_3"),'type'] <- "stele"

rdsf <- read_tsv("../../info2")
temp<-left_join(datasets@meta.data,rdsf,by="group")
datasets@meta.data$species<-temp$species
datasets@meta.data$medium<-temp$medium
#dataset <- subset(x=datasets, medium=="Y")
datasets.list <- SplitObject(datasets, split.by="type")

dev <- read.csv("02_only_genes.csv",header=T)
dev$gene <- dev$id

find_medium_diff <- function(datasets,gene)
    {
    temp <- datasets.list[[1]][['SCT']][gene,] %>% as.data.frame
    cells <- colnames(temp[,temp>0])
    length1 <- datasets$length[cells]
    temp <- datasets.list[[2]][['SCT']][gene,] %>% as.data.frame
    cells <- colnames(temp[,temp>0])
    length2 <- datasets$length[cells]
    diff <- abs(median(length1)-median(length2))
    return(diff)
    }

find_length <- function(datasets,gene,sp)
    {
    temp <- datasets.list[[1]][['SCT']][gene,] %>% as.data.frame
    cells <- colnames(temp[,temp>0])
    length1 <- datasets$length[cells]
    temp <- datasets.list[[2]][['SCT']][gene,] %>% as.data.frame
    cells <- colnames(temp[,temp>0])
    length2 <- datasets$length[cells]
    result1 <- try(var.test(length1,length2,alternative = "two.sided"),silent=TRUE)
    if ('try-error' %in% class(result1))
        {final <- FALSE}else{
        if(result1$p.value >0.05)
            {result <- t.test(length1,length2,alternative="two.sided",paired=F,var.equal=T,conf.level=0.95)
            }else{
             result <- t.test(length1,length2,alternative="two.sided",paired=F,var.equal=F,conf.level=0.95)}
        if(result$p.value > sp)
             {final <- TRUE}else{
             final <- FALSE}}
    return(final)
    }


a <- c()
dev$length_diff <- "F"
for(gene in dev$gene)
   {if(gene %in% rownames(datasets[['SCT']]))
    {final <- find_length(datasets,gene,sp=0.05)
     diff <- find_medium_diff(datasets,gene)
     a <- c(a,final)
     genename <- gene
     dev[dev$gene==genename,]$length_diff <- diff
    }}

dev2 <- dev[a,]
write.csv(file="02_core_dev_gene2_difflength.csv",dev,quote=FALSE,row.names=FALSE)
write.csv(file="02_vertical_core_dev_genes.csv",dev2,quote=FALSE,row.names=FALSE)
