library(Seurat)
library(viridis)
library(monocle3)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
datasets <- readRDS("../../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")
tissues <- c("Meristem_2","Endodermis","Cortex","Trasition_zone")
phases <- c()
lengths <- c()
for(name in names(datasets@images))
    {
    name2 <- gsub("image_BIN40","",name)
    cor <- read.csv(paste0("../../geneswitches/files/",name2,".csv"),header=T,row.names=1)
    datasets@images[[name]]@coordinates[,c('x','y')] <- cor
    temp <- subset(x=datasets,idents=tissues)
    tip <- temp@images[[name]]@coordinates[,'y'] %>% min()
    phase <- (datasets@images[[name]]@coordinates[,'y'] -tip)%/% 10 
    length <- (datasets@images[[name]]@coordinates[,'y'] -tip)
    names(phase) <- rownames(datasets@images[[name]]@coordinates)
    names(length) <- rownames(datasets@images[[name]]@coordinates)
    phases <- c(phases,phase)
    lengths <-c(lengths,length)
    }   

datasets <- AddMetaData(object = datasets, metadata = phases,col.name = 'phase')
datasets <- AddMetaData(object = datasets, metadata = lengths,col.name = 'length')

datasets.list <- SplitObject(datasets, split.by="group")

module_info <- read.csv("phase_module.csv",row.names=1)
gene_group_df<-module_info[,c('id','module','gene_name')]

dir.create("images")
setwd("images")
meta <- list()
for (i in 1:length(datasets.list))
    {
    BRAIN <- subset(x=datasets.list[[i]],idents=tissues)
    data <- as(as.matrix(BRAIN@assays$Spatial@counts), 'sparseMatrix')
    pd <- data.frame(BRAIN@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

    cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)

    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$phase)
    agg_mat <- aggregate_gene_expression(cds, gene_group_df, cell_group_df)
    names <- colnames(agg_mat)
    agg_mat <- as.data.frame(agg_mat)
    colnames(agg_mat) <- names
    agg_mat$modules <- paste0("module_",rownames(agg_mat))

    meta[[i]] <- agg_mat %>% pivot_longer(!modules, names_to = "phase", values_to = "avg_count")
    }

for (j in unique(module_info$module))
    {
     j <- paste0("module_",j)
     all_avg <- c()
    for (i in 1:length(datasets.list))
        {
        temp <- meta[[i]][meta[[i]]$modules==j,]$avg_count
        all_avg <- c(all_avg, temp)
        }
     p <- list()
     for (i in 1:length(datasets.list))
        {
        BRAIN <- subset(x=datasets.list[[i]],idents=tissues)   
        temp0 <- left_join(BRAIN@meta.data,meta[[i]][meta[[i]]$modules==j,c('phase', 'avg_count')]%>% mutate(phase = as.numeric(phase)),by='phase')
        temp0$cell_id <- rownames(BRAIN@meta.data)
        imagename <- paste0("image_",datasets.list[[i]]$group[1])
        imagename <- gsub("-",".",imagename)
        temp1 <- datasets@images[[imagename]]@coordinates[,1:2]
        temp1$cell_id <- rownames(temp1)
        temp<-left_join(temp1,temp0,by = 'cell_id')
        p[[i]]<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = avg_count),size = 1.5)+
        scale_color_viridis(option="plasma",
                            end=max(temp$avg_count[!is.na(temp$avg_count)])/max(all_avg),na.value="grey")+
        ggtitle(imagename)+
        theme(panel.grid=element_blank())
        }
    sum<-grid.arrange(grobs=p,ncol=6)
    ggsave2(paste0(j,".pdf"),sum,dpi=600,width=13,height=15)
    }


    



