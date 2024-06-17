library(getopt)

arg <- matrix(c("input", "i","1","character","input file1",
                "info","f","1","character","information file",
		"logfc","l","1","numeric","log fold change to filter the markers"
                ),byrow=T,ncol=5)

opt = getopt(arg)

library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(viridis)
library(cowplot)
library(gridExtra)
if (FALSE)
{
opt <- list()
opt$input <- "data_pos_adj.RDS"
opt$info <- "../info2"
opt$logfc <- 0.05
}

BRAIN<-readRDS(opt$input)
rdsf <- read.table(opt$info,sep="\t",header=T)

genename<-read.table("../../Markergenes/gene_v2.txt",sep="\t",quote="",header=T)
colnames(genename) <- c("gene","name")
genename <- unique(genename)

##pretreatment
temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$chip_name<-temp$chip_name
BRAIN@meta.data$line<-temp$line
BRAIN@meta.data$species<-temp$species
BRAIN@meta.data$cell.type<-Idents(BRAIN)

DefaultAssay(BRAIN) <- "SCT"


single_neighbor <- function(datasets,cell_id)
    {cor <- strsplit(cell_id,"_")
     X <- gsub("X","",cor[[1]][2]) %>% as.numeric()
     Y <- cor[[1]][3] %>% as.numeric()
     neighbors <- c(paste0(cor[[1]][1],"_X",X-1,"_",Y),
                   paste0(cor[[1]][1],"_X",X+1,"_",Y),
                   paste0(cor[[1]][1],"_X",X,"_",Y-1),
                   paste0(cor[[1]][1],"_X",X,"_",Y+1))
     surroundings <- c(paste0(cor[[1]][1],"_X",X-1,"_",Y+1),
                   paste0(cor[[1]][1],"_X",X+1,"_",Y+1),
                   paste0(cor[[1]][1],"_X",X-1,"_",Y-1),
                   paste0(cor[[1]][1],"_X",X+1,"_",Y-1))
     #removve non_existed spot
     neighbors <- neighbors[!is.na(Idents(datasets)[neighbors])]
     surroundings <- surroundings[!is.na(Idents(datasets)[surroundings])]
     return(c(neighbors,surroundings))    
    }


find_neighbor <- function(datasets,tissues,iter)
    {t1 <- c()
     t2 <- c()
    for(cell_id in WhichCells(datasets,ident=tissues[1]))
        {neighbors <- single_neighbor(datasets,cell_id)
         dis <- Idents(datasets)[neighbors] %>% as.character()
         if(tissues[2] %in% dis)
              {t1 <- c(t1,cell_id)
              for (cell_id2 in neighbors)
                   {temp <- neighbors[Idents(datasets)[neighbors]==tissues[2]]
                    t2 <- c(t2,temp)}}}
     t2 <- unique(t2)
     t0 <- c(t1,t2)
     i=1
     while(i < iter)
         {t11 <- c()
          t22 <- c()
          for(cell_id in t1)
              {neighbors <- single_neighbor(datasets,cell_id)
               temp <- neighbors[Idents(datasets)[neighbors]==tissues[1]]
               t11 <-c(t11,temp)}
          for(cell_id in t2)
              {neighbors <- single_neighbor(datasets,cell_id)
               temp <- neighbors[Idents(datasets)[neighbors]==tissues[2]]
               t22 <- c(t22,temp)}
         t0 <- c(t0,t11,t22)
         t0 <- unique(t0)
         t1 <- t11
         t2 <- t22
         i=i+1}
     return(t0)}

plot_root_shape <- function(datasets, imagename, cellss,tissues)
    {
    all_tissues <- c("Root_cap","Meristem_1","Meristem_2", "Stele_1","Stele_2", "Stele_3", "Endodermis","Cortex","Meristem_3","Epidermis_exodermis_sclerenchyma","Trasition_zone")
    mypalette <- mypalette<-c("#8ECFC9","#E64B35","#FA7F6F","#00441B","#7E6148", "#3C5488", "#EAFC1C","#74C476","#82B0D2","#4e8781","#FFBE7A")
    sub_datasets<-subset(x=datasets,cells=cellss)
    exp<-data.frame(cell.type=Idents(sub_datasets),cellid=names(Idents(sub_datasets)))     
    temp0 <- datasets@images[[imagename]]@coordinates[,1:2]
    temp0$cellid <- rownames(temp0)
    temp<-left_join(temp0,exp,by = 'cellid')
    temp$cell.type<- ordered(temp$cell.type, levels = tissues)
    plot2<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = cell.type),size = 1)+
        scale_colour_discrete(type=mypalette[all_tissues %in% tissues],na.value="grey")+
        ggtitle(imagename)+
        theme(panel.grid=element_blank())+
        coord_fixed() 
    #ggsave2("root_shape.pdf",plot2,width=4,height=7)
    return(plot2)
    }

tissues <- unique(Idents(BRAIN))
tissues <- tissues[!is.na(tissues)]
all_tissue <- combn(tissues,2)
for (i in 1:length(all_tissue))
    {tar_ti <- all_tissue[,i]
     cells <- find_neighbor(datasets=BRAIN,tissues=tar_ti,3)
     if(length(cells)>1)
         {p<-list()
         for(imagename in names(BRAIN@images))
              {p[[imagename]]<- plot_root_shape(datasets=BRAIN, imagename=imagename, cellss=cells,tissues=tar_ti)}
         sum<-grid.arrange(grobs=p,ncol=7)
         ggsave2(paste0(tar_ti[1],"_",tar_ti[2],"_dimspatial.pdf"),sum,width=40,height=20)
         sub_BRAIN <- subset(x=BRAIN,cells=cells)
         temp  <- try(FindMarkers(sub_BRAIN,ident.1 = tar_ti[1],ident.2 = tar_ti[2],min.pct = -Inf,logfc.threshold = opt$logfc),silent=FALSE)
         if(!'try-error' %in% class(temp))
             {trans_genes<-FindMarkers(sub_BRAIN, 
                          ident.1 = tar_ti[1], 
                          ident.2 = tar_ti[2], 
                          min.pct = -Inf,
                          logfc.threshold = opt$logfc,
                          verbose = FALSE)
              trans_genes<-trans_genes[trans_genes$p_val<0.01,]
              trans_genes$gene <- substring(rownames(trans_genes),1,12)
              trans_genes.info <- left_join(trans_genes,genename,by='gene')
              write.csv(trans_genes.info,file=paste0(tar_ti[1],"_",tar_ti[2],"_markers.csv"),quote=F)}}
       }
               






