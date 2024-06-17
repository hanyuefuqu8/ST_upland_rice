library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(monocle3)
library(readr)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(gridExtra)

cds <- readRDS("data.RDS")
BRAIN0<-readRDS("../../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")
CR_selected <- read.csv("../../CR_select/CR_selected.csv",header=T,quote="",row.names=1)
cr_selected <- rownames(CR_selected)
BRAIN<-subset(x=BRAIN0,cells=cr_selected)
rdsf <- read.table("../../info2",header=T)
temp <- left_join(BRAIN@meta.data,rdsf,by='group')
BRAIN$species <- temp$species
BRAINs1 <- subset(BRAIN,subset=species==unique(colData(cds)$species))
pseudotime <- read.csv("pseudotime.csv",header=T,row.names=1,quote="")
BRAINs1$pseudotime <- pseudotime


BRAIN.list <- SplitObject(BRAIN0, split.by = "group")
for( i in 1:length(BRAIN.list))
   {BRAIN.list[[i]] <- SCTransform(BRAIN.list[[i]], vst.flavor = "v2",assay='Spatial',verbose = FALSE)}
BRAIN2<-merge(BRAIN.list[[1]], y=BRAIN.list[-1])
BRAIN2<-PrepSCTFindMarkers(BRAIN2)

genename<-read.table("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/upland_rice/all_root2/harmony/Markergenes/gene_v3.txt",sep="\t",quote="",header=T)
colnames(genename) <- c("gene","name")
genename <- unique(genename)

single_neighbor <- function(datasets,cell_id)
    {dist <- sort(abs(BRAINs1$pseudotime-BRAINs1@meta.data[cell_id,'pseudotime']))
     cells <- names(dist[2:3])
    return(cells)}
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
    {all_tissues <- c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Vascular", "Parenchyma_1","Parenchyma_2", "Pericycle_like", "Ground_tissue","Ground_tissue_2","Ground_tissue_3")
    mypalette<-c(brewer.pal(4,"Reds")[4:1],"#00441B","#3C5488","#8491B4","#EAFC1C","#74C476","#c49374","#b9c474")    
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

tissues <- unique(Idents(BRAINs1))
all_tissue <- combn(tissues,2)
all_tissues <- c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Vascular", "Parenchyma_1","Parenchyma_2", "Pericycle_like", "Ground_tissue","Ground_tissue_2","Ground_tissue_3")
mypalette<-c(brewer.pal(4,"Reds")[4:1],"#00441B","#3C5488","#8491B4","#EAFC1C","#74C476","#c49374","#b9c474")

for (i in 1:length(all_tissue))
    {tar_ti <- all_tissue[,i]
     cells <- find_neighbor(datasets=BRAINs1,tissues=tar_ti,3)
     ##plot these cells spatial
     if(length(cells)>1)
         {p<-list()
         for(imagename in names(BRAIN@images))
              {p[[imagename]]<- plot_root_shape(datasets=BRAIN, imagename=imagename, cellss=cells,tissues=tar_ti)}
         sum<-grid.arrange(grobs=p,ncol=7)
         ggsave2(paste0(tar_ti[1],"_",tar_ti[2],"_dimspatial.pdf"),sum,width=40,height=20)}
     ##plot these cells umap
     select_cell <- rep("N",length(Cells(BRAINs1)))
     names(select_cell) <- Cells(BRAINs1)
     select_cell[cells] <- BRAINs1@meta.data[cells,'cell.type']%>% as.character
     colData(cds)$select <- select_cell 
     colData(cds)$select <- ordered(colData(cds)$select, levels = c(as.character(tar_ti), "N")) 
     cols <- c(mypalette[all_tissues %in% tar_ti],"grey")
     p <- plot_cells(cds,color_cells_by = "select",cell_size = 1)+ scale_colour_discrete(type=cols)
     ggsave2(paste0(tar_ti[1],"_",tar_ti[2],"_umap.pdf"),p)
     ##divide the spot into 2 parts
     sub_BRAIN <- subset(x=BRAINs1,cells=cells)
     sub_BRAIN$pseudotime2 <- as.numeric(sub_BRAIN$pseudotime >= median(sub_BRAIN$pseudotime))
     select_cell[cells] <- sub_BRAIN$pseudotime2
     colData(cds)$select <- select_cell
     colData(cds)$select <- ordered(colData(cds)$select, levels = c("0","1", "N"))
     cols <- c(mypalette[all_tissues %in% tar_ti],"grey")
     p <- plot_cells(cds,color_cells_by = "select",cell_size = 1)+ scale_colour_discrete(type=cols)
     ggsave2(paste0(tar_ti[1],"_",tar_ti[2],"_umap_2parts.pdf"),p)
     Idents(sub_BRAIN) <- sub_BRAIN$pseudotime2
     cell1 <- sub_BRAIN@meta.data[sub_BRAIN@meta.data$pseudotime2=="0",] %>% rownames
     cell2 <- sub_BRAIN@meta.data[sub_BRAIN@meta.data$pseudotime2=="1",] %>% rownames
     temp  <- try(FindMarkers(BRAIN2,ident.1 = cell1,ident.2 = cell2,min.pct = -Inf,logfc.threshold = 0.05),silent=FALSE)
         if(!'try-error' %in% class(temp))
             {trans_genes<-FindMarkers(BRAIN2, 
                          ident.1 = cell1, 
                          ident.2 = cell2, 
                          min.pct = -Inf,
                          logfc.threshold = 0.05,
                          verbose = FALSE)
              trans_genes<-trans_genes[trans_genes$p_val<0.01,]
              trans_genes$gene <- substring(rownames(trans_genes),1,12)
              trans_genes.info <- left_join(trans_genes,genename,by='gene')
              write.csv(trans_genes.info,file=paste0(tar_ti[1],"_",tar_ti[2],"_markers.csv"),quote=F)}}

