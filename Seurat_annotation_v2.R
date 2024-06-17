library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gdata)
library(RColorBrewer)
data<-readRDS("data.RDS")

data<-RenameIdents(data, 
`0` = "Ground_tissue",
 `1` = "Pericycle_like",
 `2` = "Parenchyma_2",
 `3` = "Crown_root_primordia",
 `4` = "Ground_tissue",
 `5` = "Xylem",
 `6` = "transition_tissue_1",
 `7` = "Ground_tissue",
 `8` = "Ground_tissue_2",
 `9` = "Parenchyma_1", 
`10` = "Ground_tissue",
`11` = "Ground_tissue",
`12` = "transition_tissue_2",
`13` = "Crown_root_primordia",
`14` = "Ground_tissue",
`15` = "transition_tissue_3",
`16` = "Ground_tissue",
`17` = "Parenchyma_2")
Idents(data)<- ordered(Idents(data), levels = c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Xylem", "Parenchyma_1","Parenchyma_2", "Pericycle_like", "Ground_tissue","Ground_tissue_2","Ground_tissue_3"))

mypalette<-c(brewer.pal(4,"Reds")[4:1],"#00441B","#3C5488","#8491B4","#EAFC1C","#74C476","#c49374","#b9c474")

library(gridExtra)

pdf("09_umap_3.pdf", width = 12,height=10)
plot1 <- DimPlot(object = data, reduction = "umap", label = T,label.size=5,pt.size=1,group.by = "ident",cols=mypalette)
plot1
dev.off()

data$cell.type <- Idents(data)
plot2<-VlnPlot(data,features="nCount_SCT",group.by="cell.type",pt.size = 0,log=TRUE)
ggsave2("QCVlnplot1.pdf",plot2,width=7)
plot2<-VlnPlot(data,features="nFeature_SCT",group.by="cell.type",pt.size = 0,log=TRUE)
ggsave2("QCVlnplot2.pdf",plot2,width=7)


names(data@images) <- gsub("-",".",names(data@images))
datas1 <- subset(x = data, idents = "Parenchyma_2")

datas1 <- FindNeighbors(datas1, dims = 1:20)
datas1 <- FindClusters(datas1, resolution = 0.4, verbose = FALSE)
plot2 <- DimPlot(object = datas1, reduction = "umap", label = T,label.size=10,pt.size=1,group.by = "ident")
ggsave2("subset1_recluster.pdf",plot2)

p2<-DotPlot(datas1, features = c("Os08g0162800","Os03g0170900","Os01g0919800"),assay="SCT")
ggsave2("phloem_marker_bubble.pdf",p2,width=4,height=4)
#s03g0170900,OsSUT1
#Os01g0919800,OsPIN5a

Phloem <- WhichCells(datas1,idents=2)
data@meta.data$cell.type<-as.character(data@meta.data$cell.type)
data@meta.data[rownames(data@meta.data)%in%Phloem,'cell.type'] <- "Phloem"
Idents(data)<-data$cell.type
Idents(data)<- ordered(Idents(data), levels = c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Xylem","Phloem","Parenchyma_1","Parenchyma_2", "Pericycle_like", "Ground_tissue","Ground_tissue_2","Ground_tissue_3"))

mypalette<-c(brewer.pal(4,"Reds")[4:1],"#00441B","#4DBBD5","#3C5488","#8491B4","#EAFC1C","#74C476","#c49374","#b9c474")

library(gridExtra)
p<-list()
for(i in names(data@images))
{
p[[i]]<- SpatialDimPlot(data,images=i,stroke=0,pt.size=12,label=F)+scale_fill_manual(breaks=levels(Idents(data)),values=mypalette)+ggtitle(i)+ theme(legend.position="none")
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2("dimspatial.pdf",sum,width=20,height=20)
ggsave2("dimspatial.png",sum,width=20,height=20)

pdf("09_umap_4.pdf", width = 8,height=5)
plot1 <- DimPlot(object = data, reduction = "umap", label = T,label.size=4,pt.size=0.5,group.by = "ident",cols=mypalette)
plot1
dev.off()


saveRDS(data,"data2.RDS")
