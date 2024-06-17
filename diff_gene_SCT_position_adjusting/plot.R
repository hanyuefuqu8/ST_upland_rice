library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)

dataset<- readRDS("data_pos_adj.RDS")
Idents(dataset)<- ordered(Idents(dataset), levels = c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Xylem","Phloem", "Parenchyma_1","Parenchyma_2", "Pericycle_like", "Ground_tissue","Ground_tissue_2","Ground_tissue_3"))

mypalette<-c(brewer.pal(4,"Reds")[4:1],"#00441B","#4DBBD5","#3C5488","#8491B4","#EAFC1C","#74C476","#c49374","#b9c474")
p<-list()
for(i in names(dataset@images))
{
p[[i]]<- SpatialDimPlot(dataset,images=i,stroke=0,pt.size=20,label=F)+scale_fill_manual(breaks=levels(Idents(dataset)),values=mypalette)+ggtitle(i)+theme(legend.position="none")
}
sum<-grid.arrange(grobs=p,ncol=7)
ggsave2("dimspatial_change.pdf",sum,width=40,height=20)
ggsave2("dimspatial_change.png",sum,width=40,height=20)

