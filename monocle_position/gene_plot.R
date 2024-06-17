library(Seurat)
library(viridis)
library(monocle3)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(stringr)
datasets <- readRDS("../../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")
tissues <- c("Meristem_2","Endodermis","Cortex","Trasition_zone")

load("../align.RData")
datasets <- AddMetaData(object = datasets, metadata = phases,col.name = 'phase')
datasets <- AddMetaData(object = datasets, metadata = lengths,col.name = 'length')

gene.list <- read.csv("../gene.list",header=T,quote="")
dir.create("images2")
setwd("images2")


BRAIN.list <- SplitObject(datasets, split.by = "group")
for( i in 1:length(BRAIN.list))
{
BRAIN.list[[i]] <- SCTransform(BRAIN.list[[i]], vst.flavor = "v2",assay='Spatial',verbose = FALSE)
}
BRAIN<-merge(BRAIN.list[[1]], y=BRAIN.list[-1])
BRAIN<-PrepSCTFindMarkers(BRAIN)
DefaultAssay(BRAIN) <- "SCT"
BRAIN@images <- datasets@images

dataset <- subset(x=BRAIN,idents=tissues)
dataset$group_phase <- paste0(dataset$group,"_",dataset$phase)
rdsf <- read.table("../../../info2",sep="\t",header=T)

for(i in gene.list[,'geneid'])
    {
     if (i %in% rownames(dataset[["SCT"]]@data))
     {
     exp <- AverageExpression(dataset,features=i,slot='counts',group.by = "group_phase")
     temp <- DataFrame(t(exp$SCT))
     temp$group <- str_split_fixed(rownames(temp),"_",2)[,1]
     temp$phase <- str_split_fixed(rownames(temp),"_",2)[,2]
     temp <- as.data.frame(temp)
     temp<-left_join(temp,rdsf[,c('group','species')],by="group")
     write.csv(file=paste0(gene.list[gene.list$geneid==i,'genename'],".csv"),temp,quote=F,row.names=F)
     #temp <- read.csv(paste0(gene.list[gene.list$geneid==i,'genename'],".csv"),header=T,quote="")
     p <- ggplot(data = temp, mapping = aes(x = as.numeric(phase), y = V1,group=group, color=species)) + 
     geom_line() + 
     scale_colour_manual("",values=c("Irrigated_rice"="light blue","Upland_rice"="Light Salmon"))+
     stat_summary(data=subset(temp,species == 'Irrigated_rice'),aes(group=1),fun.y=mean, geom="line", size=1.2,color="blue")+
     stat_summary(data=subset(temp,species == 'Upland_rice'),aes(group=1),fun.y=mean, geom="line", size=1.2,color="red")+
     theme_classic()
     ggsave2(file=paste0(gene.list[gene.list$geneid==i,'genename'],".pdf"),p,width=5,height=2)
     }
     }





