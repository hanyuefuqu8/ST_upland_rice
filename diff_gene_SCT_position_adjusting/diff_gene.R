library(getopt)

arg <- matrix(c("input", "i","1","character","input file1",
                "info","f","1","character","information file",
		"logfc","l","1","numeric","log fold change to filter the markers",
                "method","m","1","character","method to normalize data"
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
if(FALSE)
{
opt <- list()
opt$input <- "data_pos_adj.RDS"
opt$info <- "../info2"
opt$logfc <- 0.01
opt$method <- "SCT" 
}
BRAIN<-readRDS(opt$input)
BRAIN<-subset(BRAIN,idents="Unknown",invert=TRUE)
nacell <- Idents(BRAIN)[is.na(Idents(BRAIN))] %>% names()
BRAIN<-subset(BRAIN,cells=nacell,invert=TRUE)
rdsf <- read.table(opt$info,sep="\t",header=T)

temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$chip_name<-temp$chip_name
BRAIN@meta.data$line<-temp$line
BRAIN@meta.data$line2<-temp$line2
BRAIN@meta.data$species<-temp$species
BRAIN@meta.data$cell.type<-Idents(BRAIN)
library(RColorBrewer)
mypalette<-c(brewer.pal(4,"Reds")[4:1],"#00441B","#4DBBD5","#3C5488","#8491B4","#EAFC1C","#74C476","#c49374","#b9c474")
Idents(BRAIN)<- ordered(Idents(BRAIN), levels = c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Xylem","Phloem","Parenchyma_1","Parenchyma_2", "Pericycle_like", "Ground_tissue","Ground_tissue_2","Ground_tissue_3"))

pdf("09_umap.pdf", width = 10, height = 10)
plot1 <- DimPlot(object = BRAIN, reduction = "umap", label = F,pt.size=0.5,group.by = "ident", split.by = "chip_name",cols=mypalette)
plot2 <- DimPlot(object = BRAIN, reduction = "umap", label = F,pt.size=0.5,group.by = "ident", split.by = "species",cols=mypalette)
plot3 <- DimPlot(object = BRAIN, reduction = "umap", label = F,pt.size=0.5, group.by = "ident", split.by = "line2",cols=mypalette)
#plot2 <- SpatialDimPlot(BRAIN, stroke = 0, pt.size=2, label=T)
plot1/plot2/plot3
dev.off()

BRAIN$celltype.treat <- paste(BRAIN$species, Idents(BRAIN), sep = "_")
Idents(BRAIN) <- "celltype.treat"



plot1<-VlnPlot(BRAIN,features="nCount_Spatial",group.by="celltype.treat")
ggsave2("Violin_Spatial.png",plot1,width=14)
plot1<-VlnPlot(BRAIN,features="nCount_SCT",group.by="celltype.treat")
ggsave2("Violin_SCT.png",plot1,width=14)
plot1<-VlnPlot(BRAIN,features="nCount_Spatial",group.by="line")
ggsave2("Violin_Line_Spatial.png",plot1,width=14)
plot1<-VlnPlot(BRAIN,features="nCount_SCT",group.by="line")
ggsave2("Violin_Line_SCT.png",plot1,width=14)


if (opt$method=="RC")
{
DefaultAssay(BRAIN) <- "Spatial"
BRAIN <- NormalizeData(BRAIN, normalization.method = "RC", scale.factor = 10000)
BRAIN@meta.data$RC_nCount_Spatial<-colSums(BRAIN[['Spatial']]@data)
plot1<-VlnPlot(BRAIN,features="RC_nCount_Spatial",group.by="celltype.treat")
ggsave2("Violin_RC.png",plot1,width=14)
} else if(opt$method=="log")
{
DefaultAssay(BRAIN) <- "Spatial"
BRAIN <- NormalizeData(BRAIN, normalization.method = "LogNormalize", scale.factor = 10000)
BRAIN@meta.data$log_nCount_Spatial<-colSums(BRAIN[['Spatial']]@data)
plot1<-VlnPlot(BRAIN,features="log_nCount_Spatial",group.by="celltype.treat")
ggsave2("Violin_log.png",plot1,width=14)
} else if(opt$method=="SCT")
{
BRAIN.list <- SplitObject(BRAIN, split.by = "group")
for( i in 1:length(BRAIN.list))
{
BRAIN.list[[i]] <- SCTransform(BRAIN.list[[i]], vst.flavor = "v2",assay='Spatial',verbose = FALSE)
}
BRAIN<-merge(BRAIN.list[[1]], y=BRAIN.list[-1])
BRAIN<-PrepSCTFindMarkers(BRAIN)
DefaultAssay(BRAIN) <- "SCT"
}
plot1<-VlnPlot(BRAIN,features="nCount_SCT",group.by="line",pt.size=0.1)
ggsave2("Violin_Line_SCTv2.pdf",plot1,width=7,height=4)


p<-list()
for(k in names(table(BRAIN$cell.type)))
    {
    a<-paste0(names(table(rdsf$species))[1],"_",k)
    b<-paste0(names(table(rdsf$species))[2],"_",k)
       if(a%in%Idents(BRAIN) && b%in%Idents(BRAIN))
           {
    p[[k]]<-FindMarkers(BRAIN, 
                          ident.1 = a, 
                          ident.2 = b, 
                          min.pct = -Inf,
                          logfc.threshold = opt$logfc,
                          verbose = FALSE)
    p[[k]]$cluster<-k
    p[[k]]$gene <- rownames(p[[k]])
    #plots <- VlnPlot(BRAIN, features = rownames(p[[k+1]])[1:5], split.by = "species", group.by = "celltype",pt.size = 0, combine = FALSE)
    #sum<-wrap_plots(plots = plots, ncol = 1)
    #h<-paste0(k,"_diffgene.png")
    #ggsave2(h,sum,height=20)
    #write.table(p[[k]],"total_diffgene.csv",col.names = F, append = T,quote=F,sep=",")
           }
  }

total <- rbind(p[[1]],p[[2]])
for(i in 3:length(p))
    {total <- rbind(total,p[[i]])}
genename<-read.table("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/upland_rice/all_root2/harmony/Markergenes/gene_v3.txt",sep="\t",quote="",header=T)
colnames(genename) <- c("gene","name")
genename <- unique(genename)
total<-total[total$p_val<0.01,]
total.info <- left_join(total,genename,by='gene')
write.table(total.info,paste0(names(table(rdsf$species))[1],names(table(rdsf$species))[2],"_total_diffgene.csv"),quote=F,sep=",")





