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

BRAIN<-readRDS(opt$input)
rdsf <- read.table(opt$info,sep="\t",header=T)

temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$chip_name<-temp$chip_name
BRAIN@meta.data$line<-temp$line
BRAIN@meta.data$species<-temp$species
BRAIN@meta.data$cell.type<-Idents(BRAIN)
BRAIN$celltype.treat <- paste(BRAIN$species, Idents(BRAIN), sep = "_")
Idents(BRAIN) <- "celltype.treat"

pdf("09_umap.pdf", width = 16, height = 16)
plot1 <- DimPlot(object = BRAIN, reduction = "umap", label = T,pt.size=0.5,group.by = c("cell.type", "chip_name"))
plot2 <- DimPlot(object = BRAIN, reduction = "umap", label = T,pt.size=0.5,group.by = c("line", "species"))
#plot2 <- SpatialDimPlot(BRAIN, stroke = 0, pt.size=2, label=T)
plot1/plot2
dev.off()

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
    #plots <- VlnPlot(BRAIN, features = rownames(p[[k+1]])[1:5], split.by = "species", group.by = "celltype",pt.size = 0, combine = FALSE)
    #sum<-wrap_plots(plots = plots, ncol = 1)
    #h<-paste0(k,"_diffgene.png")
    #ggsave2(h,sum,height=20)
    write.table(p[[k]],"total_diffgene.csv",col.names = F, append = T,quote=F,sep=",")
           }
  }

data<-read.csv("total_diffgene.csv",sep=",",header=F)
#d<-quantile(data$V2,0.15)
#data<-data[data$V2<d,]

#闁荤姳绶ょ槐鏇㈡偩婵犳氨宓侀柛鎰级缂嶅棝鎮跺☉鏍у闁规椿浜弻宀冪疀閹剧懓鐝紓鍌欑劍閿氶柛鈺傤殜閹粓顢橀悙瀵割洯
BRAIN <- NormalizeData(BRAIN, normalization.method = "LogNormalize", scale.factor = 10000)
temp<-data[!duplicated(data$V1),]
a<-temp[,1]%in%rownames(BRAIN[["Spatial"]]@data)
temp<-temp[a,]
temp2<-temp[,c(1,7)]
#p<-AverageExpression(BRAIN,features=temp[,1],group.by = "ident")
p<-AverageExpression(BRAIN,features=temp[,1],slot='counts',group.by = "ident")
if (opt$method=="SCT")
{
temp3<-as.data.frame(p$SCT)
}else if(opt$method=="log")
{
temp3<-as.data.frame(p$Spatial)
}else if(opt$method=="RC")
{
temp3<-as.data.frame(p$Spatial)
}

#temp3<-temp3[,order(rownames(temp3))]
pdf("total_diffgene_heatmap.pdf",height=10,width=20)
p1<-pheatmap::pheatmap(t(temp3),
                       scale="column",
                       cluster_rows= FALSE,
                       #gaps_row = seq(from=2,by=2,length=8),
                       #cutree_cols = 6
					   )
print(p1)
dev.off()
ggsave2("total_diffgene_heatmap.png",p1,dpi=600,width=40)



