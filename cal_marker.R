library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
                "logfc","l","1","numeric","log fold change to filter the markers",
                "geneid","g","1","character","geneid table"
                ),byrow=T,ncol=5)
opt = getopt(arg)

library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(readr)
BRAIN<-readRDS(opt$input)
genename<-read.table(opt$geneid,sep="\t",quote="",header=T)
colnames(genename)<-c("gene","ID")
DefaultAssay(object = BRAIN) <- "SCT"


BRAIN<-FindVariableFeatures(BRAIN, nfeatures = 10000)
BRAIN.markers <- FindAllMarkers(object = BRAIN,
                                test.use = "wilcox",
                                 only.pos = T,
                                 min.pct = -Inf,
                                 logfc.threshold = opt$logfc)
write.table(BRAIN.markers,file="markers.txt",sep="\t",quote=F)

top10 <- BRAIN.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
tops10<-top10[top10$gene%in%rownames(BRAIN[['SCT']]@scale.data),]
topi10<-left_join(top10,genename,by="gene")
topi10<-topi10[!duplicated(topi10$gene),]


pdf ("marker_heatmap.pdf",width=9)
DoHeatmap(BRAIN, features = tops10$gene)
dev.off()
pdf("marker_dotplot.pdf",height=5,width=8)
p1<-DotPlot(BRAIN, features = unique(topi10$gene),assay="SCT")+
           scale_x_discrete(labels=topi10$ID)+
		   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
		   scale_colour_gradient(low = "#c6e6e8", high = "#fa5d19")
p1
dev.off()


dir.create("images")
setwd("images")
for(j in unique(topi10$gene))
{
library(gridExtra)
p<-list()
exp<-data.frame(FetchData(object = BRAIN, vars = j))


for(i in names(table(BRAIN$group)))
{

imagename<-paste0("image_",i)

temp<-merge(BRAIN@images[[imagename]]@coordinates[,1:2],exp,by = 'row.names')
colnames(temp)[4]<-"expression"
#h <- as.numeric(max(temp$y) - min(temp$y) + 1)
#w <- as.numeric(max(temp$x) - min(temp$x) + 1)
p[[i]]<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = expression),size = 6)+
  scale_color_gradient(limit=c(0,max(temp$expression)),low = "blue",high = "#eaed18")+
  ggtitle(i)

}
sum<-grid.arrange(grobs=p,ncol=4)
name<-topi10[topi10$gene==j,'ID'][1,1]
ggsave2(paste0(j,".pdf"),sum,width=25,height=20)

}

