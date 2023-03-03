library(dplyr)
library(gdata)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(stringr)
library(ggrepel)

BRAIN<-readRDS("../data2.RDS")
rdsf <- read.table("../diff_gene_RC/info",sep="\t",header=T)
temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$chip_name<-temp$chip_name
BRAIN@meta.data$line<-temp$line
BRAIN@meta.data$species<-temp$species
BRAIN@meta.data$cell.type<-Idents(BRAIN)
Idents(BRAIN)<-BRAIN@meta.data$species

BRAIN.list <- SplitObject(BRAIN, split.by = "group")
for( i in 1:length(BRAIN.list))
{
BRAIN.list[[i]] <- SCTransform(BRAIN.list[[i]], vst.flavor = "v2",assay='Spatial',verbose = FALSE)
}
BRAIN<-merge(BRAIN.list[[1]], y=BRAIN.list[-1])
BRAIN<-PrepSCTFindMarkers(BRAIN)
DefaultAssay(BRAIN) <- "SCT"

dataset<-read.table("all_sig_diff.xls",header=T,fill=T,row.names='X',sep='\t',quote="")
BRAIN<-BRAIN[dataset$ID,]
mk<-FindMarkers(BRAIN, 
                          ident.1 = "Irrigated_rice", 
                          ident.2 = "Upland_rice", 
                          min.pct = -Inf,
                          logfc.threshold = 0.01,
                          verbose = FALSE)


write.table(file="temp.txt",mk,sep="\t")
dataset<-mk
cut_off_pvalue = 0.01  #统计显著性
cut_off_logFC = 0.05         #差异倍数值
dataset$change = ifelse(dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_log2FC) >= cut_off_logFC, 
                          ifelse(dataset$avg_log2FC > cut_off_logFC ,'Irrigated_rice','Upland_rice'),
                          'Stable')

dataset$X<-rownames(dataset)
Gene<-c("Os04g0476100")
p <- ggplot(dataset, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=change)) +
  geom_point(alpha=0.8, size=1) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="avg_log2FC",
       y="-log10(p_val_adj)")+
  theme_bw()+
 #标记基因
    geom_text_repel(
    data = dataset[dataset$X%in%Gene,],
    aes(label = X),
    size = 1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
ggsave2("volcano.pdf",p,height=7,width=4)

pdf("stat.pdf",height=4,width=4)
p2 <- ggplot(dataset,aes(x = avg_log2FC))+geom_histogram(binwidth=0.1,fill = 'steelblue',colour = 'darkred')+ 
 theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
p2
dev.off()

library(aplot)
pdf("aligh_volcano.pdf",height=10,width=4)
p3<-p%>%insert_top(p2,height=0.3)
p3
dev.off()





