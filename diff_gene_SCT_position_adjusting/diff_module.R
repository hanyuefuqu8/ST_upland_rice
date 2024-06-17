library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(readr)
library(monocle3)
BRAIN<-readRDS("data_pos_adj.RDS")
genename<-read.csv("Irrigated_riceUpland_rice_total_diffgene.csv",quote="",header=T,row.names=1)
DefaultAssay(object = BRAIN) <- "SCT"

BRAIN$cell.type <- Idents(BRAIN)
rdsf <- read.table("../info2",sep="\t",header=T)
temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$species<-temp$species
BRAIN$celltype.treat <- paste(BRAIN$species, Idents(BRAIN), sep = "_")

BRAIN.list <- SplitObject(BRAIN, split.by = "group")
for( i in 1:length(BRAIN.list))
{
BRAIN.list[[i]] <- SCTransform(BRAIN.list[[i]], vst.flavor = "v2",assay='Spatial',verbose = FALSE)
}
BRAIN<-merge(BRAIN.list[[1]], y=BRAIN.list[-1])
BRAIN<-PrepSCTFindMarkers(BRAIN)
DefaultAssay(BRAIN) <- "SCT"


#data <- as(as.matrix(BRAIN@assays$Spatial@counts), 'sparseMatrix')
data <- as(as.matrix(BRAIN@assays$SCT@data), 'sparseMatrix')
pd <- data.frame(BRAIN@meta.data,Idents(BRAIN))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)
cds <- preprocess_cds(cds, num_dim = 100,norm_method = "none")
cds <- align_cds(cds, alignment_group = "group")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)


pr_deg_ids <- unique(genename$gene)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-3)
write.csv(gene_module_df, "gene_module_df.csv",quote=F)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$celltype.treat)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
count1 <- c(1:ncol(agg_mat))
count2 <- c(1:ncol(agg_mat))
count2[count1%%2==0]<-c((length(count1)/2+1):length(count1))
count2[count1%%2==1]<-c(1:(length(count1)/2))
agg_mat <- agg_mat[,count2]
p<-pheatmap::pheatmap(agg_mat,scale="row", clustering_method="ward.D2",  cluster_cols = FALSE)
ggsave2(file="pr_deg_heatmap.pdf",p,device='pdf',width=6,height=6)


