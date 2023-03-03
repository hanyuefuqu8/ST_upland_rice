library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
                "info", "f","1","character","information file",
		        "node", "n","1","character","cluster contains node",
		        "graph","p","1","character","choose knn or principal_graph"
                ),byrow=T,ncol=5)
opt = getopt(arg)



library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(monocle3)
library(readr)
library(cowplot)
if(FALSE)
{
opt<-list()
opt$input<-'../sub.RDS'
opt$info<-'../info2'
opt$node<-'Meristem_1'
opt$graph<-'knn'
}
BRAIN<-readRDS(opt$input)
rdsf <- read_tsv(opt$info)
temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$phase<-temp$phase
BRAIN@meta.data$cell.type<-paste0("phase",BRAIN@meta.data$phase,"_",Idents(BRAIN))

#subset1
#BRAINs1<-subset(x=BRAIN,idents=c("Procambium","Xylem","Preprocambium","Stele initials"))
data <- as(as.matrix(BRAIN@assays$Spatial@counts), 'sparseMatrix')
pd <- data.frame(BRAIN@meta.data,Idents(BRAIN))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)
#remove batch effect
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "group")


#reduce dimensions
cds <- reduce_dimension(cds)

#determine how many dimensions to use
pdf("00_dim.pdf")
plot_pc_variance_explained(cds)
dev.off()

#check batch
pdf("00_batch.pdf")
plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE)
dev.off()

pdf("01_UMAP.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "Idents.BRAIN.",group_label_size = 5,cell_size = 1, label_cell_groups=FALSE)
dev.off()

#clustering
cds <- cluster_cells(cds)
pdf("02_partition.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "partition")
dev.off()

#learn graph
cds <- learn_graph(cds)
pdf("03_trajectorise.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "Idents.BRAIN.",
           label_groups_by_cluster=TRUE,
           cell_size = 1,
           label_leaves=FALSE,
           label_branch_points=FALSE,
		   label_cell_groups=FALSE)
dev.off()

pdf("03_trajectorise_phase.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "phase",
           label_groups_by_cluster=TRUE,
           cell_size = 1,
           label_leaves=FALSE,
           label_branch_points=FALSE,
		   label_cell_groups=FALSE)
dev.off()

pdf("03_trajectorise_cell.type.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=TRUE,
           cell_size = 0.7,
           label_leaves=FALSE,
           label_branch_points=FALSE,
		   label_cell_groups=FALSE)
dev.off()

get_earliest_principal_node <- function(cds, time_bin=opt$node){

  cell_ids <- which(colData(cds)[,"Idents.BRAIN."] == time_bin)
 
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]
 
  root_pr_nodes
}

cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

pdf("trajectorise_set_root.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "pseudotime",
		   cell_size = 1,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()

pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
BRAIN@meta.data$pseudotime<-pseudotime


#change the first name of images
#b<-paste0("image_",gsub("-",".",names(table(BRAIN$group))))
#c<- names(BRAIN@images)
#names(BRAIN@images)[1]<-b[!b%in%c]


library(gridExtra)
p<-list()
for(i in names(table(BRAIN$group)))
{
imagename<-paste0("image_",gsub("-",".",i))
temp<-cbind(BRAIN@meta.data[BRAIN@meta.data$group==i,c('group','pseudotime')],BRAIN@images[[imagename]]@coordinates[,1:2])
#h <- as.numeric(max(temp$y) - min(temp$y) + 1)
#w <- as.numeric(max(temp$x) - min(temp$x) + 1)
p[[i]]<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = pseudotime),size = 6)+
  scale_color_gradient(
  #limit=c(0,max(BRAIN@meta.data$pseudotime)),
  low = "blue",high = "#eaed18")+
  ggtitle(i)

}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2("pseudotime_spatial.png",sum,dpi=600,width=25,height=20)


#finding genes changing along pseudotime
#The data frame pr_graph_test_res has the Moran's I test results for each gene in the cell_data_set. 
#If you'd like to rank the genes by effect size, sort this table by the morans_Icolumn, which ranges from -1 to +1. 
#A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells have very similar values of a gene's expression. 
#Significant values much less than zero are generally rare.
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph=opt$graph, cores=4)
top_genes<-ciliated_cds_pr_test_res %>% top_n(n = 12, wt = morans_I)
pdf("top_genes.pdf")
plot_cells(cds, genes=top_genes$gene_short_name,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()


for (i in top_genes$gene_short_name)
{
library(gridExtra)
p<-list()
for(j in names(table(BRAIN$group)))
{
imagename<-paste0("image_",gsub("-",".",j))
p[[j]]<- SpatialFeaturePlot(BRAIN,images=imagename,features=i,stroke=0,pt.size=12)+ggtitle(imagename)

}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2(paste0(i,"sub_dimspatial.png"),sum,dpi=600,width=20,height=20)

}

pr_deg <- subset(ciliated_cds_pr_test_res, q_value < 0.01)
pr_deg_ids<-rownames(pr_deg)
write.csv(file="pr_deg.csv",pr_deg,sep=",",quote=FALSE)

#collect the trajectory-variable genes into modules:
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-3)
write.table(gene_module_df, "gene_module_df.txt",col.names = F)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell.type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pdf("pr_deg_heatmap.pdf")
plot1<-pheatmap::pheatmap(agg_mat,scale="row", clustering_method="ward.D2")
plot1
dev.off()
ggsave2("pr_deg_heatmap.png",plot1,dpi=600)

write.csv(file="agg_mat.csv",agg_mat,sep=",",quote=FALSE)


