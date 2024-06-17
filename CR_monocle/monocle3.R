library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(monocle3)
library(readr)
library(cowplot)
library(viridis)
library(RColorBrewer)
if(TRUE)
{
opt<-list()
opt$input<-'../diff_gene_SCT_position_adjusting/data_pos_adj.RDS'
opt$info<-'../info2'
opt$node1<-'Pericycle_like'
opt$node2<-'Crown_root_primoridia'
opt$graph<-'principal_graph'
}
info <- read.table("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/upland_rice/all_root2/harmony/Markergenes/gene_v3.txt",quote="",header=T,sep
="\t")
colnames(info) <- c("id","gene_name")

#tissues <- c("Pericycle_like","transition_tissue_1","transition_tissue_2","transition_tissue_3","Crown_root_primordia")
CR_selected <- read.csv("../CR_select/CR_selected.csv",header=T,quote="",row.names=1)
cr_selected <- rownames(CR_selected)
BRAIN0<-readRDS(opt$input)
BRAIN<-subset(x=BRAIN0,cells=cr_selected)

Idents(BRAIN)<- ordered(Idents(BRAIN), levels = c("Crown_root_primordia","transition_tissue_1","transition_tissue_2","transition_tissue_3","Pericycle_like", "Ground_tissue_2"))
mypalette<-c(brewer.pal(4,"Reds")[4:1],"#EAFC1C","#c49374")
rdsf <- read.table(opt$info,header=T)
temp <- left_join(BRAIN@meta.data,rdsf,by='group')
BRAIN$species <- temp$species

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
#cr_cds <- cds[,(colData(cds)$Idents.BRAIN.%in%tissues)]
#cds <- cr_cds
#
#check batch
pdf("00_batch.pdf")
plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE,cell_size = 1)
dev.off()

pdf("01_UMAP.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "Idents.BRAIN.",group_label_size = 5,cell_size = 1, label_cell_groups=FALSE)+scale_colour_discrete(type=mypalette)
dev.off()

#clustering
#cds <- cluster_cells(cds)
#pdf("02_clustering.pdf")
#plot_cells(cds,color_cells_by="cluster")
#dev.off()
#clusters_factors <- clusters(cds,"UMAP")

cdso <- cds
BRAINo <- BRAIN
for(spec in unique(BRAINo$species))
   {dir.create(spec)
    setwd(spec)
    BRAIN <- subset(x=BRAINo,subset=species==spec)
    cds <- cdso[,(colData(cdso)$species==spec)]
    #learn graph
    cds <- cluster_cells(cds)
    cds <- learn_graph(cds)
    pdf("03_trajectorise.pdf", width = 7, height = 7)
    p1 <-plot_cells(cds,
           color_cells_by = "Idents.BRAIN.",
           label_groups_by_cluster=TRUE,
           cell_size = 1,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_cell_groups=FALSE)+scale_colour_discrete(type=mypalette)
    print(p1)
    dev.off()
    get_earliest_principal_node1 <- function(cds, time_bin=opt$node1)
       {cell_ids <- which(colData(cds)[,"Idents.BRAIN."] == time_bin)
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
        root_pr_nodes}
    get_earliest_principal_node2 <- function(cds, time_bin=opt$node2)
       {cell_ids <- which(colData(cds)[,"Idents.BRAIN."] == time_bin)
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
        root_pr_nodes}
    cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node1(cds))
    pdf("trajectorise_set_root.pdf", width = 7, height = 7)
    p2 <- plot_cells(cds,color_cells_by = "pseudotime",cell_size = 1,label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
    print(p2)
    dev.off()
    pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
    BRAIN@meta.data$pseudotime<-pseudotime
    write.csv(file="pseudotime.csv",BRAIN@meta.data[,'pseudotime',drop=F],quote=F)
    BRAIN$pseudotime_phase <- pseudotime%/%2
    colData(cds)$pseudotime_phase <- pseudotime%/%2
    library(gridExtra)
    p<-list()
    for(i in names(table(BRAIN$group)))
       {imagename<-paste0("image_",gsub("-",".",i))
        temp<-cbind(BRAIN@meta.data[BRAIN@meta.data$group==i,c('group','pseudotime')],BRAIN@images[[imagename]]@coordinates[,1:2])
        temp2<-cbind(BRAIN0@meta.data[BRAIN0@meta.data$group==i,'group',drop=F],BRAIN0@images[[imagename]]@coordinates[,1:2])
        temp2$cell<-rownames(temp2)
        temp$cell<-rownames(temp)
        temp3 <- merge(temp2,temp,all.x=T)
        p[[i]]<-ggplot(temp3) + geom_point(mapping=aes(x=x, y=y, colour = pseudotime),size = 3)+
        scale_color_viridis(option="plasma",end=max(temp$pseudotime[!is.infinite(temp$pseudotime)])/max(BRAIN$pseudotime[!is.infinite(BRAIN$pseudotime)]),na.value="grey")+
        theme_bw() + theme(panel.grid=element_blank())+
        ggtitle(i)+coord_fixed()}
    sum<-grid.arrange(grobs=p,ncol=4)
    ggsave2("pseudotime_spatial.pdf",sum,width=25,height=20)
    #finding genes changing along pseudotime
    #The data frame pr_graph_test_res has the Moran's I test results for each gene in the cell_data_set. 
    #If you'd like to rank the genes by effect size, sort this table by the morans_Icolumn, which ranges from -1 to +1. 
    #A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells have very similar values of a gene's expression. 
    #Significant values much less than zero are generally rare.
    ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph=opt$graph, cores=4)
    pr_deg <- subset(ciliated_cds_pr_test_res, q_value < 1e-6)
    pr_deg_ids<-rownames(pr_deg)
    pr_deg$id <- rownames(pr_deg)
    write.csv(file="pr_deg.csv",pr_deg,sep=",",quote=FALSE)
    
    #collect the trajectory-variable genes into modules:
    gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-3)
    write.csv(gene_module_df, "gene_module_df.csv",col.names = F,quote=F)
    results1 <- merge(gene_module_df,pr_deg,by="id",all.x=TRUE)
    results <- merge(results1,info,by="id",all.x=TRUE)
    name3<-paste0("phase_module",".csv")
    write.csv(results, name3,quote=F)
    plot <- plot_cells(cds, genes=gene_module_df%>%filter(module=="1"), show_trajectory_graph=TRUE,label_leaves=FALSE)+coord_fixed()
    for (mol in 2:max(gene_module_df$module%>%as.numeric)) 
        {p <-  plot_cells(cds, genes=gene_module_df%>%filter(module==mol), show_trajectory_graph=TRUE,label_leaves=FALSE)+coord_fixed()
        plot <- plot+p}
    ggsave2(file="pr_deg_module.pdf",plot,width=10,height=10)   
    saveRDS(file="data.RDS",cds) 
    setwd("../")}




