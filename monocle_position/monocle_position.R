library(Seurat)
library(ggplot2)
library(cowplot)
library(monocle3)
library(dplyr)
library(readr)

datasets <- readRDS("../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")

tissues <- c("Meristem_2","Endodermis","Cortex","Trasition_zone")

load("align.RData")

datasets <- AddMetaData(object = datasets, metadata = phases,col.name = 'phase')
datasets <- AddMetaData(object = datasets, metadata = lengths,col.name = 'length')

## choose the interested clusters
datasets<-subset(x=datasets,idents=tissues)

rdsf <- read_tsv("../info2")
temp<-left_join(datasets@meta.data,rdsf,by="group")
datasets@meta.data$species<-temp$species
datasets@meta.data$medium<-temp$medium
#dataset <- subset(x=datasets, medium=="Y")
dataset <- datasets
datasets.list <- SplitObject(dataset, split.by="species")

info <- read.table("../Markergenes/gene_v3.txt",quote="",header=T,sep="\t")
colnames(info) <- c("id","gene_name")
for(i in 1:length(datasets.list))
    {
    name1 <- names(datasets.list[i])
    dir.create(name1)
    setwd(name1)
    BRAINs <- datasets.list[[i]]
    BRAINs$cell.type <- Idents(BRAINs)
    meta <- BRAINs@meta.data
    attach(meta)
    len1 <- aggregate(length~cell.type+group,data=meta,FUN=min)
    len2 <- aggregate(length~cell.type+group,data=meta,FUN=max)
    len1$maxlength <- len2$length
    len1$range <- len1$maxlength-len1$length
    p <- ggplot() +
    geom_linerange(data = len1 %>% arrange(-range), 
                   mapping=aes(x = cell.type, ymin = length, ymax = maxlength, 
                             lwd = 1, color = range)) +
    scale_color_continuous(high = "#cffcd0", low = "#4fc290") +
    coord_flip() +
    theme_classic()    
    ggsave2("range.png",p)
    detach(meta)

    data <- as(as.matrix(BRAINs@assays$Spatial@counts), 'sparseMatrix')
    pd <- data.frame(BRAINs@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

    cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)

    #remove batch effect
    cds <- preprocess_cds(cds, num_dim = 100)
#    cds <- align_cds(cds, alignment_group = "group")
    cds <- reduce_dimension(cds)


    #Finding genes that change as a function of phase and treatment
    gene_fits <- fit_models(cds, model_formula_str = "~length")
    fit_coefs <- coefficient_table(gene_fits)
    fit_coefs_table <-  fit_coefs %>% filter(term != "(Intercept)") %>% filter (q_value < 0.05) %>%
      select(gene_id, term, q_value, estimate)

    pr_deg_ids <- unique(fit_coefs_table$gene_id)

    gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
    fit_coefs_table$id <- fit_coefs_table$gene_id
    results1 <- merge(gene_module_df,fit_coefs_table,by="id",all.x=TRUE)
    results <- merge(results1,info,by="id",all.x=TRUE)
    name3<-paste0("phase_module",".csv")
    write.csv(results, name3,quote=F)

    #outputheatmap
	#colData(cds)$phase<- ordered(colData(cds)$phase, levels = c("2","1","5"))
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),cell_group=colData(cds)$phase)
    agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
    name2<-paste0(name1,"_phase_module",".pdf")
    require(ggplotify)
    p2<-pheatmap::pheatmap(agg_mat,
                   scale="row", clustering_method="ward.D2",cluster_cols= FALSE)
    p1 <- as.ggplot(p2)
    p0<-cowplot::plot_grid(p, p1, ncol=1,rel_heights=c(1,3))
    ggsave2(file=name2,p0,device='pdf',width=6,height=6)
    setwd("../")
    }
