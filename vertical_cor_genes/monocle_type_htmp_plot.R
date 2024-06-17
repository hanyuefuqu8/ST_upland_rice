library(Seurat)
library(ggplot2)
library(cowplot)
library(monocle3)
library(dplyr)
library(readr)
library(pheatmap)


datasets <- readRDS("../../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")

tissues <- c("Meristem_2","Endodermis","Cortex","Trasition_zone","Meristem_1","Stele_1","Stele_2","Stele_3")

load("../../monocle_position/align.RData")
datasets <- AddMetaData(object = datasets, metadata = phases,col.name = 'phase')
datasets <- AddMetaData(object = datasets, metadata = lengths,col.name = 'length')
datasets$phase2 <- datasets$length%/%1

## choose the interested clusters
datasetso <- datasets
datasets<-subset(x=datasets,idents=tissues)
datasets$type <- "cortex"
datasets@meta.data[datasets@meta.data$cell.type %in% c("Meristem_1","Stele_1","Stele_2","Stele_3"),'type'] <- "stele"

rdsf <- read_tsv("../../info2")
temp<-left_join(datasets@meta.data,rdsf,by="group")
datasets@meta.data$species<-temp$species
datasets@meta.data$medium<-temp$medium
#dataset <- subset(x=datasets, medium=="Y")
datasets.list <- SplitObject(datasets, split.by="type")
datasets.list1 <- SplitObject(datasets, split.by="species")

got_mtx <- function(datasets,gene)
    {temp <- AverageExpression(datasets,assays='SCT',slot='counts',features=gene,group.by="phase2")
     temp <- temp$SCT
     temp <- temp[,colnames(temp) %>% as.numeric %>% order]
     return(temp)}

got_mtx2 <- function(datasets,gene)
    {temp <- AverageExpression(datasets,assays='SCT',slot='counts',features=gene,group.by="phase2")
     temp <- temp$SCT
     temp <- temp[,colnames(temp) %>% as.numeric %>% order]
     colnames(temp) <- paste0(datasets$species[1],"__",colnames(temp))
     return(temp)}

###############################################got vertical genes#####################################################
dev <- read.csv("02_vertical_core_dev_genes.csv",header=T)
#dev$gene <- dev$id
gene <- dev$gene
temp1 <- got_mtx(datasets.list[[1]],gene)

temp2 <- got_mtx(datasets.list[[2]],gene)

temp0 <- cbind(temp1,temp2)
temp0 <- temp0[,colnames(temp0) %>% as.numeric %>% order]

colnames(temp2) <- paste0(datasets.list[[2]]$type[1],"_",colnames(temp2))
colnames(temp1) <- paste0(datasets.list[[1]]$type[1],"_",colnames(temp1))

###########################################got selected genes #########################################################
sele <- read.csv("gene.list",header=T)
genes <- sele$geneid
temp7 <- got_mtx2(datasets.list1[[1]],genes)
temp8 <- got_mtx2(datasets.list1[[2]],genes)
###########################################got vertical only genes ####################################################
dev2 <- read.csv("02_vertical_only_genes2.csv",header=T)
gene1 <- dev2[dev2$only_type =="irri_pseudo_only"|dev2$only_type =="irri_spatial_only",'id']
gene2 <- dev2[dev2$only_type =="up_pseudo_only"|dev2$only_type =="up_spatial_only",'id']


temp3 <- got_mtx2(datasets.list1[[1]],gene1)
temp4 <- got_mtx2(datasets.list1[[2]],gene1)
temp5 <- got_mtx2(datasets.list1[[2]],gene2)
temp6 <- got_mtx2(datasets.list1[[1]],gene2)

#####these codes were modified from https://tool.biomooc.com/R_scripts/index.html

############### color bar from monocle::plot_pseudotime_heatmap ##############
# 这个 color bar 来自于monocle，画热图效果很好。该渐变色原理和细节请参考 http://www.biomooc.com/R/R-color.html#2_3
# fn1
table.ramp = function (n, mid = 0.5, sill = 0.5, base = 1, height = 1) {
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill/2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill/2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base/2)) + 1
  base.max <- round((n - 1) * (mid + base/2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}

# fn2
rgb.tables=function (n, red = c(0.75, 0.25, 1), green = c(0.5, 0.25, 1),  blue = c(0.25, 0.25, 1)) {
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}

# fn3
blue2green2red=function (n) {
  rgb.tables(n, red = c(0.8, 0.2, 1), green = c(0.5, 0.4, 0.8), blue = c(0.2, 0.2, 1))
}

bks <- seq(-3.1,3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
#hmcols
# view the color bar
#barplot(rep(1, length(hmcols)), col=hmcols, border = NA, space=0, axes=F)


length0 <- colnames(temp0)
colnames(temp0) <- make.unique(colnames(temp0))

annotate_col <- data.frame(
    length = colnames(temp0),
    Length = length0,
    row.names = 1)

library(stringr)
#annotate_col$type <- str_split_fixed(annotate_col$Length, "_", 2)[,1]
#annotate_col$Length <- str_split_fixed(annotate_col$Length, "_", 2)[,2] %>% as.numeric

testSpan2=function(temp,annotate_col, span=0.2,type=1){
  df4=apply(temp, 1, function(x){
    predict(loess(x ~ seq(1, length(x)), span=span))
  })
  df4=as.data.frame(t(df4))
  colnames(df4)=colnames(temp)
 
  ## heatmap again
  if (type==1)
    {p=pheatmap(df4 , border_color = NA, scale='row', 
             clustering_method='ward.D2',
             cluster_cols = F,
             show_colnames = F,
             #show_row_names = F,
             annotation_col = annotate_col,
             #annotation_row = annote_row,
             #annotation_colors = annote_color,
             #gaps_col = c(60),
             cutree_rows = 3,
             color=hmcols)} else{
      p=pheatmap(df4 , border_color = NA, scale='row', 
             clustering_method='ward.D2',
             cluster_cols = F,
             cluster_rows= F,
             show_colnames = F,
             #show_row_names = F,
             annotation_col = annotate_col,
             #annotation_row = annote_row,
             #annotation_colors = annote_color,
             #gaps_col = c(60),
             cutree_rows = 3,
             color=hmcols)}
  return(p)
}

p0 <- testSpan2(temp0, annotate_col, span=2, type=1)

annotate_col <- data.frame(
    length = colnames(temp1),
    Length = colnames(temp1),
    row.names = 1)
annotate_col$type <- str_split_fixed(annotate_col$Length, "_", 2)[,1]
annotate_col$Length <- str_split_fixed(annotate_col$Length, "_", 2)[,2] %>% as.numeric
nr=rownames(temp1)[p0$tree_row[["order"]]]

p1 <- testSpan2(temp1[nr,], annotate_col, span=2, type=2) 


annotate_col <- data.frame(
    length = colnames(temp2),
    Length = colnames(temp2),
    row.names = 1)
annotate_col$type <- str_split_fixed(annotate_col$Length, "_", 2)[,1]
annotate_col$Length <- str_split_fixed(annotate_col$Length, "_", 2)[,2] %>% as.numeric

p2 <- testSpan2(temp2[nr,], annotate_col, span=2, type=2) 
select_gene <- scan("selected_gene.list","")
#p2_2 <- p2 +rowAnnotation(link = anno_mark(at = which(rownames(p2) %in% selected_gene), 
#                                      labels = selected_gene, labels_gp = gpar(fontsize = 10))

require(ggplotify)
p1 <- as.ggplot(p1)
p2 <- as.ggplot(p2)
p0<-cowplot::plot_grid(p1, p2, ncol=2)
ggsave2("monocle_htmp_for_v_core.pdf",p0,width=20,height=3)

got_plot <- function(temp,type0)
     {annotate_col <- data.frame(length = colnames(temp),Length = colnames(temp),row.names = 1)
      annotate_col$type <- str_split_fixed(annotate_col$Length, "__", 2)[,1]
      annotate_col$Length <- str_split_fixed(annotate_col$Length, "__", 2)[,2] %>% as.numeric
      p <- testSpan2(temp, annotate_col, span=2, type=type0)
      return(p)}

p7 <- got_plot(temp7,2)
p8 <- got_plot(temp8,2)
p7 <- as.ggplot(p7);p8 <- as.ggplot(p8)
p0<-cowplot::plot_grid(p7, p8, ncol=2)
ggsave2("monocle_htmp_for_selected_genes.pdf",p0,width=10,height=10)


p3 <- got_plot(temp3,1)
nr=rownames(temp3)[p3$tree_row[["order"]]]
p4 <- got_plot(temp4[nr,],2)
p3 <- as.ggplot(p3);p4 <- as.ggplot(p4)
p0<-cowplot::plot_grid(p3, p4, ncol=2)
ggsave2("monocle_htmp_for_v_only_irri.pdf",p0,width=10,height=3)

p5 <- got_plot(temp5,1)
nr=rownames(temp5)[p5$tree_row[["order"]]]
p6 <- got_plot(temp6[nr,],2)
p5 <- as.ggplot(p5);p6 <- as.ggplot(p6)
p0<-cowplot::plot_grid(p5, p6, ncol=2)
ggsave2("monocle_htmp_for_v_only_up.pdf",p0,width=10,height=3)








