library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
                "outdir","o","1","character","outdir",
                "sample","s","1","character","sample,default=Maize",
                "tissue","t","1","character","tissue,default=Embro",
                "kfilt","k","1","integer","k.filter for merge",
                "dims","d","1","integer","dims option for FindNeighbors,default=15",
                "resolution","r","1","numeric","resolution option for FindClusters,[0.4-1.2],default=1",
                "help","h","0","logical", "Usage: Rscript runiDrop.R -i <input> -o <outdir> [-s SAMPLE -t TISSUE]",
                "minCG","m","1","integer","minimum gene number for each cell, i.e. nFeature_RNA, default=200",
                "rds","f","0","logical", "Save the RDS file"
               #"list","l","1","character", "interested gene"
                ),byrow=T,ncol=5)
opt = getopt(arg)
if(!is.null(opt$help) || is.null(opt$input)){
        cat(paste(getopt(arg, usage = T), "\n"))
 q()
}
if (is.null(opt$sample)){
        opt$sample <- "Maize"
}
if (is.null(opt$tissue)){
        opt$tissue <- "Embro"
}
if (is.null(opt$outdir)){
        opt$outdir <- "output"
}
if (is.null(opt$dims)){
        opt$dims <- 15
}
if (is.null(opt$resolution)){
        opt$resolution <- 1
}
if (is.null(opt$minCG)){
 opt$minCG <- 200
}


library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
file=scan(opt$input,"")
dir.create(opt$outdir);
setwd(opt$outdir);
BRAIN<-list()
n=0
w=c()
h=c()
id=c()
vf<-c()

for(i in file)
{
n=n+1
brain.data <- read.csv(i, header = T, sep=",")
rownames(brain.data) <- brain.data[,1]
brain.data <- brain.data[,2:dim(brain.data)[2]]
meta <- data.frame(Barcodes = colnames(brain.data), Sample = opt$sample, Tissue = opt$tissue, Annotation = NA, Celltype = NA)
BRAIN[[n]] <- CreateSeuratObject(counts = brain.data,
                              project = opt$tissue,
                              min.cells = 0,
                              min.features = 0,
                              assay = 'Spatial')

#adding location information
locations<-data.frame(substring(colnames(brain.data),2))
locations[,1:2]<-str_split_fixed(locations$substring.colnames.brain.data...2.,"_",2)
names(locations)<-c("x","y")
rownames(x = locations) <-colnames(brain.data)
locations$x<-as.numeric(locations$x)
locations$y<-as.numeric(locations$y)
w[n]<-as.numeric((max(locations$x) - min(locations$x) + 1)/10)
h[n]<-as.numeric((max(locations$y) - min(locations$y) + 1)/10)
BRAIN[[n]][['image']] <- new(Class = 'SlideSeq',assay = "Spatial",
    coordinates = locations)
    BRAIN[[n]] <- subset(BRAIN[[n]], subset = nFeature_Spatial > opt$minCG)
    name<-strsplit(i,"/",fixed=T)[[1]][2]
    BRAIN[[n]]@meta.data$group<-gsub(".csv","",name)
    id[n]<-gsub(".csv","",name)
    BRAIN[[n]] <- SCTransform(BRAIN[[n]], assay = "Spatial", verbose = FALSE) 
    #BRAIN[[n]] <- FindVariableFeatures(BRAIN[[n]], selection.method = "vst", nfeatures = 3000)
    vf<-c(vf,VariableFeatures(BRAIN[[n]]))  
}

#ifnb.list <- BRAIN
#features <- SelectIntegrationFeatures(object.list = ifnb.list)
#ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
#pancreas.anchors <- FindIntegrationAnchors(object.list =  ifnb.list,normalization.method = "SCT",anchor.features = features, dims = 1:15,k.anchor = 5,k.filter = opt$kfilt)
#pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors,normalization.method = "SCT",dims = 1:15)

pancreas.integrated<-merge(BRAIN[[1]], y=BRAIN[-1],add.cell.ids = id)
VariableFeatures(pancreas.integrated) <- vf
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures( pancreas.integrated))

pdf("05_ElbowPlot.pdf", width = 10)
ElbowPlot(object = pancreas.integrated)
dev.off()

pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:opt$dims)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = opt$resolution)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:opt$dims)

library(RColorBrewer)
mypalette<-c(brewer.pal(11,"BrBG"),brewer.pal(11,"PiYG"),brewer.pal(11,"PRGn"))


pdf("09_umap.pdf", width = 40,height=20)
plot1 <- DimPlot(object = pancreas.integrated, reduction = "umap", label = T,pt.size=1,group.by = "ident",cols=mypalette)
plot2 <- DimPlot(object = pancreas.integrated, reduction = "umap", pt.size=1,group.by = "group",split.by = "group")
plot3 <- DimPlot(object = pancreas.integrated, reduction = "umap", pt.size=1,group.by = "group")
wrap_plots(plot1+plot2+plot3)
dev.off()

#change the first name of images
b<-paste0("image_",names(table(pancreas.integrated$group)))
c<- names(pancreas.integrated@images)
names(pancreas.integrated@images)[1]<-b[!b%in%c]

#展示抠出来的图
library(gridExtra)
p<-list()
for(i in names(pancreas.integrated@images))
{
p[[i]]<- SpatialDimPlot(pancreas.integrated,images=i,stroke=0,pt.size=12,label=T,label.size=3)+scale_fill_manual(breaks=levels(Idents(pancreas.integrated)),values=mypalette)+ggtitle(i)
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2("dimspatial.png",sum,dpi=600,width=20,height=20)



BRAIN.list <- SplitObject(pancreas.integrated, split.by = "group")
for(i in 1:length(file))
{
   BRAIN.list[[i]]@images<- BRAIN.list[[i]]@images[i]
   pdf(paste0(BRAIN.list[[i]]@meta.data$group,"_dimspatial.pdf"))
   plot1 <- SpatialDimPlot(BRAIN.list[[i]],stroke=0,pt.size=8,label=T,label.size=3)+scale_fill_manual(breaks=levels(Idents(pancreas.integrated)),values=mypalette)
   print(plot1)
   dev.off()

   pdf(paste0(BRAIN.list[[i]]@meta.data$group,"_dimsplit.pdf"), width = 10,height=10)
   plot1 <-SpatialDimPlot(BRAIN.list[[i]], cells.highlight = CellsByIdentities(object = BRAIN.list[[i]], idents = levels(Idents(BRAIN.list[[i]]))), facet.highlight = TRUE, ncol = 5, stroke = 0, pt.size=1)
   print(plot1)
   dev.off() 
}


saveRDS(pancreas.integrated,"data.RDS")













