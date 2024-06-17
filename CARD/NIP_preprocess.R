library(Seurat)
library(dplyr)

refNip<-Read10X(data.dir = ".",strip.suffix = TRUE)
dataNip<- CreateSeuratObject(counts = refNip, project = "Japonica_root_tip", min.cells = 20, min.features = 200)
dataNip <- NormalizeData(dataNip, normalization.method = "LogNormalize", scale.factor = 10000)
dataNip <- FindVariableFeatures(dataNip, selection.method = "vst", nfeatures = 2000)
dataNip <- ScaleData(dataNip)
#dataNip <- RunPCA(dataNip, features = VariableFeatures(object = dataNip))
#dataNip <- FindNeighbors(dataNip, dims = 1:20)
#dataNip <- FindClusters(dataNip, resolution = 0.3)
data<-read.csv("Nip.expression.csv.gz",sep=",",header=T)
data<-data[,1:2]
dataNip$X <- rownames(dataNip@meta.data)
temp<-left_join(dataNip@meta.data,data,by="X")
dataNip$cluster<-temp$cluster
Idents(dataNip)<-dataNip$cluster
#saveRDS(file="Nip.RDS",dataNip)
save(scmat,sc_meta,file="Nip.RData")


scmat <- GetAssayData(object = dataNip, assay = "RNA", slot="counts")
sc_meta <- dataNip@meta.data[,c('cluster','orig.ident')]

data <- readRDS("../data.RDS")
names(data@images) <- gsub("-",".", names(data@images))
data.list <- SplitObject(data, split.by = "group")

for(i in 1:length(data.list))
    {
     name <- gsub("-",".", data.list[[i]]$group[1])
     image_name <- paste0("image_",name)
     spatial_count <- GetAssayData(object = data.list[[i]], assay = "Spatial", slot="counts")
     spatial_location <- data.list[[i]]@images[[image_name]]@coordinates[,c('x','y')]
     dir.create(name)
     save(spatial_count,spatial_location,file=paste0(name,"/spatial.RData"))
     }
