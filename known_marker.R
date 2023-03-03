library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
                "csv","v","1","character","file containing expression matrix csv files list",
                "genefile","g","1","character","gene_tissue table"
                ),byrow=T,ncol=5)
opt = getopt(arg)

library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Giotto)
library(monocle3)
library(viridis)

library(Giotto)
my_instructions = createGiottoInstructions(python_path = '/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/stlearn/bin/python')
BRAIN<-readRDS(opt$input)
BRAIN$cell_type<-Idents(BRAIN)
markergene<-read.table(opt$genefile,header=T,sep="\t",quote="")
colnames(markergene)<-c("gene","Cell_type","genename")
markergene<-markergene[markergene$gene!="",]
file<-scan(opt$csv,"")


#dotplot
for (i in names(table(markergene$Cell_type)))
{
temp<-markergene[markergene$Cell_type==i,]
temp<-temp[!duplicated(temp$gene),]
a<-temp[,'gene']%in%rownames(BRAIN[["SCT"]]@data)
temp<-temp[a,]
if(nrow(temp)!=0)
{
pdf(paste0(gsub(" ","_",i),"_dotplot.pdf"),height=20,width=20)
p1<-DotPlot(BRAIN, features = unique(temp[,'gene']),assay="SCT")+
           scale_x_discrete(label=temp[,'genename'])+
		   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
		   coord_flip()
print(p1)
dev.off()
}
}


#using giotto to plot
n=1
for(f in file)
{ 
name<-unlist(strsplit(f,'/'))
name1<-name[length(name)]
name3<-unlist(strsplit(name1,"\\."))[1]
data<-read.csv(f,sep=",")
rownames(data)<-data[,1]
data<-data[,-1]

locations<-data.frame(substring(names(data),2))
locations[,1:2]<-str_split_fixed(locations$substring.names.data...2.,"_",2)
names(locations)<-c("V1","V2")

locations$V1<-as.numeric(locations$V1)
locations$V2<-as.numeric(locations$V2)
h <- as.numeric((max(locations$V2) - min(locations$V2) + 1)/20)
w <- as.numeric((max(locations$V1) - min(locations$V1) + 1)/20)

my_giotto_object = createGiottoObject(raw_exprs = data,
                                      spatial_locs = locations,
                                      instructions = my_instructions
                                      )
    
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object, scalefactor = 6000, verbose = T)   
#鍙栫粏鑳炵被鍨嬬殑鍩哄洜澶т簬涓€涓殑
names<-table(markergene$Cell_type)[table(markergene$Cell_type)>1]%>%names()
genes<-markergene[markergene$Cell_type%in%names,]$Cell_type
names(genes)<- markergene[markergene$Cell_type%in%names,]$gene
my_giotto_object = createMetagenes(my_giotto_object, gene_clusters = genes, name = 'cluster_metagene')
    p<-spatCellPlot(my_giotto_object,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = names,
             point_size = 1.5,point_border_stroke=0, cow_n_col = 5)   
ggsave2(paste0( name3,"_1_merge.pdf"),p,device="pdf",width=w*5,height=h*5)

name2<-table(markergene$Cell_type)[table(markergene$Cell_type)<2]%>%names()
gene2<-markergene[markergene$Cell_type%in%name2,]$gene
    if (sum( gene2%in%my_giotto_object@gene_ID))
        {
p2<-spatGenePlot(my_giotto_object,
                 genes = unique(gene2), point_size = 3,
                 cow_n_col = 5
                 )
ggsave2(paste0( name3,"_2_merge.pdf"),p2,device="pdf",width=w*5,height=h*1.5) 
        }
    n=n+1       
    }
	
#using monocle3 to plot
data <- as(as.matrix(BRAIN@assays$Spatial@counts), 'sparseMatrix')
pd <- data.frame(BRAIN@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell_type)
gene_group_df<-markergene[markergene$gene!="",]

agg_mat <- aggregate_gene_expression(cds, gene_group_df, cell_group_df)

#pdf("known_marker_heatmap.pdf")
plot1<-pheatmap::pheatmap(agg_mat,scale="row",color=rocket(7) ,clustering_method="ward.D2")
#plot1
#dev.off()
ggsave2("known_marker_heatmap.pdf",plot1,width=4,height=4)



