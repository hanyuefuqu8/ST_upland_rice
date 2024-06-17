library(CARD)
library(cowplot)

load(file="spatial.RData")
load(file="../Nip.RData")

CARD_obj = createCARDObject(
	sc_count = scmat,
	sc_meta = sc_meta,
	spatial_count = spatial_count,
	spatial_location = spatial_location,
	ct.varname = "cluster",
	ct.select = unique(sc_meta$cluster),
	sample.varname = "orig.ident",
	minCountGene = 100,
	minCountSpot = 5)           

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
    "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(
	proportion = CARD_obj@Proportion_CARD,
	spatial_location = CARD_obj@spatial_location, 
 	colors = colors, 
  	radius = NULL) ### You can choose radius = NULL or your own radius number
ggsave2(file="pie.pdf",p1+coord_fixed(),width=10,height=5)

## select the cell type that we are interested
ct.visualize = colnames(CARD_obj@Proportion_CARD)
## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
	proportion = CARD_obj@Proportion_CARD,        
	spatial_location = CARD_obj@spatial_location, 
	ct.visualize = ct.visualize,                 ### selected cell types to visualize
	colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
	NumCols = 4,                                 ### number of columns in the figure panel
        pointSize = 1.0)+coord_fixed()               ### point size in ggplot2 scatterplot  
ggsave2(file="cell_type.pdf",p2,width=10,height=5)

## visualize the spatial distribution of two cell types on the same plot
p3 = CARD.visualize.prop.2CT(
proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
spatial_location = CARD_obj@spatial_location,                      ### spatial location information
ct2.visualize = ct.visualize[c(1,3)],              ### two cell types you want to visualize
colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))+coord_fixed()      ### two color scales        

ggsave2(file="2_cell_types.pdf",p3,width=10,height=5)           

#Visualize the cell type proportion correlation
p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
ggsave2(file="correlation.pdf",p4)

#Refined spatial map
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 3000,ineibor = 10,exclude = NULL)
## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
	y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)
library(ggplot2)
p5 <- ggplot(location_imputation, 
       aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    legend.position="bottom",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))+
    coord_fixed()
ggsave2("refined_grid.pdf")
p6 <- CARD.visualize.prop(
	proportion = CARD_obj@refined_prop,                         
	spatial_location = location_imputation,            
	ct.visualize = ct.visualize,                    
	colors = c("lightblue","lightyellow","red"),    
	NumCols = 4,pointSize = 0.2 )+coord_fixed()                           
ggsave2(file="refined_cell_types.pdf",p6,width=10,height=5)




#Extension of CARD in a reference-free version: CARDfree
#markerlist <- read.table("../PCMDB_root_expe.txt",sep="\t",header=FALSE,quote="")
markerlist <- read.table("../WJW_root_marker",sep="\t",header=FALSE,quote="")
markerlist <- markerlist[markerlist$V2!="",]
markerList <- list()
for(i in unique(markerlist[,1]))
    {
    markers <- markerlist[markerlist$V1==i,'V2']
    if(length(markers)>=20) 
    markerList[[i]] <- markers
    }


CARDfree_obj = createCARDfreeObject(
	markerList = markerList,
	spatial_count = spatial_count,
	spatial_location = spatial_location,
	minCountGene = 100,
	minCountSpot =5) 
CARDfree_obj = CARD_refFree(CARDfree_obj)

ct.visualize=colnames(CARDfree_obj@Proportion_CARD)
p10<-CARD.visualize.prop(
        proportion = CARDfree_obj@Proportion_CARD,
        spatial_location = CARDfree_obj@spatial_location,
        ct.visualize = ct.visualize,                 ### selected cell types to visualize
        colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
        NumCols = 4,                                 ### number of columns in the figure panel
        pointSize = 1.0)+coord_fixed()               ### point size in ggplot2 scatterplot
ggsave2(file="cell_type_for_marker.pdf",p10,width=10,height=5)









