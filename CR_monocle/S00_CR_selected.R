library(Seurat)
library(ggplot2)
library(cowplot)
library(monocle3)
library(dplyr)
library(readr)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(purrr)
datasets <- readRDS("../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")
rdsf <- read.table("../info2",sep="\t",header=T)

tissues1 <- c("Pericycle_like","transition_tissue_1","transition_tissue_2","transition_tissue_3","Crown_root_primordia","Ground_tissue_2")
tissues2 <- c("Parenchyma_2","Xylem","Phloem","Parenchyma_1")
Steles <- read.csv("Stele_select.csv",header=T,row.names=1,quote="")
steles <- rownames(Steles)

single_neighbor <- function(datasets,cell_id)
    {cor <- strsplit(cell_id,"_")
     X <- gsub("X","",cor[[1]][2]) %>% as.numeric()
     Y <- cor[[1]][3] %>% as.numeric()
     neighbors <- c(paste0(cor[[1]][1],"_X",X-1,"_",Y),
                   paste0(cor[[1]][1],"_X",X+1,"_",Y),
                   paste0(cor[[1]][1],"_X",X,"_",Y-1),
                   paste0(cor[[1]][1],"_X",X,"_",Y+1))
     surroundings <- c(paste0(cor[[1]][1],"_X",X-1,"_",Y+1),
                   paste0(cor[[1]][1],"_X",X+1,"_",Y+1),
                   paste0(cor[[1]][1],"_X",X-1,"_",Y-1),
                   paste0(cor[[1]][1],"_X",X+1,"_",Y-1))
     #removve non_existed spot
     neighbors <- neighbors[!is.na(Idents(datasets)[neighbors])]
     surroundings <- surroundings[!is.na(Idents(datasets)[surroundings])]
     return(c(neighbors,surroundings))
    }

CR_split1 <- function(datasets,tissues1)
    {c0 <- list()
     n=1
     temp <- c()
     total_cells <- WhichCells(datasets,ident=tissues1[c(5,6)])
     rest_cell <- total_cells
     while(n==1||!length(rest_cell)==0)
         {temp0 <- rest_cell[1]
          temp1 <- rest_cell[1]
          while(length(temp0)!=0)
             {temp00 <- c()
              for(cell_id in temp0)
                 {neighbors <- single_neighbor(datasets,cell_id)
                  sel <- neighbors[Idents(datasets)[neighbors]%in%tissues1[c(5,6)]]
                  sel <- sel[!sel %in% temp1]
                  temp00 <- c(temp00,sel)}
              temp00 <- unique(temp00)
              temp0 <- temp00
              temp1 <- c(temp1,temp0)}
          c0[[n]]<-temp1
          n=n+1
          rest_cell <- total_cells[!total_cells %in% unlist(c0)]}
      return(c0)}

CR_split2 <- function(datasets,c0,tissues1,tissues2)
     {c1 <- list()
      n=1
      for(i in 1:length(c0))
         {if(length(c0[[i]])>20)
             {c1[[n]] <- c0[[i]]
              temp00 <- c()
              for(cell_id in c0[[i]])
                  {neighbors <- single_neighbor(datasets,cell_id)
                   sel <- neighbors[Idents(datasets)[neighbors]%in%tissues1[-c(5,6)]]
                   temp00 <- c(temp00,sel)}
              temp00 <- unique(temp00)
              temp1 <- temp00
              temp0 <- temp00
              terminate=0
              while(terminate<3&&length(temp00)!=0)
                 {terminate=0
                  temp00 <- c()
                  for(cell_id in temp0)
                     {neighbors <- single_neighbor(datasets,cell_id)
                      if(Idents(datasets)[cell_id]==tissues1[1] && (1 %in% (Idents(datasets)[neighbors]%in%tissues2)))
                       {terminate <- terminate+1}
                      sel <- neighbors[Idents(datasets)[neighbors]%in%tissues1[-c(5,6)]]
                      sel <- sel[!sel %in% temp1]
                      temp00 <- c(temp00,sel)}
                   temp00 <- unique(temp00)
                   temp0 <- temp00
                   temp1 <- c(temp1,temp0)}
               c1[[n]] <- c(c1[[n]],temp1)
               n=n+1}}
      rep <- table(unlist(c1))[table(unlist(c1))>1] %>% names()
      for (i in 1:length(c1))
          {c1[[i]]<-c1[[i]][!c1[[i]]%in%rep]}
      return(c1)}

CR_split3 <- function(datasets,c1,tissues1,iter)
    {c2 <- list()
     for(i in 1:length(c1))
        {c0 <- c1[[i]][Idents(datasets)[c1[[i]]]!=tissues1[1]]
         temp0 <- c0
         temp1 <- c0
         m=1
         while(m<=iter)
           {temp00 <- c()
            for(cell_id in temp0)
               {neighbors <- single_neighbor(datasets,cell_id)
                sel <- neighbors[Idents(datasets)[neighbors]%in%tissues1[1]]
                sel <- sel[!sel %in% temp1]
                temp00 <- c(temp00,sel)}
            temp00 <- unique(temp00)
            temp0 <- temp00
            temp1 <- c(temp1,temp0)
            m=m+1}
         c2[[i]] <- temp1}
      rep <- table(unlist(c2))[table(unlist(c2))>1] %>% names()
      for (i in 1:length(c2))
          {c2[[i]]<-c2[[i]][!c2[[i]]%in%rep]}
      return(c2)}


Pericycle_neighbor <- function(datasets,c1,stele_cells)
    {c0 <- c()
     c2 <- unlist(c1)[Idents(datasets)[unlist(c1)]!="Pericycle_like"]
     for(cell_id in unlist(c2))
        {neighbors <- single_neighbor(datasets,cell_id)
         if (sum(neighbors %in% stele_cells)>=1)
             {c0 <- c(c0,cell_id)}}
    return(c0)}

plot_root_shape <- function(datasets,datas, imagename, cellss,tissues)
    {exp<-data.frame(cell.type=Idents(datas),cellid=names(Idents(datas)))
    temp0 <- datasets@images[[imagename]]@coordinates[,1:2]
    temp0$cellid <- rownames(temp0)
    temp<-left_join(temp0,exp,by = 'cellid')
    temp$cell.type<- ordered(temp$cell.type, levels = tissues)
    mypalette <- viridis(length(unique(Idents(datas))))
    names(mypalette) <- unique(Idents(datas))
    mypalette <- c(mypalette ,na.value  = "grey")
    temp %>%
    group_by(cell.type) %>%
    do(model = kmeans(.[c('x', 'y')], 1)) %>% ### kmeans 计算一个中心点
    ungroup() %>% group_by(cell.type) %>%
    do(map_df(.$model, broom::tidy)) %>% ### 整理模型数据
    ungroup() %>% select(cell.type,x,y ) %>% data.frame() %>%
    dplyr::rename(x.center=x,y.center=y,cell.type=cell.type) -> label.data
    plot2<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = cell.type),size = 1)+
        scale_colour_manual(values = mypalette)+
        ggtitle(imagename)+
        theme(panel.grid=element_blank())+
        #geom_label(data = label.data, aes(label = cell.type,x = x.center,y = y.center))+
        coord_fixed()
    #ggsave2("root_shape.pdf",plot2,width=4,height=7)
    return(plot2)
    }


c0 <- CR_split1(datasets=datasets,tissues1=tissues1)
datas <- subset(datasets,cells=unlist(c0))
p<-list()
for(imagename in names(datasets@images))
    {p[[imagename]]<- plot_root_shape(datasets=datasets, datas=datas,imagename=imagename, cellss=unlist(c0),tissues=tissues1)}
sum<-grid.arrange(grobs=p,ncol=7)
ggsave2("selected_tissue_dimspatial1.pdf",sum,width=40,height=20)

c1 <- CR_split2(datasets=datasets,c0=c0,tissues1=tissues1,tissues2=tissues2)
datas <- subset(datasets,cells=unlist(c1))
p<-list()
for(imagename in names(datasets@images))
    {p[[imagename]]<- plot_root_shape(datasets=datasets, datas=datas,imagename=imagename, cellss=unlist(c1),tissues=tissues1)}
sum<-grid.arrange(grobs=p,ncol=7)
ggsave2("selected_tissue_dimspatial2.pdf",sum,width=40,height=20)


c2 <- CR_split3(datasets=datasets,c1=c1,tissues1=tissues1,iter=3)
datas <- subset(datasets,cells=unlist(c2))
p<-list()
for(imagename in names(datasets@images))
    {p[[imagename]]<- plot_root_shape(datasets=datasets, datas=datas,imagename=imagename, cellss=unlist(c2),tissues=tissues1)}
sum<-grid.arrange(grobs=p,ncol=7)
ggsave2("selected_tissue_dimspatial3.pdf",sum,width=40,height=20)

#remove little vascular pericycle
for(i in 1:length(c2))
    {c0 <- c2[[i]]
     peri <- c0[Idents(datasets)[c0]=="Pericycle_like"]
     lv <- peri[!peri%in%steles]
    c2[[i]] <- c0[!c0%in%lv]}

datas <- subset(datasets,cells=unlist(c2))
datas@meta.data$cell.type<-as.character(datas@meta.data$cell.type)
for (i in 1:length(c2))
    {datas@meta.data[rownames(datas@meta.data)%in%c2[[i]],'cell.type'] <- i}
Idents(datas)<-datas$cell.type
p<-list()
for(imagename in names(datasets@images))
    {p[[imagename]]<- plot_root_shape(datasets=datasets, datas=datas,imagename=imagename, cellss=unlist(c2),tissues=c(1:length(c2)))}
sum<-grid.arrange(grobs=p,ncol=7)
ggsave2("selected_CR.pdf",sum,width=40,height=20)

c0 <- Pericycle_neighbor(datasets=datasets,c1=c2,stele_cells=steles)
datass <- subset(datasets,cells=c0)
p<-list()
for(imagename in names(datasets@images))
    {p[[imagename]]<- plot_root_shape(datasets=datasets, datas=datass,imagename=imagename, cellss=c0,tissues=tissues1)}
sum<-grid.arrange(grobs=p,ncol=7)
ggsave2("Pericycle_neighbor_dimspatial.pdf",sum,width=40,height=20)
write.csv(file="Pericycle_neighbor.csv", datas@meta.data[c0,'cell.type',drop=FALSE],quote=F)
write.csv(file="CR_selected.csv",datas@meta.data[,'cell.type',drop=FALSE],quote=F)
