library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
dataset <- readRDS("../data2.RDS")

Change_cell_type <- function(cell_id)
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
    neighbors <- neighbors[!is.na(Idents(dataset)[neighbors])]
    surroundings <- surroundings[!is.na(Idents(dataset)[surroundings])]
    dis1 <- table(Idents(dataset)[neighbors])
    dis2 <- table(Idents(dataset)[surroundings])
    dis3 <- dis1+dis2
    cell.type <- Idents(dataset)[cell_id] %>% as.character()
    if(!(Idents(dataset)[cell_id] %in% Idents(dataset)[surroundings]||Idents(dataset)[cell_id] %in% Idents(dataset)[neighbors]))
    {if(dis3[dis3!=0]%>%length()==1)
        {cell.type <- names(dis3[dis3!=0])
        }else if (dis3[dis3!=0]%>%length()==2){
        if(dis1[dis1!=0]%>%length()==1)
             {if(1 %in% dis2)
                  {cell.type <- Idents(dataset)[neighbors][1] %>% as.character()}
             } else {
             cell.type <- "Unknown"
             }
        }else{
        cell.type <- "Unknown"}
    }
    return(cell.type)
    }


for (cell_id in Cells(dataset))
    {
#     if (!cell_id %in% WhichCells(dataset,ident=c("Stele_1","Endodermis")))
#     {
     dataset@meta.data[cell_id,'cell.type'] <- Change_cell_type(cell_id)
#     }
    }

Idents(dataset) <- dataset$cell.type

saveRDS(dataset, "data_pos_adj.RDS")




