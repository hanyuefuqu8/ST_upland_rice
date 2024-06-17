library(dplyr)
library(Seurat)
dataset <- readRDS("../diff_gene_SCT_position_adjusting/data_pos_adj.RDS")
info <- read.table("../info2",header=T)

species <- unique(info$species)
tissues <- unique(Idents(dataset))
tissues <- tissues[!is.na(tissues)]


all_tissue <- combn(tissues,2)
result <- list()
n=1
for(j in 1:length(species))
    {for (i in 1:ncol(all_tissue))
         {tar_ti <- all_tissue[,i]
          temp1 <- try(read.csv(paste0(species[j],"/",tar_ti[1],"_",tar_ti[2],"_markers.csv"),header=T,quote="",row.names=1),silent=TRUE)
          temp2 <- try(read.csv(paste0(species[j],"/",tar_ti[2],"_",tar_ti[1],"_markers.csv"),header=T,quote="",row.names=1),silent=TRUE)
         if(!'try-error' %in% class(temp1)||!'try-error' %in% class(temp2))
         {if(!'try-error' %in% class(temp1))
              {temp <- temp1
               temp <- temp[temp$p_val<0.05,]
               if(nrow(temp)>0)
                    {temp$change <- "N"
                     if(nrow(temp[temp$avg_log2FC>0,])>0)
                           {temp[temp$avg_log2FC>0,]$change <- "Down"}
                     if(nrow(temp[temp$avg_log2FC<0,])>0)
                           {temp[temp$avg_log2FC<0,]$change <- "Up"}}}
          else{temp <- temp2
              temp <- temp[temp$p_val<0.05,]
               if(nrow(temp)>0)
                    {temp$change <- "N"
                     if(nrow(temp[temp$avg_log2FC>0,])>0)
                           {temp[temp$avg_log2FC>0,]$change <- "Up"}
                     if(nrow(temp[temp$avg_log2FC<0,])>0)
                           {temp[temp$avg_log2FC<0,]$change <- "Down"}}}
          if(nrow(temp)>0){
          temp$tissue <- paste0(tar_ti[1],"_",tar_ti[2])
          temp$species <- species[j]
          temp <- temp[,c('gene','name','tissue','species','change')]
          result[[n]] <- temp
          n=n+1}}}}

all <- rbind(result[[1]],result[[2]])
for (i in 3:length(result))
    {all <- rbind(all,result[[i]])}

stat <- data.frame(tissue=c(),
                  upland_rice_up=c(),
                  irrigated_rice_up=c(),
                  intersect_up=c(),
                  upland_rice_down=c(),
                  irrigated_rice_down=c(),
                  intersect_down=c())
for (i in 1:ncol(all_tissue))
    {tar_ti <- all_tissue[,i]
    temp <- all[all$tissue==paste0(tar_ti[1],"_",tar_ti[2]),]
    gene1 <- temp[temp$species=="Irrigated_rice"&temp$change=="Up",'gene']
    gene2 <- temp[temp$species=="Upland_rice"&temp$change=="Up",'gene']
    genes <- intersect(gene1,gene2)
    gene11 <- temp[temp$species=="Irrigated_rice"&temp$change=="Down",'gene']
    gene22 <- temp[temp$species=="Upland_rice"&temp$change=="Down",'gene']
    geness <- intersect(gene11,gene22)
    stat[i,'tissue'] <- paste0(tar_ti[1],"_",tar_ti[2])
    stat[i,'upland_rice_up'] <- length(setdiff(gene2,gene1))
    stat[i,'irrigated_rice_up'] <- length(setdiff(gene1,gene2))
    stat[i,'intersect_up'] <- length(genes)
    stat[i,'upland_rice_down'] <- length(setdiff(gene22,gene11))
    stat[i,'irrigated_rice_down'] <- length(setdiff(gene11,gene22))
    stat[i,'intersect_down'] <- length(geness)
    temp[temp$gene %in% c(genes,geness),'species']="both"
    temp <- unique(temp)
    write.csv(file=paste0(tar_ti[1],"_",tar_ti[2],"all_gene_.csv"),temp,quote=F,row.names=F)
    all <- all[!all$tissue==paste0(tar_ti[1],"_",tar_ti[2]),]
    all <- rbind(all,temp)
    }  

write.csv(file="stata.csv",stat,quote=F,row.names=F)
write.csv(file="all_gene.csv",all,quote=F,row.names=F)


sink("gene_analysis_result.txt", append = TRUE)
print(ncol(all))
print(paste0("union:", length(all$gene %>% unique())))
print(paste0("intersect:",length(all[all$species=="both",]$gene %>% unique())))
print(paste0("Irrigated_rice:",length(all[all$species=="Irrigated_rice",]$gene %>% unique())))
print(paste0("Upland_rice:",length(all[all$species=="Upland_rice",]$gene %>% unique())))
sink()






