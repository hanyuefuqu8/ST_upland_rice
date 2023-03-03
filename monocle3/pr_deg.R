library(dplyr)
library(VennDiagram)
library(ggplot2)
line<-c("GS69","GS106","GS74","GS112")

p<-list()
for (i in line)
{
path<-paste0(i,"/pr_deg.csv")
data<-read.csv(path,header=T)
print(dim(data))
data$group<-i
data$log_q<-0-log10(data$q_value)
data<-data[data$log_q>25,]
print(head(data))
p[[i]]<-data[,'X']
if(i!=line[1])
{
datas<-rbind(datas,data)
} else
{
datas<-data
}
}


datas<-datas[datas$log_q>25,]
tmp<-table(datas[,'X'])
tmp<-tmp[tmp>3]
datas$freq<-"not_repeatable"
datas[datas$X %in% names(tmp),'freq']<-"repeatable"
datas[datas$log_q >300,'log_q']=300

for (i in unique(datas$X))
{
av_log_q<-mean(datas[datas$X==i,'log_q'])
datas[datas$X==i,'log_q']<-av_log_q
}
temp<-unique(datas[,c('X','log_q','freq')])

pdf("q_val_pr_deg.pdf",width=3,height=3)
ggplot(temp, aes(x = log_q, fill = freq)) +
  geom_histogram(position = "identity",alpha=0.4)
dev.off()
write.table(file="all_pr_deg.txt",temp,sep='\t')

pdf("Venn_pr_deg.pdf")
# 四维韦恩图
venn.plot <- venn.diagram(
  x = list(GS69=p[['GS69']],GS106=p[['GS106']],GS74=p[['GS74']],GS112=p[['GS112']]),
  filename = NULL,
  col = "black",
  lty = "dotted",
  lwd = 4,
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 2.5,
  cat.fontfamily = "serif"
);
grid.draw(venn.plot);
grid.newpage();
dev.off()

