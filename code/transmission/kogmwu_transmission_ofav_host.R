#### packages ####

# install.packages("KOGMWU")
library(KOGMWU)


#### pairwise fates (lpv) ####

# loading KOG annotations
gene2kog=read.table("Ofaveolata_iso2kogClass.tab",sep="\t", fill=T)
head(gene2kog)

diseased_healthy=load('diseased_healthy_lpv.RData')
diseased_healthy # names of datasets in the package
lpv.diseased_healthy=kog.mwu(diseased_healthy.p,gene2kog) 
lpv.diseased_healthy 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased_healthy"=lpv.diseased_healthy,"diseased_healthy"=lpv.diseased_healthy))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_transmission_ofav_host_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()


#### pairwise fates (fc) ####

diseased_healthy=load('diseased_healthy_fc.RData')
diseased_healthy # names of datasets in the package
fc.diseased_healthy=kog.mwu(diseased_healthy.fc,gene2kog) 
fc.diseased_healthy 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased_healthy"=fc.diseased_healthy,"diseased_healthy"=fc.diseased_healthy))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_transmission_ofav_host_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()