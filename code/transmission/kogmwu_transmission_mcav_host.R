#### packages ####

# install.packages("KOGMWU")
library(KOGMWU)


#### pairwise fates (lpv) ####

# loading KOG annotations
gene2kog=read.table("Mcavernosa2015_iso2kogClass.tab",sep="\t", fill=T)
head(gene2kog)

nai_healthy=load('nai_healthy_lpv.RData')
nai_healthy # names of datasets in the package
lpv.nai_healthy=kog.mwu(nai_healthy.p,gene2kog) 
lpv.nai_healthy 

diseased_healthy=load('diseased_healthy_lpv.RData')
diseased_healthy # names of datasets in the package
lpv.diseased_healthy=kog.mwu(diseased_healthy.p,gene2kog) 
lpv.diseased_healthy 

diseased_nai=load('diseased_nai_lpv.RData')
diseased_nai # names of datasets in the package
lpv.diseased_nai=kog.mwu(diseased_nai.p,gene2kog) 
lpv.diseased_nai

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased_healthy"=lpv.diseased_healthy,"nai_healthy"=lpv.nai_healthy,"diseased_nai"=lpv.diseased_nai))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_transmission_mcav_host_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_transmission_mcav_host_corr_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="diseased_healthy",y="nai_healthy",ktable)
corrPlot(x="diseased_healthy",y="diseased_nai",ktable)
corrPlot(x="nai_healthy",y="diseased_nai",ktable)
dev.off()


#### pairwise fates (fc) ####

nai_healthy=load('nai_healthy_fc.RData')
nai_healthy # names of datasets in the package
fc.nai_healthy=kog.mwu(nai_healthy.fc,gene2kog) 
fc.nai_healthy 

diseased_healthy=load('diseased_healthy_fc.RData')
diseased_healthy # names of datasets in the package
fc.diseased_healthy=kog.mwu(diseased_healthy.fc,gene2kog) 
fc.diseased_healthy 

diseased_nai=load('diseased_nai_fc.RData')
diseased_nai # names of datasets in the package
fc.diseased_nai=kog.mwu(diseased_nai.fc,gene2kog) 
fc.diseased_nai

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased_healthy"=fc.diseased_healthy,"nai_healthy"=fc.nai_healthy,"diseased_nai"=fc.diseased_nai))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_transmission_mcav_host_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_transmission_mcav_host_corr_fc.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="diseased_healthy",y="nai_healthy",ktable)
corrPlot(x="diseased_healthy",y="diseased_nai",ktable)
corrPlot(x="nai_healthy",y="diseased_nai",ktable)
dev.off()