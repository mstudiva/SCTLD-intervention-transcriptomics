#### packages ####

# install.packages("KOGMWU")
library(KOGMWU)


#### pairwise treatments (lpv) ####

# loading KOG annotations
gene2kog=read.table("Mcavernosa2015_iso2kogClass.tab",sep="\t", fill=T)
head(gene2kog)

healthy1_healthy0=load('healthy1_healthy0_lpv.RData')
healthy1_healthy0 # names of datasets in the package
lpv.healthy1_healthy0=kog.mwu(healthy1_healthy0.p,gene2kog) 
lpv.healthy1_healthy0 

diseased0_healthy0=load('diseased0_healthy0_lpv.RData')
diseased0_healthy0 # names of datasets in the package
lpv.diseased0_healthy0=kog.mwu(diseased0_healthy0.p,gene2kog) 
lpv.diseased0_healthy0 

treated1_healthy1=load('treated1_healthy1_lpv.RData')
treated1_healthy1 # names of datasets in the package
lpv.treated1_healthy1=kog.mwu(treated1_healthy1.p,gene2kog) 
lpv.treated1_healthy1

treated1_healthy0=load('treated1_healthy0_lpv.RData')
treated1_healthy0 # names of datasets in the package
lpv.treated1_healthy0=kog.mwu(treated1_healthy0.p,gene2kog) 
lpv.treated1_healthy0 

treated1_diseased0=load('treated1_diseased0_lpv.RData')
treated1_diseased0 # names of datasets in the package
lpv.treated1_diseased0=kog.mwu(treated1_diseased0.p,gene2kog) 
lpv.treated1_diseased0

diseased0_healthy1=load('diseased0_healthy1_lpv.RData')
diseased0_healthy1 # names of datasets in the package
lpv.diseased0_healthy1=kog.mwu(diseased0_healthy1.p,gene2kog) 
lpv.diseased0_healthy1

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased0_healthy0"=lpv.diseased0_healthy0,"healthy1_healthy0"=lpv.healthy1_healthy0,"treated1_healthy1"=lpv.treated1_healthy1,"treated1_healthy0"=lpv.treated1_healthy0,"treated1_diseased0"=lpv.treated1_diseased0,"diseased0_healthy1"=lpv.diseased0_healthy1))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_intervention_mcav_host_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_intervention_mcav_host_corr_lpv.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="diseased0_healthy0",y="healthy1_healthy0",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_healthy1",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_healthy0",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="diseased0_healthy0",y="diseased0_healthy1",ktable)
corrPlot(x="healthy1_healthy0",y="treated1_healthy1",ktable)
corrPlot(x="healthy1_healthy0",y="treated1_healthy0",ktable)
corrPlot(x="healthy1_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="healthy1_healthy0",y="diseased0_healthy1",ktable)
corrPlot(x="treated1_healthy1",y="treated1_healthy0",ktable)
corrPlot(x="treated1_healthy1",y="treated1_diseased0",ktable)
corrPlot(x="treated1_healthy1",y="diseased0_healthy1",ktable)
corrPlot(x="treated1_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="treated1_healthy0",y="diseased0_healthy1",ktable)
dev.off()


#### pairwise treatments (fc) ####

healthy1_healthy0=load('healthy1_healthy0_fc.RData')
healthy1_healthy0 # names of datasets in the package
fc.healthy1_healthy0=kog.mwu(healthy1_healthy0.fc,gene2kog) 
fc.healthy1_healthy0 

diseased0_healthy0=load('diseased0_healthy0_fc.RData')
diseased0_healthy0 # names of datasets in the package
fc.diseased0_healthy0=kog.mwu(diseased0_healthy0.fc,gene2kog) 
fc.diseased0_healthy0 

treated1_healthy1=load('treated1_healthy1_fc.RData')
treated1_healthy1 # names of datasets in the package
fc.treated1_healthy1=kog.mwu(treated1_healthy1.fc,gene2kog) 
fc.treated1_healthy1

treated1_healthy0=load('treated1_healthy0_fc.RData')
treated1_healthy0 # names of datasets in the package
fc.treated1_healthy0=kog.mwu(treated1_healthy0.fc,gene2kog) 
fc.treated1_healthy0 

treated1_diseased0=load('treated1_diseased0_fc.RData')
treated1_diseased0 # names of datasets in the package
fc.treated1_diseased0=kog.mwu(treated1_diseased0.fc,gene2kog) 
fc.treated1_diseased0

diseased0_healthy1=load('diseased0_healthy1_fc.RData')
diseased0_healthy1 # names of datasets in the package
fc.diseased0_healthy1=kog.mwu(diseased0_healthy1.fc,gene2kog) 
fc.diseased0_healthy1

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased0_healthy0"=fc.diseased0_healthy0,"healthy1_healthy0"=fc.healthy1_healthy0,"treated1_healthy1"=fc.treated1_healthy1,"treated1_healthy0"=fc.treated1_healthy0,"treated1_diseased0"=fc.treated1_diseased0,"diseased0_healthy1"=fc.diseased0_healthy1))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_intervention_mcav_host_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_intervention_mcav_host_corr_fc.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="diseased0_healthy0",y="healthy1_healthy0",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_healthy1",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_healthy0",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="diseased0_healthy0",y="diseased0_healthy1",ktable)
corrPlot(x="healthy1_healthy0",y="treated1_healthy1",ktable)
corrPlot(x="healthy1_healthy0",y="treated1_healthy0",ktable)
corrPlot(x="healthy1_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="healthy1_healthy0",y="diseased0_healthy1",ktable)
corrPlot(x="treated1_healthy1",y="treated1_healthy0",ktable)
corrPlot(x="treated1_healthy1",y="treated1_diseased0",ktable)
corrPlot(x="treated1_healthy1",y="diseased0_healthy1",ktable)
corrPlot(x="treated1_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="treated1_healthy0",y="diseased0_healthy1",ktable)
dev.off()


#### filtered treatments (lpv) ####

diseased0_healthy0=load('diseased0_healthy0_lpv.RData')
diseased0_healthy0 # names of datasets in the package
lpv.diseased0_healthy0=kog.mwu(diseased0_healthy0.p,gene2kog) 
lpv.diseased0_healthy0 

treated1_diseased0=load('treated1_diseased0_lpv.RData')
treated1_diseased0 # names of datasets in the package
lpv.treated1_diseased0=kog.mwu(treated1_diseased0.p,gene2kog) 
lpv.treated1_diseased0 

treated1_healthy1=load('treated1_healthy1_lpv.RData')
treated1_healthy1 # names of datasets in the package
lpv.treated1_healthy1=kog.mwu(treated1_healthy1.p,gene2kog) 
lpv.treated1_healthy1

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased0_healthy0"=lpv.diseased0_healthy0,"treated1_diseased0"=lpv.treated1_diseased0,"treated1_healthy1"=lpv.treated1_healthy1))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_intervention_mcav_host_filtered_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_intervention_mcav_host_filtered_corr_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="diseased0_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_healthy1",ktable)
corrPlot(x="treated1_diseased0",y="treated1_healthy1",ktable)
dev.off()


#### filtered treatments (fc) ####

diseased0_healthy0=load('diseased0_healthy0_fc.RData')
diseased0_healthy0 # names of datasets in the package
fc.diseased0_healthy0=kog.mwu(diseased0_healthy0.fc,gene2kog) 
fc.diseased0_healthy0 

treated1_diseased0=load('treated1_diseased0_fc.RData')
treated1_diseased0 # names of datasets in the package
fc.treated1_diseased0=kog.mwu(treated1_diseased0.fc,gene2kog) 
fc.treated1_diseased0 

treated1_healthy1=load('treated1_healthy1_fc.RData')
treated1_healthy1 # names of datasets in the package
fc.treated1_healthy1=kog.mwu(treated1_healthy1.fc,gene2kog) 
fc.treated1_healthy1

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("diseased0_healthy0"=fc.diseased0_healthy0,"treated1_diseased0"=fc.treated1_diseased0,"treated1_healthy1"=fc.treated1_healthy1))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_intervention_mcav_host_filtered_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_intervention_mcav_host_filtered_corr_fc.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="diseased0_healthy0",y="treated1_diseased0",ktable)
corrPlot(x="diseased0_healthy0",y="treated1_healthy1",ktable)
corrPlot(x="treated1_diseased0",y="treated1_healthy1",ktable)
dev.off()