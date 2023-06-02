#### PACKAGES ####

# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install("DESeq2",dependencies=T)
# BiocManager::install("arrayQualityMetrics",dependencies=T)  # requires Xquartz, xquartz.org
# BiocManager::install("BiocParallel")

# install.packages("pheatmap")
# install.packages("VennDiagram")
# install.packages("gplots")
# install.packages("vegan")
# install.packages("plotrix")
# install.packages("ape")
# install.packages("ggplot2")
# install.packages("rgl")
# install.packages("adegenet")


#### DATA IMPORT ####
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)

#read in counts
counts = read.table("allcounts_transmission_mcav_host.txt")

# how many genes we have total?
nrow(counts) 
ncol(counts)

# how does the data look? 
head(counts)
# removing the parent sample
counts <- subset(counts, select = -c(41))

keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData)
ncol(countData)
write.csv(countData, file="countData.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna)
ncol(counts4wgcna)
write.csv(counts4wgcna, file="counts4wgcna.csv")

# importing a design .csv file
design = read.csv("design_transmission_mcav.csv", head=TRUE)
design <- design[design$genotype != "Parent", ]
design
str(design)


#### MODEL DESIGN and OUTLIERS ####

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype+fate)

# reorders fate factor according to "control" vs "treatment" levels
dds$fate <- factor(dds$fate, levels = c("healthy", "nai", "diseased"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("fate"),force=T)
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers:
outs=c(39)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype+fate)
dds$fate <- factor(dds$fate, levels = c("healthy", "nai", "diseased"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="initial.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,4],design[,6],sep=".")
# renames the column names
colnames(vsd)=snames

save(vsd,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ genotype+fate)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=TRUE)
# vsd.wg=assay(rlog(wg), blind=TRUE)
head(vsd.wg)
colnames(vsd.wg)=snames
save(vsd.wg,design,file="data4wgcna.RData")


#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_transmission_mcav_host.pdf", width=15, height=15)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
# library(rgl)
library(ape)

conditions=design
conditions$fate <- factor(conditions$fate, levels = c("healthy", "nai", "diseased"))

# creating a PCoA eigenvalue matrix
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors
# copy this table for % variation explained by each axis (Relative_eig column)
dds.pcoa$values

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues 
pdf(file="PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's

# plotting PCoA by fate and treatment
pdf(file="PCoA_transmission_mcav_host.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("green","orange","red")[as.numeric(as.factor(conditions$fate))],pch=c(15,17,19)[as.numeric(as.factor(conditions$treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Fate")
ordispider(scores, conditions$fate, label=F, col=c("green","orange","red"))
legend("topright", legend=c("healthy", "NAI", "diseased"), fill = c("green","orange","red"), bty="n")
legend("topleft", legend=c("control","sctld"), pch=c(15,19), bty="n")
plot(scores[,1], scores[,2],col=c("green","black","red")[as.numeric(as.factor(conditions$treatment))],pch=c(15,17,19)[as.numeric((as.factor(conditions$fate)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$treatment, label=F, col=c("green","black","red"))
legend("topleft", legend=c("control", "sctld"), fill = c("green","red"), bty="n")
legend("topright", legend=c("healthy","NAI","diseased"), pch=c(15,17,19), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
ad=adonis(t(vsd)~genotype+fate,data=conditions,method="manhattan",permutations=1e6)
ad

# creating pie chart to represent ANOVA results
cols=c("blue","orange","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:3],labels=row.names(ad$aov.tab)[1:4],col=cols,main="genotype vs fate")
dev.off()


#### DESEQ ####

# with multi-factor, multi-level design - using LRT
load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# model for the effect of fate: (>2 factor levels => LRT)
dds$fate <- factor(dds$fate, levels = c("healthy","nai","diseased"))
dds_fate=DESeq(dds,test="LRT",reduced=~genotype, parallel=TRUE)

# saving all models
save(dds,dds_fate,file="realModels.RData")


#### DEGs and CONTRASTS ####

load("realModels.RData")
library(DESeq2)

# fate factor
fate=results(dds_fate) 
summary(fate) 
degs_fate=row.names(fate)[fate$padj<0.1 & !(is.na(fate$padj))]

# genotype factor
genotype=results(dds) 
summary(genotype) 
degs_genotype=row.names(genotype)[genotype$padj<0.1 & !(is.na(genotype$padj))]

# fate contrasts
diseased_healthy=results(dds,contrast=c("fate","diseased","healthy"))
summary(diseased_healthy)
degs_diseased_healthy=row.names(diseased_healthy)[diseased_healthy$padj<0.1 & !(is.na(diseased_healthy$padj))]

nai_healthy=results(dds,contrast=c("fate","nai","healthy"))
summary(nai_healthy)
degs_nai_healthy=row.names(nai_healthy)[nai_healthy$padj<0.1 & !(is.na(nai_healthy$padj))]

diseased_nai=results(dds,contrast=c("fate","diseased","nai"))
summary(diseased_nai)
degs_diseased_nai=row.names(diseased_nai)[diseased_nai$padj<0.1 & !(is.na(diseased_nai$padj))]

save(fate, genotype, diseased_healthy, nai_healthy, diseased_nai,file="pvals.RData")

# density plots: are my DEGs high-abundant or low-abundant?
load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs_genotype]),col="blue")
lines(density(means[degs_fate]),col="orange")
legend("topright", title = "Factor", legend=c("genotype","fate"), fill = c("blue","orange"))
dev.off()


#### VENN DIAGRAMS ####

load("pvals.RData")
library(DESeq2)

candidates=list("genotype"=degs_genotype, "fate"=degs_fate)

# install.packages("VennDiagram")
library(VennDiagram)

# overall factors, full model
fullmodel_venn=venn.diagram(
	x = candidates,
	filename=NULL,
	col = "transparent",
	fill = c("blue", "orange"),
	alpha = 0.5,
	label.col = c("darkblue", "white", "darkred"),
	cex = 3,
	fontfamily = "sans",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col =c("darkblue", "darkred"),
	cat.cex = 3,
	cat.fontfamily = "sans",
	cat.dist = c(0.06, 0.06),
	cat.pos = 3
	)
pdf(file="Venn_transmission_mcav_host.pdf", height=6, width=6)
grid.draw(fullmodel_venn)
dev.off()

pairwise=list("diseased_healthy"=degs_diseased_healthy,"nai_healthy"=degs_nai_healthy, "diseased_nai"=degs_diseased_nai)

# overall factors, full model
pairwise.venn=venn.diagram(
  x = pairwise,
  filename=NULL,
  col = "transparent",
  fill = c("blue", "orange", "lightblue"),
  alpha = 0.5,
  label.col = c("darkblue", "white", "darkred", "white", "white", "white", "cornflowerblue"),
  cex = 3,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkblue", "darkred", "cornflowerblue"),
  cat.cex = 3,
  cat.fontfamily = "sans",
  cat.dist = c(0.06, 0.06, -0.06),
  cat.pos = 3
)
pdf(file="Venn_transmission_mcav_host_pairwise.pdf", height=8, width=8)
grid.draw(pairwise.venn)
dev.off()


#### GO/KOG EXPORT ####

load("realModels.RData")
load("pvals.RData")

# fold change (fc) can only be used for binary factors, such as control/treatment, or specific contrasts comparing two factor levels
# log p value (lpv) is for multi-level factors, including binary factors

# genotype factor
# signed log p-values: -log(pvalue)* direction:
source=genotype[!is.na(genotype$pvalue),]
genotype.p=data.frame("gene"=row.names(source))
genotype.p$lpv=-log(source[,"pvalue"],10)
genotype.p$lpv[source$stat<0]=genotype.p$lpv[source$stat<0]*-1
head(genotype.p)
write.csv(genotype.p,file="genotype_lpv.csv",row.names=F,quote=F)
save(genotype.p,file="genotype_lpv.RData")

# fate factor
# signed log p-values: -log(pvalue)* direction:
source=fate[!is.na(fate$pvalue),]
fate.p=data.frame("gene"=row.names(source))
fate.p$lpv=-log(source[,"pvalue"],10)
fate.p$lpv[source$stat<0]=fate.p$lpv[source$stat<0]*-1
head(fate.p)
write.csv(fate.p,file="fate_lpv.csv",row.names=F,quote=F)
save(fate.p,file="fate_lpv.RData")

# fate contrasts

# diseased vs healthy
# log2 fold changes:
source=diseased_healthy[!is.na(diseased_healthy$pvalue),]
diseased_healthy.fc=data.frame("gene"=row.names(source))
diseased_healthy.fc$lfc=source[,"log2FoldChange"]
head(diseased_healthy.fc)
write.csv(diseased_healthy.fc,file="diseased_healthy_fc.csv",row.names=F,quote=F)
save(diseased_healthy.fc,file="diseased_healthy_fc.RData")

# signed log p-values: -log(pvalue)* direction:
diseased_healthy.p=data.frame("gene"=row.names(source))
diseased_healthy.p$lpv=-log(source[,"pvalue"],10)
diseased_healthy.p$lpv[source$stat<0]=diseased_healthy.p$lpv[source$stat<0]*-1
head(diseased_healthy.p)
write.csv(diseased_healthy.p,file="diseased_healthy_lpv.csv",row.names=F,quote=F)
save(diseased_healthy.p,file="diseased_healthy_lpv.RData")

# nai vs healthy
# log2 fold changes:
source=nai_healthy[!is.na(nai_healthy$pvalue),]
nai_healthy.fc=data.frame("gene"=row.names(source))
nai_healthy.fc$lfc=source[,"log2FoldChange"]
head(nai_healthy.fc)
write.csv(nai_healthy.fc,file="nai_healthy_fc.csv",row.names=F,quote=F)
save(nai_healthy.fc,file="nai_healthy_fc.RData")

# signed log p-values: -log(pvalue)* direction:
nai_healthy.p=data.frame("gene"=row.names(source))
nai_healthy.p$lpv=-log(source[,"pvalue"],10)
nai_healthy.p$lpv[source$stat<0]=nai_healthy.p$lpv[source$stat<0]*-1
head(nai_healthy.p)
write.csv(nai_healthy.p,file="nai_healthy_lpv.csv",row.names=F,quote=F)
save(nai_healthy.p,file="nai_healthy_lpv.RData")

# diseased vs nai
# log2 fold changes:
source=diseased_nai[!is.na(diseased_nai$pvalue),]
diseased_nai.fc=data.frame("gene"=row.names(source))
diseased_nai.fc$lfc=source[,"log2FoldChange"]
head(diseased_nai.fc)
write.csv(diseased_nai.fc,file="diseased_nai_fc.csv",row.names=F,quote=F)
save(diseased_nai.fc,file="diseased_nai_fc.RData")

# signed log p-values: -log(pvalue)* direction:
diseased_nai.p=data.frame("gene"=row.names(source))
diseased_nai.p$lpv=-log(source[,"pvalue"],10)
diseased_nai.p$lpv[source$stat<0]=diseased_nai.p$lpv[source$stat<0]*-1
head(diseased_nai.p)
write.csv(diseased_nai.p,file="diseased_nai_lpv.csv",row.names=F,quote=F)
save(diseased_nai.p,file="diseased_nai_lpv.RData")


#### CHERRY PICKING ####

diseased_healthy.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) -> cherrypicking
write.csv(cherrypicking, file = "trans_mcav_cherrypicking.csv")


#### GENE BOXPLOTS ####

library(DESeq2)
library(ggpubr)
load("realModels.RData")

# exporting counts of specific genes from immune-related searches
Mcavernosa14879 <- plotCounts(dds, gene="Mcavernosa14879", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa14879, file = "Mcavernosa14879.csv")

Mcavernosa98966 <- plotCounts(dds, gene="Mcavernosa98966", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa98966, file = "Mcavernosa98966.csv")

Mcavernosa64646 <- plotCounts(dds, gene="Mcavernosa64646", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa64646, file = "Mcavernosa64646.csv")

Mcavernosa69648 <- plotCounts(dds, gene="Mcavernosa69648", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa69648, file = "Mcavernosa69648.csv")

Mcavernosa9810 <- plotCounts(dds, gene="Mcavernosa9810", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa9810, file = "Mcavernosa9810.csv")

Mcavernosa43816 <- plotCounts(dds, gene="Mcavernosa43816", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa43816, file = "Mcavernosa43816.csv")

Mcavernosa54510 <- plotCounts(dds, gene="Mcavernosa54510", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa54510, file = "Mcavernosa54510.csv")

Mcavernosa47647 <- plotCounts(dds, gene="Mcavernosa47647", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa47647, file = "Mcavernosa47647.csv")

Mcavernosa50735 <- plotCounts(dds, gene="Mcavernosa50735", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa50735, file = "Mcavernosa50735.csv")

Mcavernosa61972 <- plotCounts(dds, gene="Mcavernosa61972", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa61972, file = "Mcavernosa61972.csv")

Mcavernosa14843 <- plotCounts(dds, gene="Mcavernosa14843", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa14843, file = "Mcavernosa14843.csv")

Mcavernosa16729 <- plotCounts(dds, gene="Mcavernosa16729", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa16729, file = "Mcavernosa16729.csv")

Mcavernosa184695 <- plotCounts(dds, gene="Mcavernosa184695", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa184695, file = "Mcavernosa184695.csv")

Mcavernosa49548 <- plotCounts(dds, gene="Mcavernosa49548", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa49548, file = "Mcavernosa49548.csv")

Mcavernosa59808 <- plotCounts(dds, gene="Mcavernosa59808", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa59808, file = "Mcavernosa59808.csv")

Mcavernosa9650 <- plotCounts(dds, gene="Mcavernosa9650", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa9650, file = "Mcavernosa9650.csv")

Mcavernosa102943 <- plotCounts(dds, gene="Mcavernosa102943", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa102943, file = "Mcavernosa102943.csv")

Mcavernosa12949 <- plotCounts(dds, gene="Mcavernosa12949", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa12949, file = "Mcavernosa12949.csv")

Mcavernosa29964 <- plotCounts(dds, gene="Mcavernosa29964", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa29964, file = "Mcavernosa29964.csv")

Mcavernosa20827 <- plotCounts(dds, gene="Mcavernosa20827", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa20827, file = "Mcavernosa20827.csv")

Mcavernosa71973 <- plotCounts(dds, gene="Mcavernosa71973", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa71973, file = "Mcavernosa71973.csv")

Mcavernosa126229 <- plotCounts(dds, gene="Mcavernosa126229", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa126229, file = "Mcavernosa126229.csv")

Mcavernosa27667 <- plotCounts(dds, gene="Mcavernosa27667", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa27667, file = "Mcavernosa27667.csv")

Mcavernosa21718 <- plotCounts(dds, gene="Mcavernosa21718", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa21718, file = "Mcavernosa21718.csv")

Mcavernosa96261 <- plotCounts(dds, gene="Mcavernosa96261", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa96261, file = "Mcavernosa96261.csv")

Mcavernosa366959 <- plotCounts(dds, gene="Mcavernosa366959", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa366959, file = "Mcavernosa366959.csv")
