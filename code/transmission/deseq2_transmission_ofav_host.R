#### PACKAGES ####

# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
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
counts = read.table("allcounts_transmission_ofav_host.txt")

# how many genes we have total?
nrow(counts) 
ncol(counts)

# how does the data look? 
head(counts)

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
design = read.csv("design_transmission_ofav.csv", head=TRUE)
design
str(design)


#### MODEL DESIGN and OUTLIERS ####

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype+treatment)

# reorders treatment factor according to "control" vs "treatment" levels
dds$treatment <- factor(dds$treatment, levels = c("control", "sctld"))

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
outs=c(5,24,34)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ genotype+treatment)
dds$treatment <- factor(dds$treatment, levels = c("control", "sctld"))

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
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ genotype+treatment)
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
pdf(file="heatmap_transmission_ofav_host.pdf", width=15, height=15)
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
pdf(file="PCoA_transmission_ofav_host.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("green","orange","red")[as.numeric(as.factor(conditions$fate))],pch=c(15,19)[as.numeric(as.factor(conditions$treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Fate")
ordispider(scores, conditions$fate, label=F, col=c("green","orange","red"))
legend("topright", legend=c("healthy", "NAI", "diseased"), fill = c("green","orange","red"), bty="n")
legend("topleft", legend=c("control","sctld"), pch=c(15,19), bty="n")
plot(scores[,1], scores[,2],col=c("green","red")[as.numeric(as.factor(conditions$treatment))],pch=c(15,17,19)[as.numeric((as.factor(conditions$fate)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
ordispider(scores, conditions$treatment, label=F, col=c("green","red"))
legend("topleft", legend=c("control", "sctld"), fill = c("green","red"), bty="n")
legend("topright", legend=c("healthy","NAI","diseased"), pch=c(15,17,19), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
ad=adonis(t(vsd)~genotype+treatment,data=conditions,method="manhattan",permutations=1e6)
ad

# creating pie chart to represent ANOVA results
cols=c("blue","orange","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:3],labels=row.names(ad$aov.tab)[1:4],col=cols,main="genotype vs treatment")
dev.off()


#### DESEQ ####

# with multi-factor, multi-level design - using LRT
load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# model for the effect of treatment: (>2 factor levels => LRT)
dds$treatment <- factor(dds$treatment, levels = c("control", "sctld"))
dds_treat=DESeq(dds,test="LRT",reduced=~genotype, parallel=TRUE)

# saving all models
save(dds,dds_treat,file="realModels.RData")


#### DEGs and CONTRASTS ####

load("realModels.RData")
library(DESeq2)

# treatment factor
treatment=results(dds_treat) 
summary(treatment) 
degs_treat=row.names(treatment)[treatment$padj<0.1 & !(is.na(treatment$padj))]

# genotype factor
genotype=results(dds) 
summary(genotype) 
degs_genotype=row.names(genotype)[genotype$padj<0.1 & !(is.na(genotype$padj))]

# treatment contrasts
diseased_healthy=results(dds,contrast=c("treatment","sctld","control"))
summary(diseased_healthy)
degs_diseased_healthy=row.names(diseased_healthy)[diseased_healthy$padj<0.1 & !(is.na(diseased_healthy$padj))]

save(treatment, genotype, diseased_healthy, file="pvals.RData")

# density plots: are my DEGs high-abundant or low-abundant?
load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs_genotype]),col="blue")
lines(density(means[degs_treat]),col="orange")
legend("topright", title = "Factor", legend=c("genotype","treatment"), fill = c("blue","orange"))
dev.off()


#### VENN DIAGRAMS ####

load("pvals.RData")
library(DESeq2)

candidates=list("genotype"=degs_genotype, "treatment"=degs_treat)

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
pdf(file="Venn_transmission_ofav_host.pdf", height=6, width=6)
grid.draw(fullmodel_venn)
dev.off()

pairwise=list("genotype"=degs_genotype, "treatment"=degs_treat, "diseased_healthy"=degs_diseased_healthy)

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
pdf(file="Venn_transmission_ofav_host_pairwise.pdf", height=8, width=8)
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

# treatment factor
# signed log p-values: -log(pvalue)* direction:
source=treatment[!is.na(treatment$pvalue),]
treatment.p=data.frame("gene"=row.names(source))
treatment.p$lpv=-log(source[,"pvalue"],10)
treatment.p$lpv[source$stat<0]=treatment.p$lpv[source$stat<0]*-1
head(treatment.p)
write.csv(treatment.p,file="treatment_lpv.csv",row.names=F,quote=F)
save(treatment.p,file="treatment_lpv.RData")

# log2 fold changes:
treatment.fc=data.frame("gene"=row.names(source))
treatment.fc$lfc=source[,"log2FoldChange"]
head(treatment.fc)
write.csv(treatment.fc,file="treatment_fc.csv",row.names=F,quote=F)
save(treatment.fc,file="treatment_fc.RData")


# treatment contrasts

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

#### CHERRY PICKING ####

diseased_healthy.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../annotate/ofav/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) -> cherrypicking
write.csv(cherrypicking, file = "trans_ofav_cherrypicking.csv")
