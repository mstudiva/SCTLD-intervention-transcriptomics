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
counts = read.table("allcounts_intervention_mcav_host.txt")

# how many genes we have total?
nrow(counts) 
ncol(counts)

# how does the data look? 
head(counts)

# filtering out low-count genes
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
design = read.csv("design_intervention.csv", head=TRUE)
design$time <- as.factor(design$time)
design
str(design)


#### MODEL DESIGN and OUTLIERS ####

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ treatment.time)

# reorders fate factor according to "control" vs "treatment" levels
dds$treatment.time <- factor(dds$treatment.time, levels = c("control.0","control.1","sctld.0","sctld.1"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("treatment.time"),force=T)
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers:
outs=c(1,2,9,12,21,25,33,37)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ treatment.time)
dds$treatment.time <- factor(dds$treatment.time, levels = c("control.0","control.1","sctld.0","sctld.1"))

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
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ treatment.time)
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
pdf(file="heatmap_intervention_host.pdf", width=15, height=15)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
# library(rgl)
library(ape)

conditions=design
conditions$treatment.time <- factor(conditions$treatment.time, levels = c("control.0","control.1","sctld.0","sctld.1"))

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

# plotting PCoA by fate and time
pdf(file="PCoA_intervention_host.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("red","green","orange")[as.numeric(as.factor(conditions$fate))],pch=c(1,19)[as.numeric(as.factor(conditions$time))], xlab="Coordinate 1", ylab="Coordinate 2", main="Fate")
ordispider(scores, conditions$fate, label=F, col=c("red","green","orange"))
legend("topright", legend=c("healthy", "diseased", "treated"), fill = c("green","red","orange"), bty="n")
legend("topleft", legend=c("t0","t1"), pch=c(1,19), bty="n")
plot(scores[,1], scores[,2],col=c("grey","black")[as.numeric(as.factor(conditions$time))],pch=c(15,17,25)[as.numeric((as.factor(conditions$fate)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Time")
ordispider(scores, conditions$time, label=F, col=c("grey","black"))
legend("topleft", legend=c("t0", "t1"), fill = c("grey","black"), bty="n")
legend("topright", legend=c("healthy","diseased","treated"), pch=c(15,17,25), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
ad=adonis2(t(vsd)~time*treatment,data=conditions,method="manhattan",permutations=1e6)
ad

# creating pie chart to represent ANOVA results
cols=c("blue","orange","lightblue","grey80")
pdf(file="ANOVA_pie.pdf", width=6, height=6)
pie(ad$aov.tab$R2[1:4],labels=row.names(ad$aov.tab)[1:4],col=cols,main="time vs treatment")
dev.off()


#### DAPC ####

library(adegenet)
library(parallel)
# detectCores()
library(dplyr)
library(tidyr)
library(stringr)
load("vsd.RData")

conditions=design
conditions$treatment.time <- factor(conditions$treatment.time, levels = c("sctld.0","sctld.1","control.0","control.1"))

# runs simulations on randomly-chosen datasets of 90% of the total dataset to test the number of PCs to retain
set.seed(999)
# by depth, excluding transplants
xvalDapc(t(vsd[,conditions$treatment.time!="sctld.1"]),conditions$treatment.time[conditions$treatment.time!="sctld.1"], n.rep=100, parallel="multicore", ncpus= 12)
# This tells us 15 PCs is the most successful in terms of correct assignment, but we need to test again with a smaller range of possible PCs and more reps
xvalDapc(t(vsd[,conditions$treatment.time!="sctld.1"]),conditions$treatment.time[conditions$treatment.time!="sctld.1"], n.rep=1000, n.pca=10:20, parallel="multicore", ncpus= 12)
# 15 PCs 

# now running the dapc without transplants
dp=dapc(t(vsd[,conditions$treatment.time!="sctld.1"]),conditions$treatment.time[conditions$treatment.time!="sctld.1"],n.pca=15, n.da=1)

# can we predict depth treatment for the transplants?
pred=predict.dapc(dp,newdata=(t(vsd[,conditions$treatment.time=="sctld.1"])))
pred$posterior
# look at the posterior section for assignments by probability
write.csv(pred$posterior, file = "dapc_intervention_host.csv")

# creating a new dapc object to add in transplants for plotting
amox=dp
amox$ind.coord=pred$ind.scores
amox$posterior=pred$posterior
amox$assign=pred$assign
amox$grp<-as.factor(c("sctld.1","sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1", "sctld.1"))

# now exporting side by side figures of controls vs transplants
# use Adobe Illustrator to overlay transplants curve onto controls
pdf(file="DAPC_intervention_host.pdf", width=12, height=6)
par(mfrow=c(1,2))
scatter(dp, bg="white",scree.da=FALSE,legend=TRUE,solid=0.6, col= c("red","orange","greenyellow","green"))
scatter(amox, bg="white",scree.da=FALSE,legend=FALSE,solid=0.6, col= "orange")
dev.off()

# rearranging for significance testing below
dpc=data.frame(rbind(dp$ind.coord,pred$ind.scores))
dpc$sampleid <- rownames(dpc)

# adding factor column
conditions %>%
  unite('sampleid', id, treatment.time, genotype, sep=".", remove = FALSE) %>%
  mutate(sampleid = str_replace(sampleid,'-','.'))-> conditions

conditions %>%
  select(sampleid, treatment.time) -> lookup

dpc <- left_join(dpc, lookup)

# testing significance of DFA differences with MCMCglmm
# install.packages("MCMCglmm")
library(MCMCglmm)

# sets prior distribution and creates a glm
prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V=1, nu=0.002,alpha.mu=0, alpha.V=1000)))
glm <-MCMCglmm(LD1~treatment.time,random=~sampleid, family="gaussian", data=dpc,prior=prior,nitt=75000,thin=25,burnin=5000)
summary(glm)
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)                -1.517   -1.960   -1.040     2404 <4e-04 ***
#   treatment.timesctld.1       1.319    0.685    1.978     2404 <4e-04 ***
#   treatment.timecontrol.0     2.569    1.899    3.219     2630 <4e-04 ***
#   treatment.timecontrol.1     1.972    1.272    2.573     2575 <4e-04 ***
# check to make sure you don't have autocorrelation with the reps (shown as "walks" in model traces)
plot(glm)

# calculating difference in magnitudes of t1 treated state and t1 healthy state using sampled sets of parameters
delta=abs(glm$Sol[,"treatment.timesctld.1"])-abs(glm$Sol[,"treatment.timecontrol.1"])
# 95% credible interval
HPDinterval(delta)
# lower      upper
# var1 -1.2965 0.02663611

# MCMC p-value
if (is.na(table(delta<0)[2])) {
  cat("p <",signif(1/length(delta),1))
} else { cat("p =",signif(table(delta<0)[2]/length(delta),2)) }
# p = 0.98
# non significant p value indicates treated corals are indistinguishable from healthy controls at t1


#### DESEQ ####

# with multi-factor, multi-level design - using LRT
load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# saving all models
save(dds,file="realModels.RData")


#### DEGs and CONTRASTS ####

load("realModels.RData")
library(DESeq2)

# treatment
treatment_time=results(dds) 
summary(treatment_time) 
degs_treatment_time=row.names(treatment_time)[treatment_time$padj<0.1 & !(is.na(treatment_time$padj))]

# treatment and time contrasts
healthy1_healthy0=results(dds,contrast=c("treatment.time","control.1","control.0"))
summary(healthy1_healthy0)
degs_healthy1_healthy0=row.names(healthy1_healthy0)[healthy1_healthy0$padj<0.1 & !(is.na(healthy1_healthy0$padj))]

diseased0_healthy0=results(dds,contrast=c("treatment.time","sctld.0","control.0"))
summary(diseased0_healthy0)
degs_diseased0_healthy0=row.names(diseased0_healthy0)[diseased0_healthy0$padj<0.1 & !(is.na(diseased0_healthy0$padj))]

treated1_healthy1=results(dds,contrast=c("treatment.time","sctld.1","control.1"))
summary(treated1_healthy1)
degs_treated1_healthy1=row.names(treated1_healthy1)[treated1_healthy1$padj<0.1 & !(is.na(treated1_healthy1$padj))]

treated1_healthy0=results(dds,contrast=c("treatment.time","sctld.1","control.0"))
summary(treated1_healthy0)
degs_treated1_healthy0=row.names(treated1_healthy0)[treated1_healthy0$padj<0.1 & !(is.na(treated1_healthy0$padj))]

treated1_diseased0=results(dds,contrast=c("treatment.time","sctld.1","sctld.0"))
summary(treated1_diseased0)
degs_treated1_diseased0=row.names(treated1_diseased0)[treated1_diseased0$padj<0.1 & !(is.na(treated1_diseased0$padj))]

diseased0_healthy1=results(dds,contrast=c("treatment.time","sctld.0","control.1"))
summary(diseased0_healthy1)
degs_diseased0_healthy1=row.names(diseased0_healthy1)[diseased0_healthy1$padj<0.1 & !(is.na(diseased0_healthy1$padj))]

save(treatment_time,healthy1_healthy0,diseased0_healthy0,treated1_healthy1,treated1_healthy0,treated1_diseased0,diseased0_healthy1,degs_treatment_time,degs_healthy1_healthy0,degs_diseased0_healthy0,degs_treated1_healthy1,degs_treated1_healthy0,degs_treated1_diseased0,degs_diseased0_healthy1,file="pvals.RData")

# density plots: are my DEGs high-abundant or low-abundant?
load("vsd.RData")
load("pvals.RData")

means=apply(vsd,1,mean)

pdf(file="DEG_density_treatment.time.pdf", height=5, width=5)
plot(density(means))
lines(density(means[degs_healthy1_healthy0]),col="blue")
lines(density(means[degs_diseased0_healthy0]),col="orange")
lines(density(means[degs_treated1_healthy1]),col="lightblue")
lines(density(means[degs_treated1_healthy0]),col="yellow")
lines(density(means[degs_treated1_diseased0]),col="red")
lines(density(means[degs_diseased0_healthy1]),col="purple")
legend("topright", title = "Factor", legend=c("healthy1_healthy0","diseased0_healthy0","treated1_healthy1","treated1_healthy0","treated1_diseased0","diseased0_healthy1"), fill = c("blue","orange","lightblue","yellow","red","purple"))
dev.off()


#### VENN DIAGRAMS ####

load("pvals.RData")
library(DESeq2)

pairwise=list("H1/H0"=degs_healthy1_healthy0, "T1/H1"=degs_treated1_healthy1,"D0/H0"=degs_diseased0_healthy0,"T1/D0"=degs_treated1_diseased0)

# install.packages("VennDiagram")
library(VennDiagram)

# treatment/time contrasts
venn=venn.diagram(
  x = pairwise,
  filename=NULL,
  col = "transparent",
  fill = c("#01665e", "#5ab4ac", "#8c510a","#f6e8c3"),
  alpha = 0.5,
  label.col = c("#8c510a","white","#dfc27d","white","white","black","white", "white","#01665e","white","white","white","white","#35978f","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("#01665e","#35978f","#8c510a","#dfc27d"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_intervention_host.pdf", height=10, width=12)
grid.draw(venn)
dev.off()


#### GO/KOG EXPORT ####

load("realModels.RData")
load("pvals.RData")

# fold change (fc) can only be used for binary factors, such as control/treatment, or specific contrasts comparing two factor levels
# log p value (lpv) is for multi-level factors, including binary factors

# healthy1 vs healthy0
# log2 fold changes:
source=healthy1_healthy0[!is.na(healthy1_healthy0$pvalue),]
healthy1_healthy0.fc=data.frame("gene"=row.names(source))
healthy1_healthy0.fc$lfc=source[,"log2FoldChange"]
head(healthy1_healthy0.fc)
write.csv(healthy1_healthy0.fc,file="healthy1_healthy0_fc.csv",row.names=F,quote=F)
save(healthy1_healthy0.fc,file="healthy1_healthy0_fc.RData")

# signed log p-values: -log(pvalue)* direction:
healthy1_healthy0.p=data.frame("gene"=row.names(source))
healthy1_healthy0.p$lpv=-log(source[,"pvalue"],10)
healthy1_healthy0.p$lpv[source$stat<0]=healthy1_healthy0.p$lpv[source$stat<0]*-1
head(healthy1_healthy0.p)
write.csv(healthy1_healthy0.p,file="healthy1_healthy0_lpv.csv",row.names=F,quote=F)
save(healthy1_healthy0.p,file="healthy1_healthy0_lpv.RData")


# diseased0 vs healthy0
# log2 fold changes:
source=diseased0_healthy0[!is.na(diseased0_healthy0$pvalue),]
diseased0_healthy0.fc=data.frame("gene"=row.names(source))
diseased0_healthy0.fc$lfc=source[,"log2FoldChange"]
head(diseased0_healthy0.fc)
write.csv(diseased0_healthy0.fc,file="diseased0_healthy0_fc.csv",row.names=F,quote=F)
save(diseased0_healthy0.fc,file="diseased0_healthy0_fc.RData")

# signed log p-values: -log(pvalue)* direction:
diseased0_healthy0.p=data.frame("gene"=row.names(source))
diseased0_healthy0.p$lpv=-log(source[,"pvalue"],10)
diseased0_healthy0.p$lpv[source$stat<0]=diseased0_healthy0.p$lpv[source$stat<0]*-1
head(diseased0_healthy0.p)
write.csv(diseased0_healthy0.p,file="diseased0_healthy0_lpv.csv",row.names=F,quote=F)
save(diseased0_healthy0.p,file="diseased0_healthy0_lpv.RData")


# treated1 vs healthy1
# log2 fold changes:
source=treated1_healthy1[!is.na(treated1_healthy1$pvalue),]
treated1_healthy1.fc=data.frame("gene"=row.names(source))
treated1_healthy1.fc$lfc=source[,"log2FoldChange"]
head(treated1_healthy1.fc)
write.csv(treated1_healthy1.fc,file="treated1_healthy1_fc.csv",row.names=F,quote=F)
save(treated1_healthy1.fc,file="treated1_healthy1_fc.RData")

# signed log p-values: -log(pvalue)* direction:
treated1_healthy1.p=data.frame("gene"=row.names(source))
treated1_healthy1.p$lpv=-log(source[,"pvalue"],10)
treated1_healthy1.p$lpv[source$stat<0]=treated1_healthy1.p$lpv[source$stat<0]*-1
head(treated1_healthy1.p)
write.csv(treated1_healthy1.p,file="treated1_healthy1_lpv.csv",row.names=F,quote=F)
save(treated1_healthy1.p,file="treated1_healthy1_lpv.RData")


# treated1 vs healthy0
# log2 fold changes:
source=treated1_healthy0[!is.na(treated1_healthy0$pvalue),]
treated1_healthy0.fc=data.frame("gene"=row.names(source))
treated1_healthy0.fc$lfc=source[,"log2FoldChange"]
head(treated1_healthy0.fc)
write.csv(treated1_healthy0.fc,file="treated1_healthy0_fc.csv",row.names=F,quote=F)
save(treated1_healthy0.fc,file="treated1_healthy0_fc.RData")

# signed log p-values: -log(pvalue)* direction:
treated1_healthy0.p=data.frame("gene"=row.names(source))
treated1_healthy0.p$lpv=-log(source[,"pvalue"],10)
treated1_healthy0.p$lpv[source$stat<0]=treated1_healthy0.p$lpv[source$stat<0]*-1
head(treated1_healthy0.p)
write.csv(treated1_healthy0.p,file="treated1_healthy0_lpv.csv",row.names=F,quote=F)
save(treated1_healthy0.p,file="treated1_healthy0_lpv.RData")


# treated1 vs diseased0
# log2 fold changes:
source=treated1_diseased0[!is.na(treated1_diseased0$pvalue),]
treated1_diseased0.fc=data.frame("gene"=row.names(source))
treated1_diseased0.fc$lfc=source[,"log2FoldChange"]
head(treated1_diseased0.fc)
write.csv(treated1_diseased0.fc,file="treated1_diseased0_fc.csv",row.names=F,quote=F)
save(treated1_diseased0.fc,file="treated1_diseased0_fc.RData")

# signed log p-values: -log(pvalue)* direction:
treated1_diseased0.p=data.frame("gene"=row.names(source))
treated1_diseased0.p$lpv=-log(source[,"pvalue"],10)
treated1_diseased0.p$lpv[source$stat<0]=treated1_diseased0.p$lpv[source$stat<0]*-1
head(treated1_diseased0.p)
write.csv(treated1_diseased0.p,file="treated1_diseased0_lpv.csv",row.names=F,quote=F)
save(treated1_diseased0.p,file="treated1_diseased0_lpv.RData")


# diseased0 vs healthy1
# log2 fold changes:
source=diseased0_healthy1[!is.na(diseased0_healthy1$pvalue),]
diseased0_healthy1.fc=data.frame("gene"=row.names(source))
diseased0_healthy1.fc$lfc=source[,"log2FoldChange"]
head(diseased0_healthy1.fc)
write.csv(diseased0_healthy1.fc,file="diseased0_healthy1_fc.csv",row.names=F,quote=F)
save(diseased0_healthy1.fc,file="diseased0_healthy1_fc.RData")

# signed log p-values: -log(pvalue)* direction:
diseased0_healthy1.p=data.frame("gene"=row.names(source))
diseased0_healthy1.p$lpv=-log(source[,"pvalue"],10)
diseased0_healthy1.p$lpv[source$stat<0]=diseased0_healthy1.p$lpv[source$stat<0]*-1
head(diseased0_healthy1.p)
write.csv(diseased0_healthy1.p,file="diseased0_healthy1_lpv.csv",row.names=F,quote=F)
save(diseased0_healthy1.p,file="diseased0_healthy1_lpv.RData")


#### CHERRY PICKING ####

diseased0_healthy0.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) -> cherrypicking
write.csv(cherrypicking, file = "inter_d0h0_cherrypicking.csv")

treated1_diseased0.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'NF-kappaB|peroxidas|TGF-beta|protein tyrosine kinase|fibrinogen|WD repeat-containing protein|apoptosis|extracellular matrix')) -> cherrypicking
write.csv(cherrypicking, file = "inter_t1d0_cherrypicking.csv")


#### GENE BOXPLOTS ####

library(DESeq2)
library(ggpubr)
load("realModels.RData")

# exporting counts of specific genes from immune-related searches
Mcavernosa9810 <- plotCounts(dds, gene="Mcavernosa9810", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa9810, file = "Mcavernosa9810.csv")

Mcavernosa43816 <- plotCounts(dds, gene="Mcavernosa43816", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa43816, file = "Mcavernosa43816.csv")

Mcavernosa14879 <- plotCounts(dds, gene="Mcavernosa14879", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa14879, file = "Mcavernosa14879.csv")

Mcavernosa12872 <- plotCounts(dds, gene="Mcavernosa12872", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa12872, file = "Mcavernosa12872.csv")

Mcavernosa47647 <- plotCounts(dds, gene="Mcavernosa47647", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa47647, file = "Mcavernosa47647.csv")

Mcavernosa50735 <- plotCounts(dds, gene="Mcavernosa50735", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa50735, file = "Mcavernosa50735.csv")

Mcavernosa10679 <- plotCounts(dds, gene="Mcavernosa10679", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa10679, file = "Mcavernosa10679.csv")

Mcavernosa61972 <- plotCounts(dds, gene="Mcavernosa61972", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa61972, file = "Mcavernosa61972.csv")

Mcavernosa49779 <- plotCounts(dds, gene="Mcavernosa49779", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa49779, file = "Mcavernosa49779.csv")

Mcavernosa5069 <- plotCounts(dds, gene="Mcavernosa5069", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa5069, file = "Mcavernosa5069.csv")

Mcavernosa184695 <- plotCounts(dds, gene="Mcavernosa184695", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa184695, file = "Mcavernosa184695.csv")

Mcavernosa102943 <- plotCounts(dds, gene="Mcavernosa102943", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa102943, file = "Mcavernosa102943.csv")

Mcavernosa29964 <- plotCounts(dds, gene="Mcavernosa29964", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa29964, file = "Mcavernosa29964.csv")

Mcavernosa20827 <- plotCounts(dds, gene="Mcavernosa20827", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa20827, file = "Mcavernosa20827.csv")

Mcavernosa71973 <- plotCounts(dds, gene="Mcavernosa71973", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa71973, file = "Mcavernosa71973.csv")

Mcavernosa126229 <- plotCounts(dds, gene="Mcavernosa126229", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa126229, file = "Mcavernosa126229.csv")

Mcavernosa96261 <- plotCounts(dds, gene="Mcavernosa96261", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa96261, file = "Mcavernosa96261.csv")

Mcavernosa162020 <- plotCounts(dds, gene="Mcavernosa162020", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa162020, file = "Mcavernosa162020.csv")

Mcavernosa317111 <- plotCounts(dds, gene="Mcavernosa317111", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa317111, file = "Mcavernosa317111.csv")

Mcavernosa43057 <- plotCounts(dds, gene="Mcavernosa43057", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa43057, file = "Mcavernosa43057.csv")

Mcavernosa10151 <- plotCounts(dds, gene="Mcavernosa10151", intgroup="fate", returnData=TRUE)
write.csv(Mcavernosa10151, file = "Mcavernosa10151.csv")
