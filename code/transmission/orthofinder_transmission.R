#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### ORTHOFINDER ####

# Install orthofinder on your local machine using the tutorials (https://davidemms.github.io/menu/tutorials.html)

# Copy your translated protein fasta files (_out_PRO.fas) that you want to compare into a directory called 'orthofinder'
# If you have not already filtered by the longest contig per isogroup (by using fasta2SBH.pl during transcriptome annotation), follow step 7 of tutorial 2 above

# Run the following command in Terminal: 'orthofinder -f orthofinder/'
# Check the number of genes assigned to orthogroups (e.g., 'OrthoFinder assigned 29028 genes (87.6% of total) to 8731 orthogroups')
# Ideally, it should be >80%


#### ORTHOLOGS ####

orthologs <- read.table(file = "OrthoFinder/Results_May10/Orthologues/Orthologues_Ofaveolata_out_PRO/Ofaveolata_out_PRO__v__Mcavernosa2015_out_PRO.tsv", sep = "\t", header = TRUE, quote="", fill=FALSE)


#### DESEQ IMPORT ####

ofav_diseased_healthy_lpv <- read.csv(file = "../DESeq2/ofav/diseased_healthy_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

mcav_diseased_healthy_lpv <- read.csv(file = "../DESeq2/mcav2015/diseased_healthy_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_mcav" = lpv)

mcav_diseased_nai_lpv <- read.csv(file = "../DESeq2/mcav2015/diseased_nai_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_mcav" = lpv)

mcav_nai_healthy_lpv <- read.csv(file = "../DESeq2/mcav2015/nai_healthy_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_mcav" = lpv)


#### DEG MATCHING ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds ofav and mcav gene annotations, and 5) then pulls on corresponding KOG classes

# diseased vs healthy for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  inner_join(mcav_diseased_healthy_lpv, by = c("Protein_mcav" = "gene")) %>%
  inner_join(ofav_diseased_healthy_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_mcav) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../annotate/ofav/Ofaveolata_iso2geneName.tab",
             sep = "\t",
             quote="", fill=FALSE) %>%
               mutate(gene = V1,
                      annot_ofav = V2) %>%
               dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/ofav/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> diseased_healthy
diseased_healthy$Orthogroup <- make.unique(diseased_healthy$Orthogroup, sep = "_") 

# diseased vs nai for mcav vs diseased vs healthy for ofav
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  inner_join(mcav_diseased_nai_lpv, by = c("Protein_mcav" = "gene")) %>%
  inner_join(ofav_diseased_healthy_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_mcav) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../annotate/ofav/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/ofav/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> diseased_nai
diseased_nai$Orthogroup <- make.unique(diseased_nai$Orthogroup, sep = "_") 

# nai vs healthy for mcav vs diseased vs healthy for ofav
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  inner_join(mcav_nai_healthy_lpv, by = c("Protein_mcav" = "gene")) %>%
  inner_join(ofav_diseased_healthy_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_mcav) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../annotate/ofav/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/ofav/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> nai_healthy
nai_healthy$Orthogroup <- make.unique(nai_healthy$Orthogroup, sep = "_") 


#### KOG MATCHING ####

# filtering and summarizing DEGs by KOG class for high-level comparisons between species
diseased_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_mcav >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "mcav_up" = n) -> KOG_mcav_up

diseased_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_mcav <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "mcav_down" = n) -> KOG_mcav_down

diseased_healthy %>%
  mutate(KOG_ofav = replace(KOG_ofav, KOG_ofav == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG_ofav) %>%
  rename("KOG" = KOG_ofav, "ofav_up" = n) -> KOG_ofav_up

diseased_healthy %>%
  mutate(KOG_ofav = replace(KOG_ofav, KOG_ofav == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG_ofav) %>%
  rename("KOG" = KOG_ofav, "ofav_down" = n) -> KOG_ofav_down

# joining all KOG class sums in a single dataframe
KOG_mcav_up %>%
  inner_join(KOG_ofav_up, by = "KOG") %>%
  inner_join(KOG_mcav_down, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") -> KOG_match

# melting dataframe for plotting
KOG_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_melt

# creating a custom color palette
colorCount = length(unique(KOG_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_sum <- ggplot(KOG_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Number of DEGs") +
  theme_classic()
KOG_sum
ggsave("orthofinder KOG abundance.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)


#### VENN DIAGRAMS ####

# first creating a set of up/downregulated DEGs by species
diseased_healthy %>%
  filter(lpv_mcav >= 1) %>%
  pull(Orthogroup) -> mcav_up

diseased_healthy %>%
  filter(lpv_mcav <= -1) %>%
  pull(Orthogroup) -> mcav_down

diseased_healthy %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

diseased_healthy %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

# treatment/time contrasts
venn=venn.diagram(
  x = list("Mcav up"=mcav_up, "Mcav down"=mcav_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_orthofinder.pdf", height=10, width=12)
grid.draw(venn)
dev.off()


#### COMBINE SPECIES ####

# first creating a column of combined gene names from both species, then removing unannotated genes
diseased_healthy %>%
  unite("gene_name", annot_mcav:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  diseased_healthy


#### IF FOLLOWING THE 'VSD by SPECIES' SECTION, SKIP FROM HERE DOWN ####

# # loading gene counts for each species, then replacing species-specific gene IDs with orthogroup IDs, then removing NAI samples
# load("../DESeq2/ofav/initial.RData")
# counts_ofav <- subset(countData, rownames(countData) %in% diseased_healthy$Protein_ofav)
# design_ofav <- design 
# 
# load("../DESeq2/mcav2015/initial.RData")
# counts_mcav <- subset(countData, rownames(countData) %in% diseased_healthy$Protein_mcav)
# design_mcav <- design 
# 
# # have to detach DESeq2 because it contains some functions with the same name as dplyr
# detach("package:DESeq2", unload=TRUE)
# library(tidyverse)
# 
# # combining the design data frames
# design_comb <- rbind(design_mcav, design_ofav)
# design_comb$id <- as.factor(gsub("-",".", design_comb$id))
# design_comb$full_id <- paste(design_comb$id,design_comb$species,design_comb$fate,sep=".")
# 
# # replacing gene names with orthogroups
# counts_ofav %>%
#   rownames_to_column(var = "Protein_ofav") %>%
#   mutate(Protein_ofav = if_else(Protein_ofav %in% diseased_healthy$Protein_ofav, diseased_healthy$Orthogroup, diseased_healthy$Orthogroup)) %>%
#   column_to_rownames(var = "Protein_ofav") -> counts_ofav
# 
# counts_mcav %>%
#   rownames_to_column(var = "Protein_mcav") %>%
#   mutate(Protein_mcav = if_else(Protein_mcav %in% diseased_healthy$Protein_mcav, diseased_healthy$Orthogroup, diseased_healthy$Orthogroup)) %>%
#   column_to_rownames(var = "Protein_mcav") -> counts_mcav
# 
# # combining the two count data frames, then replacing sample IDs with full IDs (ID + factors)
# counts_comb <- cbind(counts_mcav, counts_ofav)
# colnames(counts_comb) = design_comb$full_id
# 
# # removing NAI samples and reordering design dataframe for plotting
# design_comb <- design_comb[(design_comb$fate != "nai"),]
# design_comb$species_fate <- paste(design_comb$species,design_comb$fate,sep=".")
# design_comb$species_fate <- factor(design_comb$species_fate, levels = c("Mcav.healthy","Mcav.diseased","Ofav.healthy","Ofav.diseased"))
# design_comb <- design_comb[order(design_comb$species_fate),]
# 
# # reordering counts matrix according to design dataframe
# counts_comb %>%
#   select(-contains("nai")) -> counts_comb
# head(counts_comb)
# counts_comb<- counts_comb[,order(design_comb$full_id)]
# tail(counts_comb)
# 
# library(DESeq2)
# 
# # making a combined DESeq model
# rownames(design_comb) = design_comb$full_id
# dds = DESeqDataSetFromMatrix(countData=counts_comb, colData=design_comb, design=~ genotype+fate)
# 
# # creating variance stabilized array of counts for heatmap
# Vsd=varianceStabilizingTransformation(dds)
# vsd=assay(Vsd)
# 
# save(orthologs, diseased_healthy, diseased_nai, nai_healthy, mcav_up, mcav_down, ofav_up, ofav_down, counts_ofav, design_ofav, counts_mcav, design_mcav, design_comb, counts_comb, vsd, file = "orthofinder_DEGs.RData")


#### HEATMAPS ####

library(tidyverse)
library(pheatmap)

# loading back in the RData 
load("orthofinder_DEGs.RData")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations, then filtering out non-annotated genes
gene_names <- as.data.frame(cbind(diseased_healthy$Orthogroup, diseased_healthy$gene_name))

# creating a dummy variable of fake p values (since you want to plot all genes)
pval_dummy <- c(rep(0.1,619))

# heatmap time!
pdf(file="orthofinder_heatmap_gene.pdf", height=100, width=82)
uniHeatmap(vsd=vsd,gene.names=gene_names,
           metric=pval_dummy, # metric of gene significance
           cutoff=1, 
           sort=c(1:ncol(vsd)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# just orthogroup designations (not really sure why it needs the random "g", but it doesn't run without some additional text string)
ortho_names <- as.data.frame(cbind(diseased_healthy$Orthogroup, paste(c(rep("g",619)),diseased_healthy$Orthogroup,sep=" ")))

pdf(file="orthofinder_heatmap_orthogroup.pdf", height=100, width=13)
uniHeatmap(vsd=vsd,gene.names=ortho_names,
           metric=pval_dummy, # metric of gene significance
           cutoff=1, 
           sort=c(1:ncol(vsd)), # overrides sorting of columns according to hierarchical clustering
           cex=0.8,
           pdf=F,
)
dev.off()


#### VSD by SPECIES ####

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup IDs, then removing NAI samples
load("../DESeq2/ofav/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% diseased_healthy$Protein_ofav)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% diseased_healthy$Protein_ofav, diseased_healthy$Orthogroup, diseased_healthy$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  select(-contains("nai")) %>%
  as.matrix() -> vsd_ofav

load("../DESeq2/mcav2015/vsd.RData")
design_mcav <- design
vsd_mcav <- subset(vsd, rownames(vsd) %in% diseased_healthy$Protein_mcav)

vsd_mcav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_mcav") %>%
  mutate(Protein_mcav = if_else(Protein_mcav %in% diseased_healthy$Protein_mcav, diseased_healthy$Orthogroup, diseased_healthy$Orthogroup)) %>%
  column_to_rownames(var = "Protein_mcav") %>%
  select(-contains("nai")) %>%
  as.matrix() -> vsd_mcav

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_mcav, vsd_ofav)
design_comb <- rbind(design_mcav, design_ofav)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$species,design_comb$fate,sep=".")

# removing NAI samples and reordering design dataframe for plotting
design_comb <- design_comb[(design_comb$fate != "nai"),]
design_comb$species_fate <- paste(design_comb$species,design_comb$fate,sep=".")
design_comb$species_fate <- factor(design_comb$species_fate, levels = c("Mcav.healthy","Mcav.diseased","Ofav.healthy","Ofav.diseased"))
design_comb <- design_comb[order(design_comb$species_fate),]

# reordering counts matrix according to design dataframe
head(vsd_comb)
vsd_comb<- vsd_comb[,order(design_comb$full_id)]
head(vsd_comb)

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(diseased_healthy$Orthogroup, diseased_healthy$gene_name))

save(orthologs, diseased_healthy, diseased_nai, nai_healthy, mcav_up, mcav_down, ofav_up, ofav_down, design_ofav, design_mcav, design_comb, vsd_ofav, vsd_mcav, vsd_comb, file = "orthofinder_DEGs_species.RData")
load("orthofinder_DEGs_species.RData")

# heatmaps of original vsd relationships (separated by species)

# p < 0.1 (all genes)
pdf(file="orthofinder_heatmap_p0.1.pdf", height=100, width=82)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_healthy$lpv_mcav)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 0.05
pdf(file="orthofinder_heatmap_p0.05.pdf", height=75, width=80)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_healthy$lpv_mcav)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-1.3, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 0.01
pdf(file="orthofinder_heatmap_p0.01.pdf", height=50, width=80)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_healthy$lpv_mcav)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-2, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 0.001
pdf(file="orthofinder_heatmap_p0.001.pdf", height=45, width=30)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_healthy$lpv_mcav)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-3, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 1e-6
pdf(file="orthofinder_heatmap_p1e6.pdf", height=15, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_healthy$lpv_mcav)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 1e-7
pdf(file="orthofinder_heatmap_p1e7.pdf", height=10, width=22)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_healthy$lpv_mcav)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-7, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()


#### COMMON GOs ####

trans_species <- read.csv(file = "trans_species.csv")
trans_species$direction = as.factor(trans_species$direction)
trans_species$cat = factor(trans_species$cat, levels = c("MF","BP","CC"))

species_plot <- ggplot(trans_species, aes(x = species, y = name, color = direction)) +
  geom_point(aes(size = genes)) + 
  scale_color_manual(values = c("0"="blue","1"="red")) +
  facet_grid(rows = vars(cat), scales = "free", space="free_y") + 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
species_plot

ggsave("transmission GO species.pdf", plot= species_plot, width=6.5, height=6.5, units="in", dpi=300)


#### CHERRY PICKING ####

mcav_diseased_healthy_cherry <- read.csv(file = "../DESeq2/mcav2015/trans_mcav_cherrypicking.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_mcav_dh" = lpv, "annot_mcav_dh" = annot)

mcav_diseased0_healthy0_cherry <- read.csv(file = "../../intervention/DESeq2/mcav2015/inter_d0h0_cherrypicking.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_mcav_d0h0" = lpv, "annot_mcav_d0h0" = annot)

mcav_treated1_diseased0_cherry <- read.csv(file = "../../intervention/DESeq2/mcav2015/inter_t1d0_cherrypicking.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_mcav_t1d0" = lpv, "annot_mcav_t1d0" = annot)

ofav_diseased_healthy_cherry <- read.csv(file = "../DESeq2/ofav/trans_ofav_cherrypicking.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_ofav_dh" = lpv, "annot_ofav_dh" = annot)


orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  right_join(mcav_diseased_healthy_cherry, by = c("Protein_mcav" = "gene")) %>%
  write.csv(file="mcav_diseased_healthy_cherry.csv")

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  right_join(mcav_diseased0_healthy0_cherry, by = c("Protein_mcav" = "gene")) %>%
  write.csv(file="mcav_diseased0_healthy0_cherry.csv")

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  right_join(mcav_treated1_diseased0_cherry, by = c("Protein_mcav" = "gene")) %>%
  write.csv(file="mcav_treated1_diseased0_cherry.csv")

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_mcav" = 	
           Mcavernosa2015_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_mcav, sep = ",") %>%
  unique() %>%
  right_join(ofav_diseased_healthy_cherry, by = c("Protein_ofav" = "gene")) %>%
  write.csv(file="ofav_diseased_healthy_cherry.csv")
 

#### BOXPLOTS MATCHING ORTHOGROUPS ####

library(stringr)
library(rcompanion)
library(rstatix)
library(ggpubr)
library(scales)

# DEGs with matching orthogroups
# OG0005649 (Spondin 2b, extracellular matrix protein)
Mcavernosa14879 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa14879.csv")
Mcavernosa14879$fate <- factor(Mcavernosa14879$fate, levels = c("healthy", "nai", "diseased"))

Ofaveolata224187 <- read.csv(file = "../DESeq2/ofav/Ofaveolata224187.csv")
# renaming Ofav treatment to fate
colnames(Ofaveolata224187)[3] = "fate"
Ofaveolata224187$fate <- str_replace(Ofaveolata224187$fate, "control", "healthy")
Ofaveolata224187$fate <- str_replace(Ofaveolata224187$fate, "sctld", "diseased")

# ANOVA and Tukey's
Mcavernosa14879_stats <- aov(count~fate,data=Mcavernosa14879) %>%
  tukey_hsd()
Mcavernosa14879_stats

Ofaveolata224187_stats <- aov(count~fate,data=Ofaveolata224187) %>%
  tukey_hsd()
Ofaveolata224187_stats

Mcavernosa14879_plot <-
  ggboxplot(
    Mcavernosa14879,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "M. cavernosa",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(3,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa14879_stats,label="p.adj.signif",y.position=3.45,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa14879_plot

Ofaveolata224187_plot <-
  ggboxplot(
    Ofaveolata224187,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    title = "O. faveolata",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = element_blank(),
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank()) +
  scale_y_log10(limit = c(3,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata224187_stats,label="p.adj.signif",y.position=3.8,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata224187_plot

# OG0000472 (Animal haem peroxidase)
Mcavernosa98966 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa98966.csv")
Mcavernosa98966$fate <- factor(Mcavernosa98966$fate, levels = c("healthy", "nai", "diseased"))

Ofaveolata263192 <- read.csv(file = "../DESeq2/ofav/Ofaveolata263192.csv")
# renaming Ofav treatment to fate
colnames(Ofaveolata263192)[3] = "fate"
Ofaveolata263192$fate <- str_replace(Ofaveolata263192$fate, "control", "healthy")
Ofaveolata263192$fate <- str_replace(Ofaveolata263192$fate, "sctld", "diseased")

# ANOVA and Tukey's
Mcavernosa98966_stats <- aov(count~fate,data=Mcavernosa98966) %>%
  tukey_hsd()
Mcavernosa98966_stats

Ofaveolata263192_stats <- aov(count~fate,data=Ofaveolata263192) %>%
  tukey_hsd()
Ofaveolata263192_stats

Mcavernosa98966_plot <-
  ggboxplot(
    Mcavernosa98966,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "M. cavernosa",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") +
  scale_y_log10(limit = c(1.95,15000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa98966_stats,label="p.adj.signif",y.position=2,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa98966_plot

Ofaveolata263192_plot <-
  ggboxplot(
    Ofaveolata263192,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    title = "O. faveolata",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = element_blank(),
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.y=element_blank()) +
  scale_y_log10(limit = c(1.95,15000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata263192_stats,label="p.adj.signif",y.position=4.15,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata263192_plot


#### MULTIPLOT MATCHING ####

transmission_match<-ggarrange(Mcavernosa14879_plot,
                           Ofaveolata224187_plot,
                           Mcavernosa98966_plot,
                           Ofaveolata263192_plot,
                      heights = c(4,4),
                      widths = c(5,3),
                      ncol = 2,
                      nrow = 2)
transmission_match<-annotate_figure(transmission_match, top = text_grob("Spondin 2b, extracellular matrix protein", color = "black", face = "bold", size = 14), 
                bottom = text_grob("Animal haem peroxidase", color = "black", face = "bold", size = 14))
transmission_match

ggsave("transmission_match.pdf", transmission_match, width=8, height=8,dpi = 300)


#### BOXPLOTS MCAV ####

# non-matching DEGs by species
# Mcavernosa64646 (Negative regulation of cysteine-type endopeptidase activity involved in execution phase of apoptosis)
Mcavernosa64646 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa64646.csv")
Mcavernosa64646$fate <- factor(Mcavernosa64646$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa64646$title <- "Negative regulation of cysteine-type endopeptidase activity"

# ANOVA and Tukey's
Mcavernosa64646_stats <- aov(count~fate,data=Mcavernosa64646) %>%
  tukey_hsd()
Mcavernosa64646_stats

Mcavernosa64646_plot <-
  ggboxplot(
    Mcavernosa64646,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="lightskyblue"), strip.text = element_text(size=9), 
        axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa64646_stats,label="p.adj.signif",y.position=0.9,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa64646_plot

# Mcavernosa69648 (Baculoviral inhibition of apoptosis protein repeat)
Mcavernosa69648 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa69648.csv")
Mcavernosa69648$fate <- factor(Mcavernosa69648$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa69648$title <- "Baculoviral inhibition of apoptosis protein repeat"

# ANOVA and Tukey's
Mcavernosa69648_stats <- aov(count~fate,data=Mcavernosa69648) %>%
  tukey_hsd()
Mcavernosa69648_stats

Mcavernosa69648_plot <-
  ggboxplot(
    Mcavernosa69648,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="lightskyblue"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa69648_stats,label="p.adj.signif",y.position=1.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa69648_plot

# Mcavernosa9810 (Extracellular matrix binding)
Mcavernosa9810 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa9810.csv")
Mcavernosa9810$fate <- factor(Mcavernosa9810$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa9810$title <- "Extracellular matrix binding"

# ANOVA and Tukey's
Mcavernosa9810_stats <- aov(count~fate,data=Mcavernosa9810) %>%
  tukey_hsd()
Mcavernosa9810_stats

Mcavernosa9810_plot <-
  ggboxplot(
    Mcavernosa9810,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
     add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="salmon"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa9810_stats,label="p.adj.signif",y.position=2.9,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa9810_plot

# Mcavernosa43816 (Spondin 2a, extracellular matrix protein)
Mcavernosa43816 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa43816.csv")
Mcavernosa43816$fate <- factor(Mcavernosa43816$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa43816$title <- "Spondin 2a, extracellular matrix protein"

# ANOVA and Tukey's
Mcavernosa43816_stats <- aov(count~fate,data=Mcavernosa43816) %>%
  tukey_hsd()
Mcavernosa43816_stats

Mcavernosa43816_plot <-
  ggboxplot(
    Mcavernosa43816,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="salmon"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa43816_stats,label="p.adj.signif",y.position=2.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa43816_plot

# Mcavernosa54510 (Sequestering of TGFbeta in extracellular matrix)
Mcavernosa54510 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa54510.csv")
Mcavernosa54510$fate <- factor(Mcavernosa54510$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa54510$title <- "Sequestering of TGFbeta in extracellular matrix"

# ANOVA and Tukey's
Mcavernosa54510_stats <- aov(count~fate,data=Mcavernosa54510) %>%
  tukey_hsd()
Mcavernosa54510_stats

Mcavernosa54510_plot <-
  ggboxplot(
    Mcavernosa54510,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="salmon"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa54510_stats,label="p.adj.signif",y.position=1.9,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa54510_plot

# Mcavernosa47647 (Activation of NF-kappaB-inducing kinase activity)
Mcavernosa47647 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa47647.csv")
Mcavernosa47647$fate <- factor(Mcavernosa47647$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa47647$title <- "Activation of NF-kappaB-inducing kinase activity"

# ANOVA and Tukey's
Mcavernosa47647_stats <- aov(count~fate,data=Mcavernosa47647) %>%
  tukey_hsd()
Mcavernosa47647_stats

Mcavernosa47647_plot <-
  ggboxplot(
    Mcavernosa47647,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa47647_stats,label="p.adj.signif",y.position=2.6,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa47647_plot

# Mcavernosa50735 (Activation of NF-kappaB-inducing kinase activity)
Mcavernosa50735 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa50735.csv")
Mcavernosa50735$fate <- factor(Mcavernosa50735$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa50735$title <- "Activation of NF-kappaB-inducing kinase activity"

# ANOVA and Tukey's
Mcavernosa50735_stats <- aov(count~fate,data=Mcavernosa50735) %>%
  tukey_hsd()
Mcavernosa50735_stats

Mcavernosa50735_plot <-
  ggboxplot(
    Mcavernosa50735,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa50735_stats,label="p.adj.signif",y.position=2.2,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa50735_plot

# Mcavernosa61972 (Animal haem peroxidase)
Mcavernosa61972 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa61972.csv")
Mcavernosa61972$fate <- factor(Mcavernosa61972$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa61972$title <- "Animal haem peroxidase"

# ANOVA and Tukey's
Mcavernosa61972_stats <- aov(count~fate,data=Mcavernosa61972) %>%
  tukey_hsd()
Mcavernosa61972_stats

Mcavernosa61972_plot <-
  ggboxplot(
    Mcavernosa61972,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa61972_stats,label="p.adj.signif",y.position=1.1,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa61972_plot

# Mcavernosa14843 (Peroxidase activity)
Mcavernosa14843 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa14843.csv")
Mcavernosa14843$fate <- factor(Mcavernosa14843$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa14843$title <- "Peroxidase activity"

# ANOVA and Tukey's
Mcavernosa14843_stats <- aov(count~fate,data=Mcavernosa14843) %>%
  tukey_hsd()
Mcavernosa14843_stats

Mcavernosa14843_plot <-
  ggboxplot(
    Mcavernosa14843,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa14843_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa14843_plot

# Mcavernosa16729 (Peroxidase activity)
Mcavernosa16729 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa16729.csv")
Mcavernosa16729$fate <- factor(Mcavernosa16729$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa16729$title <- "Peroxidase activity"

# ANOVA and Tukey's
Mcavernosa16729_stats <- aov(count~fate,data=Mcavernosa16729) %>%
  tukey_hsd()
Mcavernosa16729_stats

Mcavernosa16729_plot <-
  ggboxplot(
    Mcavernosa16729,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa16729_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa16729_plot

# Mcavernosa184695 (Non-membrane spanning protein tyrosine kinase activity)
Mcavernosa184695 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa184695.csv")
Mcavernosa184695$fate <- factor(Mcavernosa184695$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa184695$title <- "Non-membrane spanning protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa184695_stats <- aov(count~fate,data=Mcavernosa184695) %>%
  tukey_hsd()
Mcavernosa184695_stats

Mcavernosa184695_plot <-
  ggboxplot(
    Mcavernosa184695,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa184695_stats,label="p.adj.signif",y.position=1.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa184695_plot

# Mcavernosa49548 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa49548 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa49548.csv")
Mcavernosa49548$fate <- factor(Mcavernosa49548$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa49548$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa49548_stats <- aov(count~fate,data=Mcavernosa49548) %>%
  tukey_hsd()
Mcavernosa49548_stats

Mcavernosa49548_plot <-
  ggboxplot(
    Mcavernosa49548,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa49548_stats,label="p.adj.signif",y.position=1.35,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa49548_plot

# Mcavernosa59808 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa59808 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa59808.csv")
Mcavernosa59808$fate <- factor(Mcavernosa59808$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa59808$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa59808_stats <- aov(count~fate,data=Mcavernosa59808) %>%
  tukey_hsd()
Mcavernosa59808_stats

Mcavernosa59808_plot <-
  ggboxplot(
    Mcavernosa59808,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
     add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa59808_stats,label="p.adj.signif",y.position=1.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa59808_plot

# Mcavernosa9650 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa9650 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa9650.csv")
Mcavernosa9650$fate <- factor(Mcavernosa9650$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa9650$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa9650_stats <- aov(count~fate,data=Mcavernosa9650) %>%
  tukey_hsd()
Mcavernosa9650_stats

Mcavernosa9650_plot <-
  ggboxplot(
    Mcavernosa9650,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
     add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa9650_stats,label="p.adj.signif",y.position=3,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa9650_plot

# Mcavernosa102943 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa102943 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa102943.csv")
Mcavernosa102943$fate <- factor(Mcavernosa102943$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa102943$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa102943_stats <- aov(count~fate,data=Mcavernosa102943) %>%
  tukey_hsd()
Mcavernosa102943_stats

Mcavernosa102943_plot <-
  ggboxplot(
    Mcavernosa102943,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa102943_stats,label="p.adj.signif",y.position=1.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa102943_plot

# Mcavernosa12949 (Non-membrane spanning protein tyrosine kinase activity)
Mcavernosa12949 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa12949.csv")
Mcavernosa12949$fate <- factor(Mcavernosa12949$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa12949$title <- "Non-membrane spanning protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa12949_stats <- aov(count~fate,data=Mcavernosa12949) %>%
  tukey_hsd()
Mcavernosa12949_stats

Mcavernosa12949_plot <-
  ggboxplot(
    Mcavernosa12949,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa12949_stats,label="p.adj.signif",y.position=2.55,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa12949_plot

# Mcavernosa29964 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa29964 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa29964.csv")
Mcavernosa29964$fate <- factor(Mcavernosa29964$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa29964$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa29964_stats <- aov(count~fate,data=Mcavernosa29964) %>%
  tukey_hsd()
Mcavernosa29964_stats

Mcavernosa29964_plot <-
  ggboxplot(
    Mcavernosa29964,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa29964_stats,label="p.adj.signif",y.position=2.85,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa29964_plot

# Mcavernosa20827 (Positive regulation of protein tyrosine kinase activity)
Mcavernosa20827 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa20827.csv")
Mcavernosa20827$fate <- factor(Mcavernosa20827$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa20827$title <- "Positive regulation of protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa20827_stats <- aov(count~fate,data=Mcavernosa20827) %>%
  tukey_hsd()
Mcavernosa20827_stats

Mcavernosa20827_plot <-
  ggboxplot(
    Mcavernosa20827,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa20827_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa20827_plot

# Mcavernosa71973 (Non-membrane spanning protein tyrosine kinase activity)
Mcavernosa71973 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa71973.csv")
Mcavernosa71973$fate <- factor(Mcavernosa71973$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa71973$title <- "Non-membrane spanning protein tyrosine kinase activity"

# ANOVA and Tukey's
Mcavernosa71973_stats <- aov(count~fate,data=Mcavernosa71973) %>%
  tukey_hsd()
Mcavernosa71973_stats

Mcavernosa71973_plot <-
  ggboxplot(
    Mcavernosa71973,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa71973_stats,label="p.adj.signif",y.position=1.85,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa71973_plot

# Mcavernosa126229 (Transforming growth factor-beta (TGF-beta) family)
Mcavernosa126229 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa126229.csv")
Mcavernosa126229$fate <- factor(Mcavernosa126229$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa126229$title <- "Transforming growth factor-beta (TGF-beta) family"

# ANOVA and Tukey's
Mcavernosa126229_stats <- aov(count~fate,data=Mcavernosa126229) %>%
  tukey_hsd()
Mcavernosa126229_stats

Mcavernosa126229_plot <-
  ggboxplot(
    Mcavernosa126229,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa126229_stats,label="p.adj.signif",y.position=2,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa126229_plot

# Mcavernosa27667 (Transforming growth factor-beta (TGF-beta) family)
Mcavernosa27667 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa27667.csv")
Mcavernosa27667$fate <- factor(Mcavernosa27667$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa27667$title <- "Transforming growth factor-beta (TGF-beta) family"

# ANOVA and Tukey's
Mcavernosa27667_stats <- aov(count~fate,data=Mcavernosa27667) %>%
  tukey_hsd()
Mcavernosa27667_stats

Mcavernosa27667_plot <-
  ggboxplot(
    Mcavernosa27667,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa27667_stats,label="p.adj.signif",y.position=1.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa27667_plot

# Mcavernosa21718 (Transmembrane serine threonine kinase forming with the TGF-beta type II serine threonine kinase receptor)
Mcavernosa21718 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa21718.csv")
Mcavernosa21718$fate <- factor(Mcavernosa21718$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa21718$title <- "Transmembrane serine threonine kinase"

# ANOVA and Tukey's
Mcavernosa21718_stats <- aov(count~fate,data=Mcavernosa21718) %>%
  tukey_hsd()
Mcavernosa21718_stats

Mcavernosa21718_plot <-
  ggboxplot(
    Mcavernosa21718,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa21718_stats,label="p.adj.signif",y.position=1.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa21718_plot

# Mcavernosa96261 (Transforming growth factor-beta (TGF-beta) family)
Mcavernosa96261 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa96261.csv")
Mcavernosa96261$fate <- factor(Mcavernosa96261$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa96261$title <- "Transforming growth factor-beta (TGF-beta) family"

# ANOVA and Tukey's
Mcavernosa96261_stats <- aov(count~fate,data=Mcavernosa96261) %>%
  tukey_hsd()
Mcavernosa96261_stats

Mcavernosa96261_plot <-
  ggboxplot(
    Mcavernosa96261,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa96261_stats,label="p.adj.signif",y.position=1.85,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa96261_plot

# creating legend for large plot
legend_mcav <- get_legend(Mcavernosa96261_plot)

Mcavernosa96261_plot <- Mcavernosa96261_plot + theme(legend.position = "none")

# Mcavernosa366959 (F-box WD repeat-containing protein 8)
Mcavernosa366959 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa366959.csv")
Mcavernosa366959$fate <- factor(Mcavernosa366959$fate, levels = c("healthy", "nai", "diseased"))
Mcavernosa366959$title <- "F-box WD repeat-containing protein 8"

# ANOVA and Tukey's
Mcavernosa366959_stats <- aov(count~fate,data=Mcavernosa366959) %>%
  tukey_hsd()
Mcavernosa366959_stats

Mcavernosa366959_plot <-
  ggboxplot(
    Mcavernosa366959,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="lightskyblue"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa366959_stats,label="p.adj.signif",y.position=0.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa366959_plot


#### MULTIPLOT MCAV ####

# arranging all plots
transmission_mcav<-ggarrange(Mcavernosa64646_plot,
                          Mcavernosa69648_plot,
                          Mcavernosa9810_plot,
                          Mcavernosa43816_plot,
                          Mcavernosa54510_plot,
                          Mcavernosa50735_plot,
                          Mcavernosa61972_plot,
                          Mcavernosa14843_plot,
                          Mcavernosa16729_plot,
                          Mcavernosa184695_plot,
                          Mcavernosa49548_plot,
                          Mcavernosa59808_plot,
                          Mcavernosa9650_plot,
                          Mcavernosa12949_plot,
                          Mcavernosa29964_plot,
                          Mcavernosa20827_plot,
                          Mcavernosa126229_plot,
                          Mcavernosa21718_plot,
                          Mcavernosa96261_plot,
                          legend_mcav,
                          heights = c(4,4,4,4),
                           widths = c(5,5,5,5,5),
                           ncol = 5,
                           nrow = 4)
transmission_mcav
ggsave("transmission_mcav.pdf", transmission_mcav, width=20, height=16,dpi = 300)


#### BOXPLOTS OFAV ####

# non-matching DEGs by species
# Ofaveolata266584 (Sequestering of TGFbeta in extracellular matrix)
Ofaveolata266584 <- read.csv(file = "../DESeq2/ofav/Ofaveolata266584.csv")
colnames(Ofaveolata266584)[3] = "fate"
Ofaveolata266584$fate <- str_replace(Ofaveolata266584$fate, "control", "healthy")
Ofaveolata266584$fate <- str_replace(Ofaveolata266584$fate, "sctld", "diseased")

Ofaveolata266584$fate <- factor(Ofaveolata266584$fate, levels = c("healthy", "diseased"))
Ofaveolata266584$title <- "Sequestering of TGFbeta in extracellular matrix"

# ANOVA and Tukey's
Ofaveolata266584_stats <- aov(count~fate,data=Ofaveolata266584) %>%
  tukey_hsd()
Ofaveolata266584_stats

Ofaveolata266584_plot <-
  ggboxplot(
    Ofaveolata266584,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="salmon"), strip.text = element_text(size=9), 
        axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata266584_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata266584_plot

# Ofaveolata264011 (Negative regulation of extracellular matrix disassembly)
Ofaveolata264011 <- read.csv(file = "../DESeq2/ofav/Ofaveolata264011.csv")
colnames(Ofaveolata264011)[3] = "fate"
Ofaveolata264011$fate <- str_replace(Ofaveolata264011$fate, "control", "healthy")
Ofaveolata264011$fate <- str_replace(Ofaveolata264011$fate, "sctld", "diseased")

Ofaveolata264011$fate <- factor(Ofaveolata264011$fate, levels = c("healthy", "diseased"))
Ofaveolata264011$title <- "Negative regulation of extracellular matrix disassembly"

# ANOVA and Tukey's
Ofaveolata264011_stats <- aov(count~fate,data=Ofaveolata264011) %>%
  tukey_hsd()
Ofaveolata264011_stats

Ofaveolata264011_plot <-
  ggboxplot(
    Ofaveolata264011,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="salmon"), strip.text = element_text(size=9)
        , axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata264011_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata264011_plot

# Ofaveolata261513 (Activation of NF-kappaB-inducing kinase activity)
Ofaveolata261513 <- read.csv(file = "../DESeq2/ofav/Ofaveolata261513.csv")
colnames(Ofaveolata261513)[3] = "fate"
Ofaveolata261513$fate <- str_replace(Ofaveolata261513$fate, "control", "healthy")
Ofaveolata261513$fate <- str_replace(Ofaveolata261513$fate, "sctld", "diseased")

Ofaveolata261513$fate <- factor(Ofaveolata261513$fate, levels = c("healthy", "diseased"))
Ofaveolata261513$title <- "Activation of NF-kappaB-inducing kinase activity"

# ANOVA and Tukey's
Ofaveolata261513_stats <- aov(count~fate,data=Ofaveolata261513) %>%
  tukey_hsd()
Ofaveolata261513_stats

Ofaveolata261513_plot <-
  ggboxplot(
    Ofaveolata261513,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata261513_stats,label="p.adj.signif",y.position=3.1,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata261513_plot

# Ofaveolata262409 (Cytoplasmic sequestering of NF-kappaB)
Ofaveolata262409 <- read.csv(file = "../DESeq2/ofav/Ofaveolata262409.csv")
colnames(Ofaveolata262409)[3] = "fate"
Ofaveolata262409$fate <- str_replace(Ofaveolata262409$fate, "control", "healthy")
Ofaveolata262409$fate <- str_replace(Ofaveolata262409$fate, "sctld", "diseased")

Ofaveolata262409$fate <- factor(Ofaveolata262409$fate, levels = c("healthy", "diseased"))
Ofaveolata262409$title <- "Cytoplasmic sequestering of NF-kappaB"

# ANOVA and Tukey's
Ofaveolata262409_stats <- aov(count~fate,data=Ofaveolata262409) %>%
  tukey_hsd()
Ofaveolata262409_stats

Ofaveolata262409_plot <-
  ggboxplot(
    Ofaveolata262409,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata262409_stats,label="p.adj.signif",y.position=2.3,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata262409_plot

# Ofaveolata253694 (Animal haem peroxidase)
Ofaveolata253694 <- read.csv(file = "../DESeq2/ofav/Ofaveolata253694.csv")
colnames(Ofaveolata253694)[3] = "fate"
Ofaveolata253694$fate <- str_replace(Ofaveolata253694$fate, "control", "healthy")
Ofaveolata253694$fate <- str_replace(Ofaveolata253694$fate, "sctld", "diseased")

Ofaveolata253694$fate <- factor(Ofaveolata253694$fate, levels = c("healthy", "diseased"))
Ofaveolata253694$title <- "Animal haem peroxidase"

# ANOVA and Tukey's
Ofaveolata253694_stats <- aov(count~fate,data=Ofaveolata253694) %>%
  tukey_hsd()
Ofaveolata253694_stats

Ofaveolata253694_plot <-
  ggboxplot(
    Ofaveolata253694,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata253694_stats,label="p.adj.signif",y.position=1.1,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata253694_plot

# Ofaveolata261831 (Peroxidase activity)
Ofaveolata261831 <- read.csv(file = "../DESeq2/ofav/Ofaveolata261831.csv")
colnames(Ofaveolata261831)[3] = "fate"
Ofaveolata261831$fate <- str_replace(Ofaveolata261831$fate, "control", "healthy")
Ofaveolata261831$fate <- str_replace(Ofaveolata261831$fate, "sctld", "diseased")

Ofaveolata261831$fate <- factor(Ofaveolata261831$fate, levels = c("healthy", "diseased"))
Ofaveolata261831$title <- "Peroxidase activity"

# ANOVA and Tukey's
Ofaveolata261831_stats <- aov(count~fate,data=Ofaveolata261831) %>%
  tukey_hsd()
Ofaveolata261831_stats

Ofaveolata261831_plot <-
  ggboxplot(
    Ofaveolata261831,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata261831_stats,label="p.adj.signif",y.position=2,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata261831_plot

# Ofaveolata262557 (Peroxidase activity)
Ofaveolata262557 <- read.csv(file = "../DESeq2/ofav/Ofaveolata262557.csv")
colnames(Ofaveolata262557)[3] = "fate"
Ofaveolata262557$fate <- str_replace(Ofaveolata262557$fate, "control", "healthy")
Ofaveolata262557$fate <- str_replace(Ofaveolata262557$fate, "sctld", "diseased")

Ofaveolata262557$fate <- factor(Ofaveolata262557$fate, levels = c("healthy", "diseased"))
Ofaveolata262557$title <- "Peroxidase activity"

# ANOVA and Tukey's
Ofaveolata262557_stats <- aov(count~fate,data=Ofaveolata262557) %>%
  tukey_hsd()
Ofaveolata262557_stats

Ofaveolata262557_plot <-
  ggboxplot(
    Ofaveolata262557,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata262557_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata262557_plot

# Ofaveolata251117 (Peroxidase activity)
Ofaveolata251117 <- read.csv(file = "../DESeq2/ofav/Ofaveolata251117.csv")
colnames(Ofaveolata251117)[3] = "fate"
Ofaveolata251117$fate <- str_replace(Ofaveolata251117$fate, "control", "healthy")
Ofaveolata251117$fate <- str_replace(Ofaveolata251117$fate, "sctld", "diseased")

Ofaveolata251117$fate <- factor(Ofaveolata251117$fate, levels = c("healthy", "diseased"))
Ofaveolata251117$title <- "Peroxidase activity"

# ANOVA and Tukey's
Ofaveolata251117_stats <- aov(count~fate,data=Ofaveolata251117) %>%
  tukey_hsd()
Ofaveolata251117_stats

Ofaveolata251117_plot <-
  ggboxplot(
    Ofaveolata251117,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata251117_stats,label="p.adj.signif",y.position=1.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata251117_plot

# Ofaveolata262519 (Peroxidase activity)
Ofaveolata262519 <- read.csv(file = "../DESeq2/ofav/Ofaveolata262519.csv")
colnames(Ofaveolata262519)[3] = "fate"
Ofaveolata262519$fate <- str_replace(Ofaveolata262519$fate, "control", "healthy")
Ofaveolata262519$fate <- str_replace(Ofaveolata262519$fate, "sctld", "diseased")

Ofaveolata262519$fate <- factor(Ofaveolata262519$fate, levels = c("healthy", "diseased"))
Ofaveolata262519$title <- "Peroxidase activity"

# ANOVA and Tukey's
Ofaveolata262519_stats <- aov(count~fate,data=Ofaveolata262519) %>%
  tukey_hsd()
Ofaveolata262519_stats

Ofaveolata262519_plot <-
  ggboxplot(
    Ofaveolata262519,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata262519_stats,label="p.adj.signif",y.position=1.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata262519_plot

# Ofaveolata240196 (Peroxidase activity)
Ofaveolata240196 <- read.csv(file = "../DESeq2/ofav/Ofaveolata240196.csv")
colnames(Ofaveolata240196)[3] = "fate"
Ofaveolata240196$fate <- str_replace(Ofaveolata240196$fate, "control", "healthy")
Ofaveolata240196$fate <- str_replace(Ofaveolata240196$fate, "sctld", "diseased")

Ofaveolata240196$fate <- factor(Ofaveolata240196$fate, levels = c("healthy", "diseased"))
Ofaveolata240196$title <- "Peroxidase activity"

# ANOVA and Tukey's
Ofaveolata240196_stats <- aov(count~fate,data=Ofaveolata240196) %>%
  tukey_hsd()
Ofaveolata240196_stats

Ofaveolata240196_plot <-
  ggboxplot(
    Ofaveolata240196,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata240196_stats,label="p.adj.signif",y.position=1.6,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata240196_plot

# Ofaveolata248576 (Animal haem peroxidase)
Ofaveolata248576 <- read.csv(file = "../DESeq2/ofav/Ofaveolata248576.csv")
colnames(Ofaveolata248576)[3] = "fate"
Ofaveolata248576$fate <- str_replace(Ofaveolata248576$fate, "control", "healthy")
Ofaveolata248576$fate <- str_replace(Ofaveolata248576$fate, "sctld", "diseased")

Ofaveolata248576$fate <- factor(Ofaveolata248576$fate, levels = c("healthy", "diseased"))
Ofaveolata248576$title <- "Animal haem peroxidase"

# ANOVA and Tukey's
Ofaveolata248576_stats <- aov(count~fate,data=Ofaveolata248576) %>%
  tukey_hsd()
Ofaveolata248576_stats

Ofaveolata248576_plot <-
  ggboxplot(
    Ofaveolata248576,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata248576_stats,label="p.adj.signif",y.position=1.8,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata248576_plot

# Ofaveolata243769 (Glutathione peroxidase activity)
Ofaveolata243769 <- read.csv(file = "../DESeq2/ofav/Ofaveolata243769.csv")
colnames(Ofaveolata243769)[3] = "fate"
Ofaveolata243769$fate <- str_replace(Ofaveolata243769$fate, "control", "healthy")
Ofaveolata243769$fate <- str_replace(Ofaveolata243769$fate, "sctld", "diseased")

Ofaveolata243769$fate <- factor(Ofaveolata243769$fate, levels = c("healthy", "diseased"))
Ofaveolata243769$title <- "Glutathione peroxidase activity"

# ANOVA and Tukey's
Ofaveolata243769_stats <- aov(count~fate,data=Ofaveolata243769) %>%
  tukey_hsd()
Ofaveolata243769_stats

Ofaveolata243769_plot <-
  ggboxplot(
    Ofaveolata243769,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata243769_stats,label="p.adj.signif",y.position=2.3,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata243769_plot

# Ofaveolata256330 (Glutathione peroxidase activity)
Ofaveolata256330 <- read.csv(file = "../DESeq2/ofav/Ofaveolata256330.csv")
colnames(Ofaveolata256330)[3] = "fate"
Ofaveolata256330$fate <- str_replace(Ofaveolata256330$fate, "control", "healthy")
Ofaveolata256330$fate <- str_replace(Ofaveolata256330$fate, "sctld", "diseased")

Ofaveolata256330$fate <- factor(Ofaveolata256330$fate, levels = c("healthy", "diseased"))
Ofaveolata256330$title <- "Glutathione peroxidase activity"

# ANOVA and Tukey's
Ofaveolata256330_stats <- aov(count~fate,data=Ofaveolata256330) %>%
  tukey_hsd()
Ofaveolata256330_stats

Ofaveolata256330_plot <-
  ggboxplot(
    Ofaveolata256330,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata256330_stats,label="p.adj.signif",y.position=2.9,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata256330_plot

# Ofaveolata254093 (Transmembrane receptor protein tyrosine kinase activity)
Ofaveolata254093 <- read.csv(file = "../DESeq2/ofav/Ofaveolata254093.csv")
colnames(Ofaveolata254093)[3] = "fate"
Ofaveolata254093$fate <- str_replace(Ofaveolata254093$fate, "control", "healthy")
Ofaveolata254093$fate <- str_replace(Ofaveolata254093$fate, "sctld", "diseased")

Ofaveolata254093$fate <- factor(Ofaveolata254093$fate, levels = c("healthy", "diseased"))
Ofaveolata254093$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Ofaveolata254093_stats <- aov(count~fate,data=Ofaveolata254093) %>%
  tukey_hsd()
Ofaveolata254093_stats

Ofaveolata254093_plot <-
  ggboxplot(
    Ofaveolata254093,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata254093_stats,label="p.adj.signif",y.position=2.1,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata254093_plot

# Ofaveolata267705 (Transmembrane receptor protein tyrosine kinase activity)
Ofaveolata267705 <- read.csv(file = "../DESeq2/ofav/Ofaveolata267705.csv")
colnames(Ofaveolata267705)[3] = "fate"
Ofaveolata267705$fate <- str_replace(Ofaveolata267705$fate, "control", "healthy")
Ofaveolata267705$fate <- str_replace(Ofaveolata267705$fate, "sctld", "diseased")

Ofaveolata267705$fate <- factor(Ofaveolata267705$fate, levels = c("healthy", "diseased"))
Ofaveolata267705$title <- "Transmembrane receptor protein tyrosine kinase activity"

# ANOVA and Tukey's
Ofaveolata267705_stats <- aov(count~fate,data=Ofaveolata267705) %>%
  tukey_hsd()
Ofaveolata267705_stats

Ofaveolata267705_plot <-
  ggboxplot(
    Ofaveolata267705,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata267705_stats,label="p.adj.signif",y.position=2.3,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata267705_plot

# Ofaveolata268849 (Non-membrane spanning protein tyrosine kinase activity)
Ofaveolata268849 <- read.csv(file = "../DESeq2/ofav/Ofaveolata268849.csv")
colnames(Ofaveolata268849)[3] = "fate"
Ofaveolata268849$fate <- str_replace(Ofaveolata268849$fate, "control", "healthy")
Ofaveolata268849$fate <- str_replace(Ofaveolata268849$fate, "sctld", "diseased")

Ofaveolata268849$fate <- factor(Ofaveolata268849$fate, levels = c("healthy", "diseased"))
Ofaveolata268849$title <- "Non-membrane spanning protein tyrosine kinase activity"

# ANOVA and Tukey's
Ofaveolata268849_stats <- aov(count~fate,data=Ofaveolata268849) %>%
  tukey_hsd()
Ofaveolata268849_stats

Ofaveolata268849_plot <-
  ggboxplot(
    Ofaveolata268849,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata268849_stats,label="p.adj.signif",y.position=2,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata268849_plot

# Ofaveolata265457 (WD repeat-containing protein 63)
Ofaveolata265457 <- read.csv(file = "../DESeq2/ofav/Ofaveolata265457.csv")
colnames(Ofaveolata265457)[3] = "fate"
Ofaveolata265457$fate <- str_replace(Ofaveolata265457$fate, "control", "healthy")
Ofaveolata265457$fate <- str_replace(Ofaveolata265457$fate, "sctld", "diseased")

Ofaveolata265457$fate <- factor(Ofaveolata265457$fate, levels = c("healthy", "diseased"))
Ofaveolata265457$title <- "WD repeat-containing protein 63"

# ANOVA and Tukey's
Ofaveolata265457_stats <- aov(count~fate,data=Ofaveolata265457) %>%
  tukey_hsd()
Ofaveolata265457_stats

Ofaveolata265457_plot <-
  ggboxplot(
    Ofaveolata265457,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right", strip.background = element_rect(fill="lightskyblue"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata265457_stats,label="p.adj.signif",y.position=1.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Ofaveolata265457_plot

# creating legend for large plot
legend_ofav <- get_legend(Ofaveolata265457_plot)

Ofaveolata265457_plot <- Ofaveolata265457_plot + theme(legend.position = "none")


#### MULTIPLOT OFAV ####

# arranging all plots
transmission_ofav<-ggarrange(Ofaveolata264011_plot,
                          Ofaveolata261513_plot,
                          Ofaveolata262409_plot,
                          Ofaveolata253694_plot,
                          NULL,
                          Ofaveolata262557_plot,
                          Ofaveolata251117_plot,
                          Ofaveolata262519_plot,
                          Ofaveolata240196_plot,
                          legend_ofav,
                          Ofaveolata248576_plot,
                          Ofaveolata256330_plot,
                          Ofaveolata267705_plot,
                          Ofaveolata265457_plot,
                          NULL,
                          heights = c(4,4,4),
                          widths = c(5,5,5,5,1),
                          ncol = 5,
                          nrow = 3)
transmission_ofav
ggsave("transmission_ofav.pdf", transmission_ofav, width=16.8, height=15,dpi = 300)
