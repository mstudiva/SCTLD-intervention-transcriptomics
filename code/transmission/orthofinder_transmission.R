#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)

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

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) then add ofav and mcav gene annotations

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
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> nai_healthy
nai_healthy$Orthogroup <- make.unique(nai_healthy$Orthogroup, sep = "_") 


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

# loading gene counts for each species, then replacing species-specific gene IDs with orthogroup IDs, then removing NAI samples
load("../DESeq2/ofav/initial.RData")
counts_ofav <- subset(countData, rownames(countData) %in% diseased_healthy$Protein_ofav)
design_ofav <- design 

load("../DESeq2/mcav2015/initial.RData")
counts_mcav <- subset(countData, rownames(countData) %in% diseased_healthy$Protein_mcav)
design_mcav <- design 

# have to detach DESeq2 because it contains some functions with the same name as dplyr
detach("package:DESeq2", unload=TRUE)
library(tidyverse)

# combining the design data frames
design_comb <- rbind(design_mcav, design_ofav)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$species,design_comb$fate,sep=".")

# replacing gene names with orthogroups
counts_ofav %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% diseased_healthy$Protein_ofav, diseased_healthy$Orthogroup, diseased_healthy$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") -> counts_ofav

counts_mcav %>%
  rownames_to_column(var = "Protein_mcav") %>%
  mutate(Protein_mcav = if_else(Protein_mcav %in% diseased_healthy$Protein_mcav, diseased_healthy$Orthogroup, diseased_healthy$Orthogroup)) %>%
  column_to_rownames(var = "Protein_mcav") -> counts_mcav

# combining the two count data frames, then replacing sample IDs with full IDs (ID + factors)
counts_comb <- cbind(counts_mcav, counts_ofav)
colnames(counts_comb) = design_comb$full_id

# removing NAI samples and reordering design dataframe for plotting
design_comb <- design_comb[(design_comb$fate != "nai"),]
design_comb$species_fate <- paste(design_comb$species,design_comb$fate,sep=".")
design_comb$species_fate <- factor(design_comb$species_fate, levels = c("Mcav.healthy","Mcav.diseased","Ofav.healthy","Ofav.diseased"))
design_comb <- design_comb[order(design_comb$species_fate),]

# reordering counts matrix according to design dataframe
counts_comb %>%
  select(-contains("nai")) -> counts_comb
head(counts_comb)
counts_comb<- counts_comb[,order(design_comb$full_id)]
tail(counts_comb)

library(DESeq2)

# making a combined DESeq model
rownames(design_comb) = design_comb$full_id
dds = DESeqDataSetFromMatrix(countData=counts_comb, colData=design_comb, design=~ genotype+fate)

# creating variance stabilized array of counts for heatmap
Vsd=varianceStabilizingTransformation(dds)
vsd=assay(Vsd)

save(orthologs, diseased_healthy, diseased_nai, nai_healthy, mcav_up, mcav_down, ofav_up, ofav_down, counts_ofav, design_ofav, counts_mcav, design_mcav, design_comb, counts_comb, vsd, file = "orthofinder_DEGs.RData")


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
