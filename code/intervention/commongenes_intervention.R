#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### DESEQ IMPORT ####

trans_diseased_healthy_lpv <- read.csv(file = "../../transmission/DESeq2/mcav2015/diseased_nai_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_dh" = lpv)

inter_diseased0_healthy0_lpv <- read.csv(file = "../DESeq2/mcav2015/diseased0_healthy0_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_d0h0" = lpv)

inter_treated1_healthy1_lpv <- read.csv(file = "../DESeq2/mcav2015/treated1_healthy1_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_t1h1" = lpv)

inter_treated1_diseased0_lpv <- read.csv(file = "../DESeq2/mcav2015/treated1_diseased0_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_t1d0" = lpv)


#### DEG MATCHING ####

# These sections of code do several things: 1) join common DEGs across experiments with -log10(pval), 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3)  adds mcav gene annotations, and 4) then pulls on corresponding KOG classes

# diseased vs healthy for both experiments
trans_diseased_healthy_lpv %>%
  inner_join(inter_diseased0_healthy0_lpv, by = "gene") %>%
  filter(abs(lpv_dh) >= 1 & abs(lpv_d0h0) >= 1) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> diseased_healthy

# treated vs healthy for intervention experiment
trans_diseased_healthy_lpv %>%
  inner_join(inter_treated1_healthy1_lpv, by = "gene") %>%
  filter(abs(lpv_dh) >= 1 & abs(lpv_t1h1) >= 1) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> treated_healthy

# treated vs diseased for intervention experiment
trans_diseased_healthy_lpv %>%
  inner_join(inter_treated1_diseased0_lpv, by = "gene") %>%
  filter(abs(lpv_dh) >= 1 & abs(lpv_t1d0) >= 1) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> treated_diseased


#### KOG MATCHING ####

# filtering and summarizing DEGs by KOG class for high-level comparisons between experiments
diseased_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_up" = n) -> KOG_trans_up

diseased_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_down" = n) -> KOG_trans_down

diseased_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_d0h0 >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "dh_up" = n) -> KOG_dh_up

diseased_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_d0h0 <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "dh_down" = n) -> KOG_dh_down

# joining all KOG class sums in a single dataframe
KOG_trans_up %>%
  inner_join(KOG_dh_up, by = "KOG") %>%
  inner_join(KOG_trans_down, by = "KOG") %>%
  inner_join(KOG_dh_down, by = "KOG") -> KOG_dh_match

# melting dataframe for plotting
KOG_dh_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_dh_melt

# creating a custom color palette
colorCount = length(unique(KOG_dh_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_dh_sum <- ggplot(KOG_dh_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Proportion of DEGs") +
  theme_classic()
KOG_dh_sum
ggsave("orthofinder KOG dh abundance.pdf", plot= KOG_dh_sum, width=8, height=6, units="in", dpi=300)


#### KOG MATCH TREATED ####

# filtering and summarizing DEGs by KOG class for high-level comparisons between experiments
treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_up" = n) -> KOG_trans_up2

treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_down" = n) -> KOG_trans_down2

treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_t1h1 >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "th_up" = n) -> KOG_th_up

treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_t1h1 <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "th_down" = n) -> KOG_th_down

# joining all KOG class sums in a single dataframe
KOG_trans_up2 %>%
  inner_join(KOG_th_up, by = "KOG") %>%
  inner_join(KOG_trans_down2, by = "KOG") %>%
  inner_join(KOG_th_down, by = "KOG") -> KOG_th_match

# melting dataframe for plotting
KOG_th_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_th_melt

# creating a custom color palette
colorCount = length(unique(KOG_th_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_th_sum <- ggplot(KOG_th_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Proportion of DEGs") +
  theme_classic()
KOG_th_sum
ggsave("orthofinder KOG th abundance.pdf", plot= KOG_th_sum, width=8, height=6, units="in", dpi=300)


#### KOG MATCH TREATED 2 ####

# filtering and summarizing DEGs by KOG class for high-level comparisons between experiments
treated_diseased %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_up" = n) -> KOG_trans_up3

treated_diseased %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_down" = n) -> KOG_trans_down3

treated_diseased %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_t1d0 >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "td_up" = n) -> KOG_td_up

treated_diseased %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_t1d0 <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "td_down" = n) -> KOG_td_down

# joining all KOG class sums in a single dataframe
KOG_trans_up3 %>%
  inner_join(KOG_td_up, by = "KOG") %>%
  inner_join(KOG_trans_down3, by = "KOG") %>%
  inner_join(KOG_td_down, by = "KOG") -> KOG_td_match

# melting dataframe for plotting
KOG_td_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_td_melt

# creating a custom color palette
colorCount = length(unique(KOG_td_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_td_sum <- ggplot(KOG_td_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Proportion of DEGs") +
  theme_classic()
KOG_td_sum
ggsave("orthofinder KOG td abundance.pdf", plot= KOG_td_sum, width=8, height=6, units="in", dpi=300)


#### VENN DIAGRAMS ####

# first creating a set of up/downregulated DEGs by diseased vs healthy for each experiment
diseased_healthy %>%
  filter(lpv_dh >= 1) %>%
  pull(gene) -> trans_up

diseased_healthy %>%
  filter(lpv_dh <= -1) %>%
  pull(gene) -> trans_down

diseased_healthy %>%
  filter(lpv_d0h0 >= 1) %>%
  pull(gene) -> dis_up

diseased_healthy %>%
  filter(lpv_d0h0 <= -1) %>%
  pull(gene) -> dis_down

# then creating a second set for treated vs healthy
treated_healthy %>%
  filter(lpv_dh >= 1) %>%
  pull(gene) -> trans_up2

treated_healthy %>%
  filter(lpv_dh <= -1) %>%
  pull(gene) -> trans_down2

treated_healthy %>%
  filter(lpv_t1h1 >= 1) %>%
  pull(gene) -> treat_up

treated_healthy %>%
  filter(lpv_t1h1 <= -1) %>%
  pull(gene) -> treat_down

# then creating a final set for treated vs diseased
treated_diseased %>%
  filter(lpv_dh >= 1) %>%
  pull(gene) -> trans_up3

treated_diseased %>%
  filter(lpv_dh <= -1) %>%
  pull(gene) -> trans_down3

treated_diseased %>%
  filter(lpv_t1d0 >= 1) %>%
  pull(gene) -> treat_up2

treated_diseased %>%
  filter(lpv_t1d0 <= -1) %>%
  pull(gene) -> treat_down2

# diseased vs healthy by experiment
venn_dh=h=venn.diagram(
  x = list("Trans up"=trans_up, "Trans down"=trans_down,"Inter up"=dis_up, "Inter down"=dis_down),
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
pdf(file="Venn_diseased_healthy.pdf", height=10, width=12)
grid.draw(venn_dh)
dev.off()

# treated vs healthy by experiment
venn_th=h=venn.diagram(
  x = list("Trans up"=trans_up2, "Trans down"=trans_down2,"Inter up"=treat_up, "Inter down"=treat_down),
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
pdf(file="Venn_treated_healthy.pdf", height=10, width=12)
grid.draw(venn_th)
dev.off()

# treated vs healthy by experiment
venn_td=h=venn.diagram(
  x = list("Trans up"=trans_up3, "Trans down"=trans_down3,"Inter up"=treat_up2, "Inter down"=treat_down2),
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
pdf(file="Venn_treated_diseased.pdf", height=10, width=12)
grid.draw(venn_td)
dev.off()


#### COMBINE TREATMENTS ####

# diseased vs healthy vs treated for both experiments
trans_diseased_healthy_lpv %>%
  inner_join(inter_diseased0_healthy0_lpv, by = "gene") %>%
  inner_join(inter_treated1_healthy1_lpv, by = "gene") %>%
  filter(abs(lpv_dh) >= 1 & abs(lpv_d0h0) >= 1 & abs(lpv_t1h1) >= 1) %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     gene_name = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>%
  left_join(read.table(file = "../../annotate/mcav2015/Mcavernosa2015_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> diseased_treated_healthy


# first creating a column of combined gene names from both species, then removing unannotated genes
diseased_treated_healthy %>%
  mutate(gene_name = str_replace(gene_name, "NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "-","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  diseased_treated_healthy


#### KOG MATCH ALL TREATMENTS ####

# filtering and summarizing DEGs by KOG class for high-level comparisons between experiments
diseased_treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_up" = n) -> KOG_trans_up4

diseased_treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_dh <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "trans_down" = n) -> KOG_trans_down4

diseased_treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_d0h0 >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "dh_up" = n) -> KOG_dh_up2

diseased_treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_d0h0 <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "dh_down" = n) -> KOG_dh_down2

diseased_treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_t1h1 >= 1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "th_up" = n) -> KOG_th_up2

diseased_treated_healthy %>%
  mutate(KOG_mcav = replace(KOG_mcav, KOG_mcav == "", NA)) %>%
  filter(lpv_t1h1 <= -1) %>%
  count(KOG_mcav) %>%
  rename("KOG" = KOG_mcav, "th_down" = n) -> KOG_th_down2

# joining all KOG class sums in a single dataframe
KOG_trans_up4 %>%
  inner_join(KOG_dh_up2, by = "KOG") %>%
  inner_join(KOG_th_up2, by = "KOG") %>%
  inner_join(KOG_trans_down4, by = "KOG") %>%
  inner_join(KOG_dh_down2, by = "KOG") %>%
  inner_join(KOG_th_down2, by = "KOG") -> KOG_all_match

# melting dataframe for plotting
KOG_all_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_all_melt

# creating a custom color palette
colorCount = length(unique(KOG_all_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_all_sum <- ggplot(KOG_all_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Proportion of DEGs") +
  theme_classic()
KOG_all_sum
ggsave("orthofinder KOG all abundance.pdf", plot= KOG_all_sum, width=10, height=6, units="in", dpi=300)


#### VSD by EXPERIMENT/TREATMENT ####

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup IDs, then removing NAI samples
load("../../transmission/DESeq2/mcav2015/vsd.RData")
design %>%
  unite("full_id", id,genotype,fate, sep = "-", remove = FALSE) %>%
  select(full_id, id, fate) %>%
  mutate(experiment = "Transmission") -> design_trans
vsd_trans <- subset(vsd, rownames(vsd) %in% diseased_treated_healthy$gene)
vsd_trans2 <- subset(vsd, rownames(vsd) %in% diseased_healthy$gene)

vsd_trans %>%
  as.data.frame() %>%
  select(-contains("nai")) %>%
  as.matrix() -> vsd_trans

vsd_trans2 %>%
  as.data.frame() %>%
  select(-contains("nai")) %>%
  as.matrix() -> vsd_trans2

load("../DESeq2/mcav2015/vsd.RData")
design %>%
  unite("full_id", id,treatment.time,genotype, sep = "-", remove = FALSE) %>%
  select(full_id, id, fate) %>%  
  mutate(experiment = "Intervention") -> design_int
vsd_int <- subset(vsd, rownames(vsd) %in% diseased_treated_healthy$gene)
vsd_int2 <- subset(vsd, rownames(vsd) %in% diseased_healthy$gene)

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_trans, vsd_int)
vsd_comb2 <- cbind(vsd_trans2, vsd_int2)
design_comb <- rbind(design_trans, design_int)
design_comb$id <- gsub("-",".", design_comb$id)
design_comb$full_id <- gsub("-",".", design_comb$full_id)

# removing NAI samples and reordering design dataframe for plotting
design_comb <- design_comb[(design_comb$fate != "nai"),]
design_comb$experiment_fate <- paste(design_comb$experiment,design_comb$fate,sep=".")
design_comb$experiment_fate <- factor(design_comb$experiment_fate, levels = c("Transmission.healthy","Transmission.diseased","Intervention.healthy","Intervention.diseased","Intervention.treated"))
design_comb <- design_comb[order(design_comb$experiment_fate),]
design_comb$label <- paste(design_comb$id,design_comb$fate,sep=".")

# reordering counts matrix according to design dataframe
head(vsd_comb)
design_order <- design_comb$full_id
vsd_comb <- vsd_comb[,design_order]
colnames(vsd_comb) <- design_comb$label
head(vsd_comb)

head(vsd_comb2)
vsd_comb2 <- vsd_comb2[,design_order]
colnames(vsd_comb2) <- design_comb$label
head(vsd_comb2)


#### HEATMAPS ####

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene to gene annotations
gene_names <- as.data.frame(cbind(diseased_treated_healthy$gene, diseased_treated_healthy$gene_name))
gene_names2 <- as.data.frame(cbind(diseased_healthy$gene, diseased_healthy$annot_mcav))

# save(orthologs, diseased_healthy, diseased_nai, nai_healthy, mcav_up, mcav_down, ofav_up, ofav_down, design_ofav, design_mcav, design_comb, vsd_ofav, vsd_mcav, vsd_comb, file = "orthofinder_DEGs_species.RData")
# load("orthofinder_DEGs_species.RData")

# heatmap of original vsd relationships (separated by experiment/treatment)

# p < 0.1 (all genes)
pdf(file="commongenes_heatmap_p0.1.pdf", height=8, width=24)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(diseased_treated_healthy$lpv_dh)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-1, 
           # sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           sort=colnames(vsd_comb),
           cex=0.8,
           pdf=F,
)
dev.off()

# just diseased vs healthy
pdf(file="commongenes_dh_heatmap_p0.1.pdf", height=48, width=26)
uniHeatmap(vsd=vsd_comb2,gene.names=gene_names2,
           metric=-(abs(diseased_healthy$lpv_dh)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-1, 
           # sort=c(1:ncol(vsd_comb2)), # overrides sorting of columns according to hierarchical clustering
           sort=colnames(vsd_comb2),
           cex=0.8,
           pdf=F,
)
dev.off()
