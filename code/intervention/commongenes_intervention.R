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

# p < 0.0001
pdf(file="commongenes_dh_heatmap_p0.0001.pdf", height=5, width=15)
uniHeatmap(vsd=vsd_comb2,gene.names=gene_names2,
           metric=-(abs(diseased_healthy$lpv_dh)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-4, 
           # sort=c(1:ncol(vsd_comb2)), # overrides sorting of columns according to hierarchical clustering
           sort=colnames(vsd_comb2),
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 1e-6
pdf(file="commongenes_dh_heatmap_p1e6.pdf", height=1.75, width=15)
uniHeatmap(vsd=vsd_comb2,gene.names=gene_names2,
           metric=-(abs(diseased_healthy$lpv_dh)), # metric of gene significance
           # metric2=-(abs(diseased_healthy$lpv_ofav)),
           cutoff=-6, 
           # sort=c(1:ncol(vsd_comb2)), # overrides sorting of columns according to hierarchical clustering
           sort=colnames(vsd_comb2),
           cex=0.8,
           pdf=F,
)
dev.off()


#### BOXPLOTS MATCHING EXPERIMENTS ####

library(stringr)
library(rcompanion)
library(rstatix)
library(ggpubr)
library(scales)

# Mcavernosa9810 (Extracellular matrix binding)
Mcavernosa9810 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa9810.csv")
Mcavernosa9810$fate <- factor(Mcavernosa9810$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa9810_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa9810.csv")
Mcavernosa9810_trans$fate <- factor(Mcavernosa9810_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa9810_stats <- aov(count~fate,data=Mcavernosa9810) %>%
  tukey_hsd()
Mcavernosa9810_stats

Mcavernosa9810_trans_stats <- aov(count~fate,data=Mcavernosa9810_trans) %>%
  tukey_hsd()
Mcavernosa9810_trans_stats

Mcavernosa9810_plot <-
  ggboxplot(
    Mcavernosa9810,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa9810_stats,label="p.adj.signif",y.position=3,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa9810_plot

Mcavernosa9810_trans_plot <-
  ggboxplot(
    Mcavernosa9810_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa9810_trans_stats,label="p.adj.signif",y.position=3.1,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa9810_trans_plot

# Mcavernosa43816 (Spondin 2a, extracellular matrix protein)
Mcavernosa43816 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa43816.csv")
Mcavernosa43816$fate <- factor(Mcavernosa43816$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa43816_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa43816.csv")
Mcavernosa43816_trans$fate <- factor(Mcavernosa43816_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa43816_stats <- aov(count~fate,data=Mcavernosa43816) %>%
  tukey_hsd()
Mcavernosa43816_stats

Mcavernosa43816_trans_stats <- aov(count~fate,data=Mcavernosa43816_trans) %>%
  tukey_hsd()
Mcavernosa43816_trans_stats

Mcavernosa43816_plot <-
  ggboxplot(
    Mcavernosa43816,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa43816_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa43816_plot

Mcavernosa43816_trans_plot <-
  ggboxplot(
    Mcavernosa43816_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa43816_trans_stats,label="p.adj.signif",y.position=3,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa43816_trans_plot

# Mcavernosa14879 (Spondin 2b, extracellular matrix protein)
Mcavernosa14879 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa14879.csv")
Mcavernosa14879$fate <- factor(Mcavernosa14879$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa14879_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa14879.csv")
Mcavernosa14879_trans$fate <- factor(Mcavernosa14879_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa14879_stats <- aov(count~fate,data=Mcavernosa14879) %>%
  tukey_hsd()
Mcavernosa14879_stats

Mcavernosa14879_trans_stats <- aov(count~fate,data=Mcavernosa14879_trans) %>%
  tukey_hsd()
Mcavernosa14879_trans_stats

Mcavernosa14879_plot <-
  ggboxplot(
    Mcavernosa14879,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa14879_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa14879_plot

Mcavernosa14879_trans_plot <-
  ggboxplot(
    Mcavernosa14879_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa14879_trans_stats,label="p.adj.signif",y.position=3.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa14879_trans_plot

# Mcavernosa47647 (Activation of NF-kappaB-inducing kinase activity)
Mcavernosa47647 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa47647.csv")
Mcavernosa47647$fate <- factor(Mcavernosa47647$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa47647_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa47647.csv")
Mcavernosa47647_trans$fate <- factor(Mcavernosa47647_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa47647_stats <- aov(count~fate,data=Mcavernosa47647) %>%
  tukey_hsd()
Mcavernosa47647_stats

Mcavernosa47647_trans_stats <- aov(count~fate,data=Mcavernosa47647_trans) %>%
  tukey_hsd()
Mcavernosa47647_trans_stats

Mcavernosa47647_plot <-
  ggboxplot(
    Mcavernosa47647,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa47647_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa47647_plot

Mcavernosa47647_trans_plot <-
  ggboxplot(
    Mcavernosa47647_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa47647_trans_stats,label="p.adj.signif",y.position=2.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa47647_trans_plot

# Mcavernosa50735 (Activation of NF-kappaB-inducing kinase activity)
Mcavernosa50735 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa50735.csv")
Mcavernosa50735$fate <- factor(Mcavernosa50735$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa50735_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa50735.csv")
Mcavernosa50735_trans$fate <- factor(Mcavernosa50735_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa50735_stats <- aov(count~fate,data=Mcavernosa50735) %>%
  tukey_hsd()
Mcavernosa50735_stats

Mcavernosa50735_trans_stats <- aov(count~fate,data=Mcavernosa50735_trans) %>%
  tukey_hsd()
Mcavernosa50735_trans_stats

Mcavernosa50735_plot <-
  ggboxplot(
    Mcavernosa50735,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa50735_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa50735_plot

Mcavernosa50735_trans_plot <-
  ggboxplot(
    Mcavernosa50735_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa50735_trans_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa50735_trans_plot


# Mcavernosa61972 (Animal haem peroxidase)
Mcavernosa61972 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa61972.csv")
Mcavernosa61972$fate <- factor(Mcavernosa61972$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa61972_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa61972.csv")
Mcavernosa61972_trans$fate <- factor(Mcavernosa61972_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa61972_stats <- aov(count~fate,data=Mcavernosa61972) %>%
  tukey_hsd()
Mcavernosa61972_stats

Mcavernosa61972_trans_stats <- aov(count~fate,data=Mcavernosa61972_trans) %>%
  tukey_hsd()
Mcavernosa61972_trans_stats

Mcavernosa61972_plot <-
  ggboxplot(
    Mcavernosa61972,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.4,100), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa61972_stats,label="p.adj.signif",y.position=1.3,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa61972_plot

Mcavernosa61972_trans_plot <-
  ggboxplot(
    Mcavernosa61972_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limit = c(0.4,100), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa61972_trans_stats,label="p.adj.signif",y.position=1.2,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa61972_trans_plot

# Mcavernosa184695 (Non-membrane spanning protein tyrosine kinase activity)
Mcavernosa184695 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa184695.csv")
Mcavernosa184695$fate <- factor(Mcavernosa184695$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa184695_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa184695.csv")
Mcavernosa184695_trans$fate <- factor(Mcavernosa184695_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa184695_stats <- aov(count~fate,data=Mcavernosa184695) %>%
  tukey_hsd()
Mcavernosa184695_stats

Mcavernosa184695_trans_stats <- aov(count~fate,data=Mcavernosa184695_trans) %>%
  tukey_hsd()
Mcavernosa184695_trans_stats

Mcavernosa184695_plot <-
  ggboxplot(
    Mcavernosa184695,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa184695_stats,label="p.adj.signif",y.position=2,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa184695_plot

Mcavernosa184695_trans_plot <-
  ggboxplot(
    Mcavernosa184695_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa184695_trans_stats,label="p.adj.signif",y.position=1.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa184695_trans_plot


# Mcavernosa29964 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa29964 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa29964.csv")
Mcavernosa29964$fate <- factor(Mcavernosa29964$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa29964_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa29964.csv")
Mcavernosa29964_trans$fate <- factor(Mcavernosa29964_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa29964_stats <- aov(count~fate,data=Mcavernosa29964) %>%
  tukey_hsd()
Mcavernosa29964_stats

Mcavernosa29964_trans_stats <- aov(count~fate,data=Mcavernosa29964_trans) %>%
  tukey_hsd()
Mcavernosa29964_trans_stats

Mcavernosa29964_plot <-
  ggboxplot(
    Mcavernosa29964,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa29964_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa29964_plot

Mcavernosa29964_trans_plot <-
  ggboxplot(
    Mcavernosa29964_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa29964_trans_stats,label="p.adj.signif",y.position=3,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa29964_trans_plot

# Mcavernosa20827 (Positive regulation of protein tyrosine kinase activity)
Mcavernosa20827 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa20827.csv")
Mcavernosa20827$fate <- factor(Mcavernosa20827$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa20827_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa20827.csv")
Mcavernosa20827_trans$fate <- factor(Mcavernosa20827_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa20827_stats <- aov(count~fate,data=Mcavernosa20827) %>%
  tukey_hsd()
Mcavernosa20827_stats

Mcavernosa20827_trans_stats <- aov(count~fate,data=Mcavernosa20827_trans) %>%
  tukey_hsd()
Mcavernosa20827_trans_stats

Mcavernosa20827_plot <-
  ggboxplot(
    Mcavernosa20827,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa20827_stats,label="p.adj.signif",y.position=2.5,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa20827_plot

Mcavernosa20827_trans_plot <-
  ggboxplot(
    Mcavernosa20827_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.1,12000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa20827_trans_stats,label="p.adj.signif",y.position=2.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa20827_trans_plot


# Mcavernosa126229 (Transforming growth factor-beta (TGF-beta) family)
Mcavernosa126229 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa126229.csv")
Mcavernosa126229$fate <- factor(Mcavernosa126229$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa126229_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa126229.csv")
Mcavernosa126229_trans$fate <- factor(Mcavernosa126229_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa126229_stats <- aov(count~fate,data=Mcavernosa126229) %>%
  tukey_hsd()
Mcavernosa126229_stats

Mcavernosa126229_trans_stats <- aov(count~fate,data=Mcavernosa126229_trans) %>%
  tukey_hsd()
Mcavernosa126229_trans_stats

Mcavernosa126229_plot <-
  ggboxplot(
    Mcavernosa126229,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank()) +
  scale_y_log10(limit = c(0.25,1000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa126229_stats,label="p.adj.signif",y.position=1.75,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa126229_plot

Mcavernosa126229_trans_plot <-
  ggboxplot(
    Mcavernosa126229_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.x=element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank()) +
  scale_y_log10(limit = c(0.25,1000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa126229_trans_stats,label="p.adj.signif",y.position=2.25,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa126229_trans_plot

# Mcavernosa96261 (Transforming growth factor-beta (TGF-beta) family)
Mcavernosa96261 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa96261.csv")
Mcavernosa96261$fate <- factor(Mcavernosa96261$fate, levels = c("healthy", "diseased", "treated"))

Mcavernosa96261_trans <- read.csv(file = "../../transmission/DESeq2/mcav2015/Mcavernosa96261.csv")
Mcavernosa96261_trans$fate <- factor(Mcavernosa96261_trans$fate, levels = c("healthy", "nai", "diseased"))

# ANOVA and Tukey's
Mcavernosa96261_stats <- aov(count~fate,data=Mcavernosa96261) %>%
  tukey_hsd()
Mcavernosa96261_stats

Mcavernosa96261_trans_stats <- aov(count~fate,data=Mcavernosa96261_trans) %>%
  tukey_hsd()
Mcavernosa96261_trans_stats

Mcavernosa96261_plot <-
  ggboxplot(
    Mcavernosa96261,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    title = "Intervention",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") +
  scale_y_log10(limit = c(1,500), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa96261_stats,label="p.adj.signif",y.position=1.65,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa96261_plot

Mcavernosa96261_trans_plot <-
  ggboxplot(
    Mcavernosa96261_trans,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "orange", "red"),
    title = "Transmission",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_log10(limit = c(1,500), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa96261_trans_stats,label="p.adj.signif",y.position=2,step.increase=0.2,inherit.aes=FALSE,size=3)
Mcavernosa96261_trans_plot


#### MULTIPLOT MATCHING ####

intervention_orthologs<-ggarrange(Mcavernosa61972_plot,
                           Mcavernosa61972_trans_plot,
                           Mcavernosa126229_plot,
                           Mcavernosa126229_trans_plot,
                           Mcavernosa96261_plot,
                           Mcavernosa96261_trans_plot,
                           heights = c(4,4,4),
                           widths = c(3.1,3),
                           ncol = 2,
                           nrow = 3)
intervention_orthologs<-annotate_figure(intervention_orthologs, top = text_grob("Animal haem peroxidase", color = "black", face = "bold", size = 14), 
                                 bottom = text_grob("Transforming growth factor-beta (TGF-beta) family", color = "black", face = "bold", size = 14))
intervention_orthologs

ggsave("intervention_orthologs.pdf", intervention_orthologs, width=8, height=10,dpi = 300)


#### BOXPLOTS INTERVENTION ####

# Mcavernosa12872 (Fibulin-like extracellular matrix protein 2)
Mcavernosa12872 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa12872.csv")
Mcavernosa12872$fate <- factor(Mcavernosa12872$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa12872$title <- "Fibulin-like extracellular matrix protein 2"

# ANOVA and Tukey's
Mcavernosa12872_stats <- aov(count~fate,data=Mcavernosa12872) %>%
  tukey_hsd()
Mcavernosa12872_stats

Mcavernosa12872_plot <-
  ggboxplot(
    Mcavernosa12872,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="salmon"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa12872_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa12872_plot

# Mcavernosa10679 (Negative regulation of I-kappaB kinase/NF-kappaB signaling)
Mcavernosa10679 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa10679.csv")
Mcavernosa10679$fate <- factor(Mcavernosa10679$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa10679$title <- "Negative regulation of I-kappaB kinase/NF-kappaB signaling"

# ANOVA and Tukey's
Mcavernosa10679_stats <- aov(count~fate,data=Mcavernosa10679) %>%
  tukey_hsd()
Mcavernosa10679_stats

Mcavernosa10679_plot <-
  ggboxplot(
    Mcavernosa10679,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa10679_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa10679_plot

# Mcavernosa49779 (Glutathione peroxidase activity)
Mcavernosa49779 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa49779.csv")
Mcavernosa49779$fate <- factor(Mcavernosa49779$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa49779$title <- "Glutathione peroxidase activity"

# ANOVA and Tukey's
Mcavernosa49779_stats <- aov(count~fate,data=Mcavernosa49779) %>%
  tukey_hsd()
Mcavernosa49779_stats

Mcavernosa49779_plot <-
  ggboxplot(
    Mcavernosa49779,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa49779_stats,label="p.adj.signif",y.position=2.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa49779_plot

# Mcavernosa5069 (Glutathione peroxidase activity)
Mcavernosa5069 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa5069.csv")
Mcavernosa5069$fate <- factor(Mcavernosa5069$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa5069$title <- "Glutathione peroxidase activity"

# ANOVA and Tukey's
Mcavernosa5069_stats <- aov(count~fate,data=Mcavernosa5069) %>%
  tukey_hsd()
Mcavernosa5069_stats

Mcavernosa5069_plot <-
  ggboxplot(
    Mcavernosa5069,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa5069_stats,label="p.adj.signif",y.position=1.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa5069_plot

# Mcavernosa102943 (Transmembrane receptor protein tyrosine kinase activity)
Mcavernosa102943 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa102943.csv")
Mcavernosa102943$fate <- factor(Mcavernosa102943$fate, levels = c("healthy", "diseased", "treated"))
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
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.text.x=element_blank(), axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa102943_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa102943_plot

# Mcavernosa71973 (Non-membrane spanning protein tyrosine kinase activity)
Mcavernosa71973 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa71973.csv")
Mcavernosa71973$fate <- factor(Mcavernosa71973$fate, levels = c("healthy", "diseased", "treated"))
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
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", strip.background = element_rect(fill="plum"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa71973_stats,label="p.adj.signif",y.position=1.75,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa71973_plot

# Mcavernosa162020 (WD repeat-containing protein 37)
Mcavernosa162020 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa162020.csv")
Mcavernosa162020$fate <- factor(Mcavernosa162020$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa162020$title <- "WD repeat-containing protein 37"

# ANOVA and Tukey's
Mcavernosa162020_stats <- aov(count~fate,data=Mcavernosa162020) %>%
  tukey_hsd()
Mcavernosa162020_stats

Mcavernosa162020_plot <-
  ggboxplot(
    Mcavernosa162020,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
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
        axis.text.y = element_text(hjust=0)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa162020_stats,label="p.adj.signif",y.position=0.85,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa162020_plot

# Mcavernosa317111 (F-box WD repeat-containing protein)
Mcavernosa317111 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa317111.csv")
Mcavernosa317111$fate <- factor(Mcavernosa317111$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa317111$title <- "F-box WD repeat-containing protein"

# ANOVA and Tukey's
Mcavernosa317111_stats <- aov(count~fate,data=Mcavernosa317111) %>%
  tukey_hsd()
Mcavernosa317111_stats

Mcavernosa317111_plot <-
  ggboxplot(
    Mcavernosa317111,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
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
        axis.title.y=element_blank(), axis.text.y = element_blank()) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa317111_stats,label="p.adj.signif",y.position=0.95,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa317111_plot

# Mcavernosa43057 (WD repeat-containing protein 41)
Mcavernosa43057 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa43057.csv")
Mcavernosa43057$fate <- factor(Mcavernosa43057$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa43057$title <- "WD repeat-containing protein 41"

# ANOVA and Tukey's
Mcavernosa43057_stats <- aov(count~fate,data=Mcavernosa43057) %>%
  tukey_hsd()
Mcavernosa43057_stats

Mcavernosa43057_plot <-
  ggboxplot(
    Mcavernosa43057,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
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
  stat_pvalue_manual(Mcavernosa43057_stats,label="p.adj.signif",y.position=2.25,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa43057_plot

# Mcavernosa10151 (WD repeat-containing protein 82)
Mcavernosa10151 <- read.csv(file = "../DESeq2/mcav2015/Mcavernosa10151.csv")
Mcavernosa10151$fate <- factor(Mcavernosa10151$fate, levels = c("healthy", "diseased", "treated"))
Mcavernosa10151$title <- "WD repeat-containing protein 82"

# ANOVA and Tukey's
Mcavernosa10151_stats <- aov(count~fate,data=Mcavernosa10151) %>%
  tukey_hsd()
Mcavernosa10151_stats

Mcavernosa10151_plot <-
  ggboxplot(
    Mcavernosa10151,
    x = "fate",
    y = "count",
    color = "grey30",
    fill = "fate",
    palette = c("green", "red", "orange"),
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'fate') + 
  theme_classic() +    facet_grid(. ~ title) +
  facet_grid(. ~ title) +
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right", strip.background = element_rect(fill="lightskyblue"), strip.text = element_text(size=9), 
        axis.title.y=element_blank(), axis.text.y = element_blank()) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                 labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Mcavernosa10151_stats,label="p.adj.signif",y.position=2.2,step.increase=0.1,inherit.aes=FALSE,size=3)
Mcavernosa10151_plot

# creating legend for large plot
legend_mcav <- get_legend(Mcavernosa10151_plot)

Mcavernosa10151_plot <- Mcavernosa10151_plot + theme(legend.position = "none")


#### MULTIPLOT INTERVENTION ####

intervention_mcav <-ggarrange(Mcavernosa102943_plot,
                          Mcavernosa71973_plot,
                          legend_mcav,
                          Mcavernosa162020_plot,
                          Mcavernosa317111_plot,
                          Mcavernosa10151_plot,
                          heights = c(4,4),
                          widths = c(5.1,5,5),
                          ncol = 3,
                          nrow = 2)
intervention_mcav
ggsave("intervention_mcav.pdf", intervention_mcav, width=10, height=8,dpi = 300)
