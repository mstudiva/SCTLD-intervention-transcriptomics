#### PACKAGES ####

library(tidyverse)
library(VennDiagram)

#### orthofinder ####

# Install orthofinder on your local machine using the tutorials (https://davidemms.github.io/menu/tutorials.html)

# Copy your translated protein fasta files (_out_PRO.fas) that you want to compare into a directory called 'orthofinder'
# If you have not already filtered by the longest contig per isogroup (by using fasta2SBH.pl during transcriptome annotation), follow step 7 of tutorial 2 above

# Run the following command in Terminal: 'orthofinder -f orthofinder/'
# Check the number of genes assigned to orthogroups (e.g., 'OrthoFinder assigned 29028 genes (87.6% of total) to 8731 orthogroups')
# Ideally, it should be >80%


#### ORTHOLOGS ####

orthologs <- read.table(file = "OrthoFinder/Results_May10/Orthologues/Orthologues_Ofaveolata_out_PRO/Ofaveolata_out_PRO__v__Mcavernosa2015_out_PRO.tsv", sep = "\t", header = TRUE)


#### DESEQ IMPORT ####

ofav_diseased_healthy_lpv <- read.csv(file = "../transmission/DESeq2/ofav/diseased_healthy_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

mcav_diseased_healthy_lpv <- read.csv(file = "../transmission/DESeq2/mcav2015/diseased_healthy_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_mcav" = lpv)

mcav_diseased_nai_lpv <- read.csv(file = "../transmission/DESeq2/mcav2015/diseased_nai_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_mcav" = lpv)

mcav_nai_healthy_lpv <- read.csv(file = "../transmission/DESeq2/mcav2015/nai_healthy_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_mcav" = lpv)


#### DEG MATCH ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.05 pval cutoff, 4) then add ofav and mcav gene annotations

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
  filter(abs(lpv_mcav) >= 1.30103 & abs(lpv_ofav) >= 1.30103) %>%
  left_join(read.table(file = "../annotate/ofav/Ofaveolata_iso2geneName.tab",
             sep = "\t") %>%
               mutate(gene = V1,
                      annot_ofav = V2) %>%
               dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t") %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> diseased_healthy
             
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
  filter(abs(lpv_mcav) >= 1.30103 & abs(lpv_ofav) >= 1.30103) %>%
  left_join(read.table(file = "../annotate/ofav/Ofaveolata_iso2geneName.tab",
                       sep = "\t") %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t") %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> diseased_nai

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
  filter(abs(lpv_mcav) >= 1.30103 & abs(lpv_ofav) >= 1.30103) %>%
  left_join(read.table(file = "../annotate/ofav/Ofaveolata_iso2geneName.tab",
                       sep = "\t") %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../annotate/mcav2015/Mcavernosa2015_iso2geneName.tab",
                       sep = "\t") %>%
              mutate(gene = V1,
                     annot_mcav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_mcav" = "gene")) -> nai_healthy


#### VENN DIAGRAMS ####

# first creating a set of up/downregulated DEGs by species
diseased_healthy %>%
  filter(lpv_mcav >= 1.30103) %>%
  pull(Orthogroup) -> mcav_up

diseased_healthy %>%
  filter(lpv_mcav <= -1.30103) %>%
  pull(Orthogroup) -> mcav_down

diseased_healthy %>%
  filter(lpv_ofav >= 1.30103) %>%
  pull(Orthogroup) -> ofav_up

diseased_healthy %>%
  filter(lpv_ofav <= -1.30103) %>%
  pull(Orthogroup) -> ofav_down


list=list("Mcav up"=mcav_up, "Mcav down"=mcav_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down)

# treatment/time contrasts
venn=venn.diagram(
  x = list,
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
