# What does the expression of genes with ASM look like?

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASM")

library(readr)
library(rowr)
library(ggplot2)

nonrepro_amr_genes_unique_list <- read_csv("Gene_Lists/nonrepro_amr_genes_unique_list.csv", 
                                           col_names = FALSE)
colnames(nonrepro_amr_genes_unique_list) <- "geneID"

repro_amr_genes_unique_list <- read_csv("Gene_Lists/repro_amr_genes_unique_list.csv", 
                                           col_names = FALSE)
colnames(repro_amr_genes_unique_list) <- "geneID"

both_amr_genes_unique_list <- read_csv("Gene_Lists/both_amr_genes.csv", 
                                        col_names = FALSE)
colnames(both_amr_genes_unique_list) <- "geneID"


# Overlap with any ASE?

ase <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/sig_ASE_repro_and_sterile.txt", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)

ase_repro <- merge(repro_amr_genes_unique_list, ase, by="geneID") # 0
ase_srerile <- merge(nonrepro_amr_genes_unique_list, ase, by="geneID") #1  LOC105666711 tyrosine-protein kinase Btk29A (immunity)

ase_both <- merge(both_amr_genes_unique_list, ase, by="geneID") #6

# General trend for genes showing ASM to be more/less ASE?
ggplot(ase_both, aes(x=status, y=propmatexpr, fill=status))+
  geom_boxplot()+
  geom_point()

## GO BACK TO ASE SCRIPT
