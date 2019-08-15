## Organise gene lists:

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Allele_specific_MERN/Methylation/Gene_Lists")

library(readr)

both <- read_csv("both_amr_genes.csv", 
                           col_names = FALSE) #97

nonrepro <- read_csv("nonrepro_amr_genes_full_list.csv", 
                                         col_names = FALSE) #195

repro <- read_csv("repro_amr_genes_full_list.csv", 
                                      col_names = FALSE) #287

# Get list of genes only found in repro / nonrepro 
repro_unique <- repro[!(repro$X1 %in% both$X1),]
write.csv(repro_unique, file="repro_amr_genes_unique_list.csv") #190


nonrepro_unique <- nonrepro[!(nonrepro$X1 %in% both$X1),]
write.csv(nonrepro_unique, file="nonrepro_amr_genes_unique_list.csv") #98
