#########################################################
# Obtain gene names from NCBI given a LOC id
#########################################################

# Make a text file with each gene id on a separate line
# Go to: https://www.ncbi.nlm.nih.gov/sites/batchentrez
# Select the gene database, load the list and search
# You can then export the output and use script below to tag the names on to your original file

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Allele_specific_MERN/Methylation")

library(readr)


original_file<-read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/Allele_specific_MERN/Methylation/amr_annotated_with_gene_and_exon.csv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

batch_results<-read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/Allele_specific_MERN/Methylation/gene_names_from_NCBI.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(batch_results)[6]<-"gene"

output<-merge(original_file, batch_results, by="gene")
output<-output[,c(1:8,15)]

write.csv(output, file="amr_annotated_with_gene_exon_name.csv")
