# Make a file of allele proportion across colonies so can plot across SNPS
# for Ben's graph suggestion

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/counts")

library(readr)
library(sqldf)
require(reshape2)
library(doBy)
library(dplyr)
library(tidyr)
library(ggplot2)


file.list = list.files(("./"),pattern="*.csv")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE)
}

samples <- lapply(file.list, read_file1)


# cut down dataframes and filter on total coverage of >10 and add columns where the more expessed
# SNP is in one column and the lower expressed in the other
# Also remove SNPs which have 0 count in the ref/alt allele, could be homozygous and miscalled

num_snps_start <- vector("numeric", 18)
num_snps_10coverage <- vector("numeric", 18)
num_snps_final <- vector("numeric", 18)

for (i in seq_along(samples)){
  num_snps_start[i] <- nrow(samples[[i]])
  samples[[i]] <- subset(samples[[i]], totalCount > 10)
  
  num_snps_10coverage[i] <- nrow(samples[[i]])
  samples[[i]] <- samples[[i]][,c(1,2,6,7)]
  
  samples[[i]]$allele1 <- ifelse(samples[[i]]$refCount < samples[[i]]$altCount
                                 , samples[[i]]$altCount, samples[[i]]$refCount)
  
  samples[[i]]$allele2 <- ifelse(samples[[i]]$refCount > samples[[i]]$altCount
                                 , samples[[i]]$altCount, samples[[i]]$refCount)
  
  samples[[i]] <- subset(samples[[i]], refCount > 0 & altCount > 0)
  num_snps_final[i] <- nrow(samples[[i]])
}


sample_names <- gsub("_allele_specific_counts_per_SNP.csv", "", file.list)
sample_names <- gsub("trim_", "", sample_names)
names(samples) <- sample_names

genes<- read_table2("genes_with_total_cpgs.txt")

samples2 <- list()

for(i in seq_along(samples)){
  df <- samples[[i]]
  output <- sqldf("SELECT sg.contig,
                  sg.position,
                  sg.allele1,
                  sg.allele2,
                  fg.geneID,
                  fg.start,
                  fg.end,
                  fg.chr
                  FROM df AS sg
                  LEFT JOIN genes AS fg 
                  ON sg.contig = fg.chr
                  AND (sg.position >= fg.start AND sg.position <= fg.end)")
  output <- output[!is.na(output$geneID),]
  samples2[[i]] <- output[,-c(6:8)]
}

names(samples2) <- sample_names

# Remove samples not tested for ASE (j818 and j824)
names(samples2)
samples2[[14]] <- NULL
names(samples2)
samples2[[15]] <- NULL
names(samples2)

repro <- samples2[c(2,4,5,7,9,11,13,16)]
names(repro)
  
sterile <- samples2[c(1,3,6,8,10,12,14,15)]
names(sterile)

sterile_all<- bind_rows(sterile, .id = "sample")
repro_all<- bind_rows(repro, .id = "sample")

sterile_combined <- summaryBy(allele1 + allele2 ~ geneID + position, data=sterile_all, FUN=sum) 
repro_combined <- summaryBy(allele1 + allele2 ~ geneID + position, data=repro_all, FUN=sum) 

sterile_combined$total <- sterile_combined$allele1.sum + sterile_combined$allele2.sum
sterile_combined$proportion <- sterile_combined$allele1.sum / sterile_combined$total
sterile_combined$status <- "sterile"

repro_combined$total <- repro_combined$allele1.sum + repro_combined$allele2.sum
repro_combined$proportion <- repro_combined$allele1.sum / repro_combined$total
repro_combined$status <- "repro"

all_data <- rbind(repro_combined, sterile_combined)

all_final <- summaryBy(proportion ~ geneID + position, data=all_data, FUN=mean) 

sig_ASE_repro_and_sterile <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/sig_ASE_repro_and_sterile.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
all_final_only_ASE <- all_final[(all_final$geneID %in% sig_ASE_repro_and_sterile$geneID), ]
#all_final_no_ASE <- all_final[!(all_final$geneID %in% sig_ASE_repro_and_sterile$geneID), ]


# Make a plot to check the allele count / ratio, (from Ben)

# Random non-ASE genes
test <- subset(all_final, geneID == "LOC100643962") # 10 snps
test1 <- subset(all_final, geneID == "LOC105667086") # 10 snps
test2 <- subset(all_final, geneID == "LOC100649484") # 3 snps
test3 <- subset(all_final, geneID == "LOC105666889") # 38 snps

# ASE genes
ase1 <- subset(all_final, geneID == "LOC100650626") 
ase2 <- subset(all_final, geneID == "LOC105666882") 
ase3 <- subset(all_final, geneID == "LOC100646400") 
ASE_genes <- as.data.frame(unique((all_final_only_ASE$geneID))) #98

genes <- ASE_genes$`unique((all_final_only_ASE$geneID))`

pdf("testing.pdf")
for (i in genes){
  ase1 <- subset(all_final, geneID == i)
 # print(i)
  print(ggplot(ase1, aes(x=position, group=geneID))+
    geom_line(aes(y=proportion.mean))+ylim(0.5,1.0))
}
dev.off()

# Graphs look crap 
# Check mean number of SNPs per gene anyway ...
snp_count <- as.data.frame(table(all_final$geneID))
mean(snp_count$Freq) # 7.5 SNPs
sd(snp_count$Freq) # 12.4
median(snp_count$Freq) # 4
range(snp_count$Freq) # 1-283

snp_count_ASE <- as.data.frame(table(all_final$geneID[all_final$geneID
                                                      %in% ASE_genes$`unique((all_final_only_ASE$geneID))`]))
mean(snp_count_ASE$Freq) # 21 SNPs
sd(snp_count_ASE$Freq) # 18
median(snp_count_ASE$Freq) # 16
range(snp_count_ASE$Freq) # 1-88
