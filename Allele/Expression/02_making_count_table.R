## Making counts for genes rather than SNPs for ASE

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

snp_counts <- melt(data.frame(num_snps_start,num_snps_10coverage,num_snps_final,sample_names))
colnames(snp_counts)<- c("sample","snp_condition","count")
snp_counts2 <- dcast(snp_counts, sample ~ snp_condition)

# SNP count bar plot
ggplot(data=snp_counts2, aes(x=sample))+
  geom_bar(aes(y=num_snps_start), stat="identity",fill="dodgerblue",show.legend=T)+
  geom_bar(aes(y=num_snps_10coverage), stat="identity",fill="skyblue2",show.legend=T)+
  geom_bar(aes(y=num_snps_final), stat="identity",fill="turquoise",show.legend=T)+
  ylab("SNP Count")+
  geom_text(data = snp_counts2, aes(label = num_snps_start, y =num_snps_start -500),
            size=4)+
  #geom_text(data=snp_counts2, aes(label=num_snps_10coverage, y=num_snps_10coverage -400),
   #         size=4)+
  geom_text(data = snp_counts2, aes(label = num_snps_final, y = num_snps_final - 500),
            size=4)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20))+
  xlab("Sample")

mean(snp_counts2$num_snps_start) # 17753
sd(snp_counts2$num_snps_start) # 6840
mean(snp_counts2$num_snps_10coverage) # 9355
sd(snp_counts2$num_snps_10coverage) # 3781
mean(snp_counts2$num_snps_final) # 9297
sd(snp_counts2$num_snps_final) # 3755

# annotate SNPs with gene ID
genes<- read_table2("genes_with_total_cpgs.txt")

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
  output <- output[,-c(6:8)]
  
  output2 <- summaryBy(allele1 + allele2 ~ geneID, data=output, FUN=sum) 
  
  output3 <- melt(output2)
  output3$geneID <- paste(output3$geneID,"_",output3$variable)
  output3$geneID <- gsub(" ","",output3$geneID)
  output3$geneID <- gsub(".sum","",output3$geneID)
  output3 <- output3[,-2]
  colnames(output3) <- NULL
  
  myfile <- file.path("../counts_per_gene/", paste(names(samples)[i],"_","allele_counts_per_gene.txt"))
  write.table(output3, file=myfile, quote=F, sep="\t", row.names=F)
}



# Make a count matrix now
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/counts_per_gene")

file.list = list.files(("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/counts_per_gene"),pattern="*.txt")
samples = as.character(sapply(file.list,function(s) gsub(" _ allele_counts_per_gene.txt","", s)))
samples

info = read.csv("../sample_information.csv")

out = do.call(rbind, lapply(1:length(samples),
                            function(i) {
                              cnts = read_delim(file.list[[i]], "\t", escape_double = F, col_names = F, trim_ws = T) # counts
                              colnames(cnts) = c("geneID","count")
                              sampinfo = info[info$sample==samples[[i]],] # sample info
                              suppressWarnings(data.frame(sampinfo, cnts))
                            } ))

out <- out %>% separate(geneID, c("geneID", "allele"), "_a")
out$sample <- paste(out$sample, "_a", out$allele)
out$sample <- gsub(" " , "", out$sample)
out <- out[,-5]

length(out[out$geneID=="LOC100631058",]$sample) # Test: there are 36 lvls good
write.csv(out, "count_matrix_ASE.csv", row.names = F)

# Count number of genes per sample
look <- out %>% 
  group_by(sample) %>%
  summarise(no_rows = length(sample))

look<- look[grep("allele1", look$sample), ]
look$sample <- gsub("_allele1", "", look$sample)


ggplot(data=look, aes(x=sample))+
  geom_bar(aes(y=no_rows), stat="identity",fill="dodgerblue",show.legend=F)+
 ylab("Gene Count")+
  geom_text(data = look, aes(label = no_rows, y =no_rows -100),
            size=4)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20))+
  xlab("Sample")

mean(look$no_rows)
sd(look$no_rows)


# Make a different count matrix now
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/counts_per_gene")

file.list = list.files(("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/counts_per_gene"),pattern="*.txt")
samples = as.character(sapply(file.list,function(s) gsub(" _ allele_counts_per_gene.txt","", s)))
samples

info = read.csv("../sample_information.csv")

out = do.call(rbind, lapply(1:length(samples),
                            function(i) {
                              cnts = read_delim(file.list[[i]], "\t", escape_double = F, col_names = F, trim_ws = T) # counts
                              colnames(cnts) = c("geneID","count")
                              sampinfo = info[info$sample==samples[[i]],] # sample info
                              suppressWarnings(data.frame(sampinfo, cnts))
                            } ))

out <- out %>% separate(geneID, c("geneID", "allele"), "_a")
out$allele <- gsub("ll","all", out$allele)
out <- subset(out, !sample == "j8_24") 
out <- out[,-1]

out1 <- subset(out, allele == "allele1")
colnames(out1)[5] <- "allele1"
out1 <- out1[,-4]

out2 <- subset(out, allele == "allele2")
colnames(out2)[5] <- "allele2"
out2 <- out2[,-4]

final <- cbind(out1, out2$allele2)
colnames(final)[5] <- "allele2"

final1<- summaryBy(allele1 + allele2 ~ geneID + colony + status, data=final, FUN=sum) 
final1$total_counts <- final1$allele1.sum + final1$allele2.sum
final1$prop.allele1 <- final1$allele1.sum / final1$total_counts
final1$prop.allele2 <- final1$allele2.sum / final1$total_counts

write.table(final1, "count_matrix_ASE_for_model.csv", row.names = F,
          col.names = T, quote = F, sep = '\t')






