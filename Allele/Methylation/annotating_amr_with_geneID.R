### Annotating allelically methylated regions with gene ID and exon if applicable ######
setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASM/Intermediate_files")
library(plyr)
library(dplyr)
library(sqldf)
library(readr)
library(ggplot2)


# Read in annotation file, containing gene name and exon info
genome_annotation<-read.csv.sql("annotation_information_without_quotes.csv",
                           sql ="select * from file", sep=",",header = T)
genome_annotation<-genome_annotation[,c(2,6,3,4,7,9)]
colnames(genome_annotation)<-c("scaffold","gene","start","end","feature","exon_number")
genome_annotation<-subset(genome_annotation, c(feature=="gene"|feature=="exon"|feature=="intron"))


# Read in the .amr files
repro_amrs <- read_delim("repro_3CpG_amr.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
nonrepro_amrs <- read_delim("nonrepro_3CpG_amr.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


# Cut down file sizes to make merge quicker later 
repro_amrs<-repro_amrs[,c(1,2,3,5)]
nonrepro_amrs<-nonrepro_amrs[,c(1,2,3,5)]


# Add column names
colnames(repro_amrs)<-c("scaffold","amr_start","amr_end","pval")
colnames(nonrepro_amrs)<-c("scaffold","amr_start","amr_end","pval")


# Now annotate the files with gene information, for the the AMR must fall completely within a gene or exon to be annotated
output_repro <- sqldf("SELECT repro_amrs.scaffold,
      repro_amrs.amr_start,
      repro_amrs.amr_end,
      repro_amrs.pval,
      genome_annotation.scaffold,
      genome_annotation.gene,
      genome_annotation.start,
      genome_annotation.end,
      genome_annotation.feature,
      genome_annotation.exon_number
      FROM repro_amrs AS repro_amrs
      LEFT JOIN genome_annotation AS genome_annotation 
      ON repro_amrs.scaffold = genome_annotation.scaffold
      AND (repro_amrs.amr_start >= genome_annotation.start AND repro_amrs.amr_end <= genome_annotation.end)") 

output_nonrepro <- sqldf("SELECT nonrepro_amrs.scaffold,
      nonrepro_amrs.amr_start,
      nonrepro_amrs.amr_end,
      nonrepro_amrs.pval,
      genome_annotation.scaffold,
      genome_annotation.gene,
      genome_annotation.start,
      genome_annotation.end,
      genome_annotation.feature,
      genome_annotation.exon_number
      FROM nonrepro_amrs AS nonrepro_amrs
      LEFT JOIN genome_annotation AS genome_annotation 
      ON nonrepro_amrs.scaffold = genome_annotation.scaffold
      AND (nonrepro_amrs.amr_start >= genome_annotation.start AND nonrepro_amrs.amr_end <= genome_annotation.end)")


# Remove duplicate rows
output_repro<-output_repro[!duplicated(output_repro),]
output_nonrepro<-output_nonrepro[!duplicated(output_nonrepro),]

# Remove unimportant information
output_repro<-output_repro[,-c(5,7,8,10)]
output_nonrepro<-output_nonrepro[,-c(5,7,8,10)]

write.csv(output_repro, file="repro_amr_genes_all_info.csv")
write.csv(output_nonrepro, file="nonrepro_amr_genes_all_info.csv")

# Count number not within a gene/exon/intron
repro_count <- output_repro[is.na(output_repro$gene),] #26 / 303
sterile_count <- output_nonrepro[is.na(output_nonrepro$gene),] #15 / 201


# Get gene list for AMRs in repro or non-repro
gene_list_repro<-as.data.frame(unique(output_repro$gene)) #287
gene_list_nonrepro<-as.data.frame(unique(output_nonrepro$gene)) #196

write.csv(gene_list_repro, file="repro_amr_genes.csv")
write.csv(gene_list_nonrepro, file="nonrepro_amr_genes.csv")


# AMR genes in common between repro and non repro (97)
colnames(gene_list_nonrepro)<-"gene"
colnames(gene_list_repro)<-"gene"
both_gene_list<-merge(gene_list_nonrepro, gene_list_repro, by="gene")
both_gene_list <- as.data.frame(both_gene_list[!is.na(both_gene_list$gene),])

write.csv(duplicated_genes, file="both_amr_genes.csv")

# 97 in common
# 190 unique genes in repro
# 98 unique in sterile

# Significant overlap? 11030 genes in genome.
phyper(97, 287, 10743, 195, lower.tail = F) # p = 2.701932e-107

library(UpSetR)

upset_repro <- as.data.frame(gene_list_repro[!is.na(gene_list_repro$gene),])
colnames(upset_repro) <- "geneID"
upset_sterile <- as.data.frame(gene_list_nonrepro[!is.na(gene_list_nonrepro$gene),])
colnames(upset_sterile) <- "geneID"

upset_repro$Reproductive <- 1
upset_sterile$Sterile <- 1

upset_data <- merge(upset_repro, upset_sterile, by="geneID", all=T)
upset_data[is.na(upset_data)] <- 0
upset_data <- upset_data[!duplicated(upset_data),]

upset(upset_data, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity")

# Make just one file for further analysis with repro_status as column
#output_repro$repro_status<-"repro"
#output_nonrepro$repro_status<-"non_repro"

#final_output<-rbind(output_nonrepro, output_repro)
#write.table(final_output, file="amr_annotated_with_gene_and_exon.csv", sep="\t", row.names=F, quote=F)


# Now lets look if they're in exons or introns within a gene (no neither)

exons_repro <- subset(output_repro, feature =="exon") 
exons_repro <- exons_repro[!duplicated(exons_repro),] #49
introns_repro <- subset(output_repro, feature =="intron") 
introns_repro <- introns_repro[!duplicated(introns_repro),] #7
both_repro <- merge(exons_repro, introns_repro, by= "gene") # 0

gene_repro <- subset(output_repro, feature =="gene")
gene_repro <- gene_repro[!duplicated(gene_repro),] 
gene_repro_not_exon <- gene_repro[(!(gene_repro$gene %in% exons_repro$gene)),]
or_intron <- gene_repro_not_exon[(!(gene_repro_not_exon$gene %in% introns_repro$gene)),]
look <- as.data.frame(unique(or_intron$gene)) # 231 (231+49+7 = 287)



exons_sterile <- subset(output_nonrepro, feature =="exon") 
exons_sterile <- exons_sterile[!duplicated(exons_sterile),] #35
introns_sterile <- subset(output_nonrepro, feature =="intron") 
introns_sterile <- introns_sterile[!duplicated(introns_sterile),] #6
both_sterile<- merge(exons_sterile, introns_sterile, by= "gene") # 0

gene_sterile <- subset(output_nonrepro, feature =="gene")
gene_sterile <- gene_sterile[!duplicated(gene_sterile),] 
gene_sterile_not_exon <- gene_sterile[(!(gene_sterile$gene %in% exons_sterile$gene)),]
or_intron_Sterile <- gene_sterile_not_exon[(!(gene_sterile_not_exon$gene %in% introns_sterile$gene)),]
look <- as.data.frame(unique(or_intron_Sterile$gene)) # 154 (154+35+6 = 196)



# Test of independence for enrichment in exons? and diff between castes?
R1 = c(49, 7)
R2 = c(35, 6)

rows   = 2

Matriz = matrix(c(R1, R2),
                nrow=rows,
                byrow=TRUE)

rownames(Matriz) = c("Repro", "Sterile")          # Naming the rows and
colnames(Matriz) = c("Exons", "Introns")   #  columns is optional.
chisq.test(Matriz) # X-squared = 9.672e-06, df = 1, p-value = 0.9975
# Meaning repro/sterile doesn't affect the ratio between exons/introns


# Are these the same exons and introns in repro and sterile?
No_Annotation <- c(231, 154)
Exon <- c(49, 35)
Intron <- c(7, 6)
dat <- data_frame(No_Annotation,Exon,Intron)
dat$Status <- c("Reproductive", "Sterile")


ggplot(data=dat, aes(x=Status))+
  geom_bar(aes(y=No_Annotation), stat="identity",fill="dodgerblue",show.legend=T)+
  geom_bar(aes(y=Exon), stat="identity",fill="skyblue2",show.legend=T)+
  geom_bar(aes(y=Intron), stat="identity",fill="turquoise",show.legend=T)+
  ylab("Allelically Methylated Region Count")+
  geom_text(data = dat, aes(label = No_Annotation, y= No_Annotation -3),
            size=6)+
  geom_text(data=dat, aes(label=Exon, y=Exon -3 ),
            size=6)+
  geom_text(data = dat, aes(label = Intron, y =Intron - 3),
            size=6)+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18))

# Extra plot to get the legend
library(reshape2)
dat1 <- melt(dat)
ggplot(data=dat1, aes(x=Status, y=value, fill=variable))+
  geom_bar(stat="identity", position="fill")+
  ylab("SNP Count")+
  ggtitle("Number of SNPs per Sample")+
  scale_fill_manual(breaks=c("No_Annotation","Exon","Intron"),
                    labels = c("No Annotation","Exon","Intron"),
                    values= c("dodgerblue","skyblue2","turquoise"))+
  theme(legend.text = element_text(size=18))




