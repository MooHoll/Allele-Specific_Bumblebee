###-------------------------------------------------------------------
# BINOMIAL GLMs TO DETERMINE ASE
###-------------------------------------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE")
library(pbapply)
library(plyr)
library(doBy)
library(export) 
library(effects)
library(devtools)
library(broom)
library(dplyr)
library(afex)
library(lsmeans)
library(car)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(readr)
library(UpSetR)
library(effects)
library (lme4)
library(ggpubr)
library(reshape2)

###------------------------------------------------------------------
# Read in count data
###------------------------------------------------------------------

df <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE/count_matrix_ASE_for_model.csv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)


###------------------------------------------------------------------
# Look at data
###------------------------------------------------------------------

# Look at the summed counts over all the SNP positions for each gene
 range(df$prop.allele1) # 0.5000000 0.9964093
 hist(df$prop.allele1,col="steelblue",
      xlab="Proportion of allele1 expression")
 
 range(df$prop.allele2) # 0.003590664 0.500000000
 hist(df$prop.allele2,col="steelblue",
      xlab="Proportion of allele1 expression")


###------------------------------------------------------------------
# Filter data to keep only genes in all crosses/conditions
###------------------------------------------------------------------


# See which genes occur in all experimental conditions/crosses
df <- as.data.frame(df)
df1 = summaryBy(allele1.sum+allele2.sum ~ geneID, data=df, FUN=length)


# Keep only genes for which we have counts in all 6 conditions (repro/sterile x3 colonies)
# This leave only 2489 genes
selectedgenes = df1$geneID[df1$allele1.sum.length == 6]
length(selectedgenes)
df2 = df[(df$geneID %in% selectedgenes),]

# Remove colony and pool by status (as missing two from colony 8?)
df3 = summaryBy(allele1.sum+allele2.sum ~ geneID + status + colony, data=df2, FUN=sum)
df3 <- df3[!is.na(df3),]

###------------------------------------------------------------------
# Run model
###------------------------------------------------------------------

set_sum_contrasts()
head(df3)
df3 <- as.data.frame(df3)
length(unique(df3$geneID)) # 2489
colnames(df3) <- c("geneID","status","colony","allele1","allele2")
df3 <- df3[!is.na(df3),]
df3$status <- as.factor(df3$status)
df3$colony <- as.factor(df3$colony)
df3$geneID <- as.factor(df3$geneID)

# Test with 1 gene
gene = "LOC100628631" 
dfsubs = subset(df3, geneID==gene)

# Take into account overdispersion using quasibinomial
fit = glm(cbind(allele1, allele2) ~ status + colony, 
           family = quasibinomial(link=logit), data=dfsubs)
summary(fit) # Wald z tests
Anova(fit,type="III") # type III Anova tests

# Proportion of allele expression according to repr status
meanpreds=data.frame(lsmeans::test(lsmeans(fit,~status + colony,type="response"), adjust = "none")) # without multiplicity correction, otherwise use adjust = "mvt"
cints=data.frame(confint(lsmeans(fit,~status + colony,type="response")))[,c(6:7)]
meanpreds=data.frame(meanpreds[,c(1:4)],cints,meanpreds[,c(6:7)])
meanpreds

# Look how proportion of expression changes given status
plot(allEffects(fit), type="response", ylab="Proportion of allele1 expression")


# Model for all genes
overallmeanpreds = do.call(rbind,
                           lapply(levels(df3$geneID),
                                    function (gene) { 
                                      dfsubs = df3[df3$geneID==gene,]
                                      fit = glm(cbind(allele1, allele2) ~ status + colony, 
                                                 family = quasibinomial(link=logit), data=dfsubs)
                                      # overall mean proportion of allele1 expression 
                                      #overallmeanpreds = data.frame(summary(suppressMessages(lsmeans(fit, ~ status + colony))))
                                      #overallmeanpreds$z = overallmeanpreds$lsmean/overallmeanpreds$SE # z scores for deviation from 50:50
                                      #overallmeanpreds$p = 2*pnorm(abs(overallmeanpreds$z),lower.tail = FALSE) # 2 sided p values for deviation from 50:50 (not corrected for mult testing)
                                      #overallmeanpreds$lsmean = plogis(overallmeanpreds$lsmean) # backtransform to original scale
                                      #overallmeanpreds$asymp.LCL = plogis(overallmeanpreds$asymp.LCL)
                                      #overallmeanpreds$asymp.UCL = plogis(overallmeanpreds$asymp.UCL)
                                      #df=df[,-5]
                                       meanpreds=data.frame(lsmeans::test(lsmeans(fit,~status + colony,type="response"), adjust = "none")) # without multiplicity correction, otherwise use adjust = "mvt"
                                      cints=data.frame(confint(lsmeans(fit,~status + colony,type="response")))[,c(6:7)]
                                      meanpreds=data.frame(meanpreds[,c(1:4)],cints,meanpreds[,c(6:7)])
                                      df=data.frame(geneID=gene,meanpreds)
                                      colnames(df)=c("geneID","status","colony","propmatexpr","SE","propmatexpr.LCL","propmatexpr.UCL","z","p")
                                      df$avgpropmatexpr=mean(df$propmatexpr)
                                      df
                                    }
                           ) )

# Benjamini-Hochberg correction for multiple testing
overallmeanpreds$padj = p.adjust(overallmeanpreds$p, method="BH") 

# Add a column showing the significant rows
overallmeanpreds$sig = (overallmeanpreds$padj<0.05)*1

# Pool data across colony
overallmeanpreds$matexprbias = (overallmeanpreds$propmatexpr>0.5)*1
df1 = summaryBy(sig + matexprbias ~ geneID+status, data=overallmeanpreds, FUN=sum)


# Pull out data that is signifiacnt in all colonies for repro / sterile
sig_all_colonies_repr = df1[df1$status=="repro",]$geneID[(df1[df1$status=="repro",]$sig.sum==3&df1[df1$status=="repro",]$matexprbias.sum==3)|
                                                                    (df1[df1$status=="repro",]$sig.sum==3&df1[df1$status=="repro",]$matexprbias.sum==0)]

sig_all_colonies_sterile = df1[df1$status=="non_repro",]$geneID[(df1[df1$status=="non_repro",]$sig.sum==3&df1[df1$status=="non_repro",]$matexprbias.sum==3)|
                                                           (df1[df1$status=="non_repro",]$sig.sum==3&df1[df1$status=="non_repro",]$matexprbias.sum==0)]


length(sig_all_colonies_repr ) # 1819
length(sig_all_colonies_sterile ) # 1832

library(venn)
venn(list(sig_all_colonies_repr,sig_all_colonies_sterile), snames=c("reproductive W","nonreproductive W"), 
     ilab=TRUE, zcolor = "style", cexil=2, cexsn=2)

# Add new columns to show if the row is significant in all colonies for repro/sterile
overallmeanpreds$sig_all_colonies_repr = 0
overallmeanpreds$sig_all_colonies_repr[overallmeanpreds$geneID %in% sig_all_colonies_repr] = 1
overallmeanpreds$sig_all_colonies_sterile = 0
overallmeanpreds$sig_all_colonies_sterile[overallmeanpreds$geneID %in% sig_all_colonies_sterile] = 1

# Pull out the repro/sterile rows which are signigicant on all colonies
signbias_repr = (overallmeanpreds$sig_all_colonies_repr==1)
signbias_nonrepr = (overallmeanpreds$sig_all_colonies_sterile==1)

# Filter these rows with more stringent proportions of expression (70/30%)
signbiasedgenes_repr = unique(overallmeanpreds$geneID[signbias_repr&(overallmeanpreds$avgpropmatexpr>0.7|overallmeanpreds$avgpropmatexpr<0.3)])
length(signbiasedgenes_repr) # 93
signbiasedgenes_nonrepr = unique(overallmeanpreds$geneID[signbias_nonrepr&(overallmeanpreds$avgpropmatexpr>0.7|overallmeanpreds$avgpropmatexpr<0.3)])
length(signbiasedgenes_nonrepr) # 97


venn(list(signbiasedgenes_repr,signbiasedgenes_nonrepr), snames=c("reproductive W","nonreproductive W"), 
     ilab=TRUE, zcolor = "style", cexil=2, cexsn=2) # 92 in common, 1 for repro and 5 for sterile

###------------------------------------------------------------------
# Make plot and output files
###------------------------------------------------------------------

plot(overallmeanpreds$avgpropmatexpr)

significant_repro <- overallmeanpreds[overallmeanpreds$avgpropmatexpr>0.7 & overallmeanpreds$sig_all_colonies_repr == 1 &
                                        overallmeanpreds$status == "repro",]
significant_sterile <- overallmeanpreds[overallmeanpreds$avgpropmatexpr>0.7 & overallmeanpreds$sig_all_colonies_sterile == 1 &
                                        overallmeanpreds$status == "non_repro",]

significant_repro <- significant_repro[significant_repro$colony == "j1",]
significant_sterile <- significant_sterile[significant_sterile$colony == "j1",]

significant_repro <- significant_repro[,-c(3,12,13,14,15)]
significant_sterile <- significant_sterile[,-c(3,12,13,14,15)]

final_sig_list <- rbind(significant_repro, significant_sterile)
#write.table(final_sig_list, file="sig_ASE_repro_and_sterile.txt", col.names = T,
 #           row.names = F, quote = F, sep = '\t')

sig_repro <- as.data.frame(significant_repro$geneID)
colnames(sig_repro) <- "geneID"
sig_sterile <- as.data.frame(significant_sterile$geneID)
colnames(sig_sterile) <- "geneID"

both <- merge(sig_repro, sig_sterile, by="geneID") # 92
unique_sig_sterile <- as.data.frame(sig_sterile[!(sig_sterile$geneID %in% 
                                                    both$geneID ),]) #5
colnames(unique_sig_sterile) <- "geneID"
unique_sig_repro<- as.data.frame(sig_repro[!(sig_repro$geneID %in% 
                                                    both$geneID ),]) #1
colnames(unique_sig_repro) <- "geneID"


overallmeanpreds$gene_status <- "non_sig"
#overallmeanpreds$gene_status[overallmeanpreds$geneID %in% unique_sig_repro$geneID] <- "repro"
#overallmeanpreds$gene_status[overallmeanpreds$geneID %in% unique_sig_sterile$geneID] <- "sterile"
#overallmeanpreds$gene_status[overallmeanpreds$geneID %in% both$geneID] <- "both"
overallmeanpreds$gene_status[overallmeanpreds$geneID %in% final_sig_list$geneID] <- "sig"

plot_data <- overallmeanpreds[!duplicated(overallmeanpreds[,c(1,10,16)]),] #2489

ggplot(plot_data, aes(x=1, y=avgpropmatexpr, colour=gene_status))+
  geom_jitter(size=3, alpha = 0.5)+
  scale_color_manual(breaks=c("non_sig", "sig"),
                     values=c("black","red"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=20),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size=18))+
  ylab("Average Allelic Expression Proportion")+
  xlab("")

## Upset overlap of ASE for repro/nonrepro
head(sig_repro)

sig_repro$Reproductive <- 1
sig_sterile$Sterile <- 1

upset_data <- merge(sig_repro, sig_sterile, by="geneID", all=T)

upset_data[is.na(upset_data)] <- 0

upset(upset_data, order.by = "freq", text.scale = 2)

# Significant overlap? 11030 genes in genome, 97 sterile, 93 repro
phyper(92, 97, 10933, 93, lower.tail = F) # p=6.491171e-226 (YES)

## Plot a few of the most ASE common to repro/sterile

plot_data1 <- overallmeanpreds[overallmeanpreds$geneID %in% final_sig_list$geneID,]
plot_data1_repro <- subset(plot_data1, status == "repro")
plot_data1_sterile <- subset(plot_data1, status == "non_repro")

a1 <- ggplot(plot_data1_repro) +
  aes(x = colony, y = propmatexpr, group = geneID, colour = geneID) +
  geom_point(size = 3)+
  geom_line(size = 2, alpha = 0.4) +
  # xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(axis.text =element_text(size=20),
        axis.title = element_blank(),
        plot.title = element_text(size=20),
        legend.position = "none")+
  ggtitle("Reproductive")+
  ylim(0.5,1.0)

a2 <- ggplot(plot_data1_sterile) +
  aes(x = colony, y = propmatexpr, group = geneID, colour = geneID) +
  geom_point(size = 3)+
  geom_line(size = 2, alpha = 0.4) +
  # xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(axis.text =element_text(size=20),
        axis.title = element_blank(),
        plot.title = element_text(size=20),
        legend.position = "none")+
  ggtitle("Sterile")+
  ylim(0.5,1.0)

all <- ggarrange(a1, a2, ncol=2, nrow=1)

annotate_figure(all, 
                left = text_grob("Average Proportion of Allelic Expression", 
                                 color = "black", rot = 90, size=20),
                bottom = text_grob("Colony", 
                                   color = "black", size =20))




gene <- subset(overallmeanpreds, geneID == "LOC105666882")

a1 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
 # xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=18),
        axis.text =element_text(size=16),
        axis.title = element_blank(),
        plot.title = element_text(size=16))+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC105666882: Mitochondrial 28S ribosomal protein")+
  ylim(0.5,1.0)

gene <- subset(overallmeanpreds, geneID == "LOC100650626")

a2<- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
 # xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=18),
        axis.text =element_text(size=16),
        axis.title = element_blank(),
        plot.title = element_text(size=16))+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100650626: glucosylceramidase-like")+
  ylim(0.5,1.0)


gene <- subset(overallmeanpreds, geneID == "LOC100646400")

a3 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
#  xlab("Colony") + ylab("Average Proportion of Allelic Expression") +  
  theme_bw()+
  theme(legend.text=element_text(size=18),
        axis.text =element_text(size=16),
        axis.title = element_blank(),
        plot.title = element_text(size=16))+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100646400: ionotropic receptor 25a")+
  ylim(0.5,1.0)


top_3 <- ggarrange(a1, a2, a3, ncol=1, nrow=3, common.legend = TRUE, legend="right")

annotate_figure(top_3, 
                left = text_grob("Average Proportion of Allelic Expression", 
                                 color = "black", rot = 90, size=20),
                bottom = text_grob("Colony", 
                                 color = "black", size =20,
                                 hjust = 1.5 ))



# Plot Repro and Sterile Only Genes

# Repro:
# LOC100645175 transcription factor 12

gene <- subset(overallmeanpreds, geneID == "LOC100645175")

p0 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
  xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size=10),
        axis.title = element_blank())+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100645175: transcription factor 12")+
  ylim(0.5,1.0)



# Sterile:
# LOC100643437 probable uridine-cytidine kinase
# LOC100645301 filamin-A
# LOC100648482 antichymotrypsin-2 (a SERPIN)
# LOC100644839  laccase-1
# LOC105665975 golgin subfamily A member 4

gene <- subset(overallmeanpreds, geneID == "LOC100643437")

p3 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
  xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size=10),
        axis.title = element_blank())+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100643437: uridine cytidine kinase")+
  ylim(0.5,1.0)



gene <- subset(overallmeanpreds, geneID == "LOC100645301")

p4 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
  xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size=10),
        axis.title = element_blank())+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100645301: filamin-A")+
  ylim(0.5,1.0)


gene <- subset(overallmeanpreds, geneID == "LOC100648482")

p6 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
  xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size=10),
        axis.title = element_blank())+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100648482: antichymotrypsin-2")+
  ylim(0.5,1.0)

gene <- subset(overallmeanpreds, geneID == "LOC100644839")

p7 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
  xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size=10),
        axis.title = element_blank())+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC100644839: laccase-1")+
  ylim(0.5,1.0)


gene <- subset(overallmeanpreds, geneID == "LOC105665975")

p8 <- ggplot(gene) +
  aes(x = colony, color = status, y = propmatexpr, group = status) +
  geom_point(size = 6)+
  geom_line(size = 2) +
  xlab("Colony") + ylab("Average Proportion of Allelic Expression") + 
  theme_bw()+
  theme(legend.text=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size=10),
        axis.title = element_blank())+
  scale_color_manual("", breaks=c("non_repro","repro"),
                     values = c("#00CC66","#0066CC"),
                     labels=c("Sterile","Reproductive"))+
  ggtitle("LOC105665975: golgin subfamily A4")+
  ylim(0.5,1.0)


unique <- ggarrange(p3, p4, p6, p7, p8, p0, ncol=3, nrow=2, common.legend = T, legend = "right")


annotate_figure(unique, 
                left = text_grob("Average Proportion of Allelic Expression", 
                                 color = "black", rot = 90, size=16),
                bottom = text_grob("Colony", 
                                   color = "black", size =16,
                                   hjust = 1.5 ))


## Plot all showing ASE repro against sterile

head(overallmeanpreds)

overallmeanpreds<- as.data.frame(overallmeanpreds)
new = summaryBy(propmatexpr ~ status + geneID, data=overallmeanpreds, FUN=mean)

casted_data <- dcast(new, geneID  ~ status)

qplot(x=non_repro,
      y=repro,
      data=casted_data,
      col=ifelse(geneID %in% final_sig_list$geneID, I("red"), I("black")),
               # ifelse(geneID %in% diff_meth_genes$geneID, "blue", "black")),
      size=I(3),
      pch = I(16),
      alpha = 0.02,
      xlab = expression(paste("Average Proportion of Alleleic Expression: Sterile")),
      ylab = expression(paste("Average Proportion of Alleleic Expression: Reproductive")),
      geom="point")+
  theme_bw()+
  scale_color_identity() +
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size=18),
        legend.position = "none")

cor.test(casted_data$non_repro, casted_data$repro, method = "spearman")


## See how ASM genes look ... 
nonrepro_amr_genes_unique_list <- read_csv("../ASM/Gene_Lists/nonrepro_amr_genes_unique_list.csv", 
                                           col_names = FALSE)
colnames(nonrepro_amr_genes_unique_list) <- "geneID"

repro_amr_genes_unique_list <- read_csv("../ASM/Gene_Lists/repro_amr_genes_unique_list.csv", 
                                        col_names = FALSE)
colnames(repro_amr_genes_unique_list) <- "geneID"

both_amr_genes_unique_list <- read_csv("../ASM/Gene_Lists/both_amr_genes.csv", 
                                       col_names = FALSE)
colnames(both_amr_genes_unique_list) <- "geneID"

#Start by making an upset plot for ASM alone
nonrepro_amr_genes_unique_list$Sterile <- 1
repro_amr_genes_unique_list$Reproductive <- 1
both_amr_genes_unique_list$Reproductive <- 1
both_amr_genes_unique_list$Sterile <- 1

both_ASM <- merge(nonrepro_amr_genes_unique_list, repro_amr_genes_unique_list,
                  by="geneID", all=T)
all_ASM <- rbind(both_ASM, both_amr_genes_unique_list)

all_ASM[is.na(all_ASM)] <- 0
all_ASM <- all_ASM[!duplicated(all_ASM$geneID),]

upset(all_ASM, text.scale = 2, order.by = "freq")

#Significant? genes in genome 11030, 195 repro, 287 sterile
phyper(97, 195, 10835, 287, lower.tail = F) # p=2.701932e-107 (YES)

# Now have a look at ASM and ASE overlap for repro/sterile
head(all_ASM)
head(upset_data)

colnames(all_ASM) <- c("geneID","Repro_AMG", "Sterile_AMG")
colnames(upset_data) <- c("geneID", "Repro_AEG","Sterile_AEG")

ASM_ASE <- merge(all_ASM, upset_data, by="geneID", all=T)
ASM_ASE[is.na(ASM_ASE)] <- 0

ASM_ASE <- ASM_ASE[!duplicated(ASM_ASE),]
upset(ASM_ASE, text.scale = 2, order.by = "freq")

# Take a look at the 6 genes showing ASE and ASM in repro/sterile and the one
# in repro AEG, sterile AEG and repro AMG

in_all <- ASM_ASE[(ASM_ASE$Repro_AMG == 1 & ASM_ASE$Sterile_AMG == 1 &
                    ASM_ASE$Repro_AEG == 1 & ASM_ASE$Sterile_AEG ==1),]
# LOC100643777 40S ribosomal protein S6
# LOC100643807 protein Jumonji
# LOC100645564 peripheral plasma membrane protein CASK
# LOC100651252 ATP-dependent RNA helicase WM6
# LOC100652132 importin-11
# LOC105666051 uncharacterized

# Significant? 11030 genes in genome, 98 ASE and 385 ASM and 6 overlap
phyper(6, 98, 10932, 385, lower.tail = F) # p=0.05519901 (NO)

# Take a look at one gen overlap
other_one <- ASM_ASE[(ASM_ASE$Repro_AMG == 1 & ASM_ASE$Sterile_AMG == 0 &
                     ASM_ASE$Repro_AEG == 1 & ASM_ASE$Sterile_AEG ==1),]

# LOC105666711 tyrosine-protein kinase Btk29A (immunity)

# Significant? 11030 genes in genome, 98 ASE and 195 ASM and 1 overlap
sum(all_ASM$Repro_AMG)
phyper(1, 98, 10932, 195, lower.tail = F) # p=0.5197835 (NO)



### Make some kind of plot of the ASM genes ASE levels ... 


head(overallmeanpreds)
head(all_ASM)


ASE_levels_of_ASM <- merge(all_ASM, overallmeanpreds, by="geneID")
ASE_levels_of_ASM <- ASE_levels_of_ASM[,c(1,2,3,4,12,14)]
ASE_levels_of_ASM <- ASE_levels_of_ASM[!duplicated(ASE_levels_of_ASM),]

ASE_levels_of_ASM$Repro_AMG[ASE_levels_of_ASM$Repro_AMG == 1] <- ASE_levels_of_ASM$avgpropmatexpr
ASE_levels_of_ASM$Sterile_AMG[ASE_levels_of_ASM$Sterile_AMG == 1] <- ASE_levels_of_ASM$avgpropmatexpr

ASE_levels_of_ASM$Repro_AMG[ASE_levels_of_ASM$Repro_AMG <0.5 ] <- NA
ASE_levels_of_ASM$Sterile_AMG[ASE_levels_of_ASM$Sterile_AMG <0.5 ] <- NA

ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "repro" 
                         & !is.na(ASE_levels_of_ASM$Repro_AMG)
                         & !is.na(ASE_levels_of_ASM$Sterile_AMG)] <- "repro_AMG_both"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "repro" 
                              & !is.na(ASE_levels_of_ASM$Repro_AMG)
                              & is.na(ASE_levels_of_ASM$Sterile_AMG)] <- "repro_AMG_repro"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "repro" 
                              & is.na(ASE_levels_of_ASM$Repro_AMG)
                              & !is.na(ASE_levels_of_ASM$Sterile_AMG)] <- "repro_AMG_sterile"

ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "non_repro" 
                              & !is.na(ASE_levels_of_ASM$Repro_AMG)
                              & !is.na(ASE_levels_of_ASM$Sterile_AMG)] <- "sterile_AMG_both"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "non_repro" 
                              & !is.na(ASE_levels_of_ASM$Repro_AMG)
                              & is.na(ASE_levels_of_ASM$Sterile_AMG)] <- "sterile_AMG_repro"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "non_repro" 
                              & is.na(ASE_levels_of_ASM$Repro_AMG)
                              & !is.na(ASE_levels_of_ASM$Sterile_AMG)] <- "sterile_AMG_sterile"


ggplot(ASE_levels_of_ASM, aes(x=plot_column, y=avgpropmatexpr, fill=status))+
  geom_boxplot(outlier.colour = "red")+
  geom_point(position = position_jitter(w = 0.2, h = 0))+
  theme_bw()+
  xlab("Allelically Methylated Gene Category")+
  ylab("Average Proportion of Allelic Expression")+
  scale_fill_manual(breaks=c("repro","non_repro"),
                    labels=c("Reproductive", "Sterile"),
                    values=c("#00CC66","#0066CC"),
                    name="Status")+
  scale_x_discrete(limits=c("repro_AMG_both", "sterile_AMG_both", "repro_AMG_repro",
                            "sterile_AMG_repro", "repro_AMG_sterile","sterile_AMG_sterile"),
                   labels=c("Both","Both", "Reproductive","Reproductive","Sterile","Sterile"))+
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18))

model_data <- ASE_levels_of_ASM 
model_data$plot_column <- gsub("repro_", "", model_data$plot_column)
model_data$plot_column <- gsub("sterile_", "", model_data$plot_column)



#Above graph possibly not the one I want ...

# The difference in ASE between castes plotted for three ASM categories
#head(model_data)
#head(all_ASM) #385
#tail(new)

#ase_data_by_caste <- dcast(new,  geneID ~ status)
# + = higher proper ASE in repro
#ase_data_by_caste$ase_diff <- ase_data_by_caste$repro - ase_data_by_caste$non_repro

#new_plot_data <- merge(all_ASM, ase_data_by_caste, by ="geneID", all.x = T)
#new_plot_data$ase_diff[is.na(new_plot_data$ase_diff)] <- 0


## Stats here...
model_data$plot_column <- as.factor(model_data$plot_column)

model1<-lm(avgpropmatexpr~plot_column*status, data=model_data)
model2<-lm(avgpropmatexpr~plot_column+status, data=model_data)
anova(model1,model2)
#   Res.Df     RSS Df   Sum of Sq      F Pr(>F) (no interaction between ASM state and repro status of ASE gene)
#1    253 0.68383                             
#2    255 0.68389 -2 -6.1683e-05 0.0114 0.9887

# Can now just do a Kruskal-wallis and dunn test
kruskal.test(model_data$avgpropmatexpr ~ model_data$plot_column) 
# Kruskal-Wallis chi-squared = 30.103, df = 2, p-value = 2.906e-07

library(FSA)
dunnTest(model_data$avgpropmatexpr ~ model_data$plot_column,
         method="bh")
#                Comparison         Z      P.unadj        P.adj
#1    AMG_both - AMG_repro  5.445310 5.171528e-08 1.551458e-07
#2  AMG_both - AMG_sterile  3.621866 2.924860e-04 4.387289e-04
#3 AMG_repro - AMG_sterile -2.688114 7.185684e-03 7.185684e-03




# Get gene lists on there own of ASE genes for GO analysis
repro <- subset(sig_ASE_repro_and_sterile, status == "repro")
sterile <- subset(sig_ASE_repro_and_sterile, status == "non_repro")

both <- merge(repro, sterile, by = "geneID")
both <- as.data.frame(both$geneID)
write.table(both, file="ASE_genes_in_both.txt", col.names = F, row.names = F, sep = "\t", quote = F)

repro_only <- as.data.frame(repro$geneID[!(repro$geneID %in% both$`both$geneID`)])
sterile_only <- as.data.frame(sterile$geneID[!(sterile$geneID %in% both$`both$geneID`)])

write.table(repro_only, file="ASE_genes_in_repro_only.txt", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(sterile_only, file="ASE_genes_in_sterile_only.txt", col.names = F, row.names = F, sep = "\t", quote = F)