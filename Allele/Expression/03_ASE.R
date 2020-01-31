###-------------------------------------------------------------------
# BINOMIAL GLMs TO DETERMINE ASE
###-------------------------------------------------------------------

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/MERN_Analysis/ASE")
library(pbapply)
library(plyr)
library(doBy)
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
 range(df$prop.allele1) # 0.5000000 0.984127
 hist(df$prop.allele1,col="steelblue",
      xlab="Proportion of allele1 expression", main = NULL,
      cex.lab=1.5, cex.axis=1.5)
 
 range(df$prop.allele2) # 0.01587302 0.500000000
 hist(df$prop.allele2,col="steelblue",
      xlab="Proportion of allele2 expression",main = NULL,
      cex.lab=1.5, cex.axis=1.5)

###------------------------------------------------------------------
# Filter data to keep only genes in all crosses/conditions
###------------------------------------------------------------------


# See which genes occur in all experimental conditions/crosses
df <- as.data.frame(df)
df1 = summaryBy(allele1.sum+allele2.sum ~ geneID, data=df, FUN=length)

# Keep only genes for which we have counts in all 6 conditions (repro/sterile x3 colonies)
# This leave only 1712 genes, need to be less stringent
#selectedgenes = df1$geneID[df1$allele1.sum.length == 6]

selectedgenes = df1$geneID[df1$allele1.sum.length > 4] 

#length(selectedgenes) # 2673 (BETTER) # This thows up issues later, try to resolve, in model it needs factors with 2 or more lvls
df2 = df[(df$geneID %in% selectedgenes),]

# Remove colony and pool by status (as missing one from colony 8? J8_24)
df3 = summaryBy(allele1.sum+allele2.sum ~ geneID + status + colony, data=df2, FUN=sum)
df3 <- df3[!is.na(df3),]

###------------------------------------------------------------------
# Run model
###------------------------------------------------------------------

set_sum_contrasts()
head(df3)
df3 <- as.data.frame(df3)
length(unique(df3$geneID)) # 1712
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
meanpreds=data.frame(emmeans::test(lsmeans(fit,~status + colony,type="response"), adjust = "none")) # without multiplicity correction, otherwise use adjust = "mvt"
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
                                       meanpreds=data.frame(emmeans::test(lsmeans(fit,~status + colony,type="response"), adjust = "none")) # without multiplicity correction, otherwise use adjust = "mvt"
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


length(sig_all_colonies_repr ) # 1708
length(sig_all_colonies_sterile ) # 1644

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

# Filter these rows with more stringent proportions of expression
signbiasedgenes_repr = unique(overallmeanpreds$geneID
                              [signbias_repr&(overallmeanpreds$avgpropmatexpr>0.65|overallmeanpreds$avgpropmatexpr<0.35)])
length(signbiasedgenes_repr) # 131
signbiasedgenes_nonrepr = unique(overallmeanpreds$geneID
                                 [signbias_nonrepr&(overallmeanpreds$avgpropmatexpr>0.65|overallmeanpreds$avgpropmatexpr<0.35)])
length(signbiasedgenes_nonrepr) # 124

hist(overallmeanpreds$propmatexpr[overallmeanpreds$sig==1])
hist(overallmeanpreds$propmatexpr[overallmeanpreds$sig_all_colonies_sterile==1])

cut_off <- c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
num_genes_repro <- c(1708,1529,534,131,58, 23,8,3,1,0)
num_genes_sterile <- c(1644,1481,538,125,59,25,9,4,2,0)

significance_threshold <- data.frame(cut_off, num_genes_repro, num_genes_sterile)
significance_threshold$total <- significance_threshold$num_genes_repro + significance_threshold$num_genes_sterile

plot(significance_threshold$cut_off ~ significance_threshold$total)

ggplot(significance_threshold, aes(x=cut_off, y=total))+
  geom_point(size=3, alpha = 0.5)+
  geom_line()+
  theme_bw()+
  theme(axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.x  = element_text(size=16),
        axis.text.y = element_text(size=16))+
  xlab("Significance threshold (proportion alleleic expression)")+
  ylab("Number of significant genes")+
  geom_vline(xintercept = 0.65, color = "red", linetype="dashed") 


venn(list(signbiasedgenes_repr,signbiasedgenes_nonrepr), snames=c("reproductive W","nonreproductive W"), 
     ilab=TRUE, zcolor = "style", cexil=2, cexsn=2) 
###------------------------------------------------------------------
# Make plot and output files
###------------------------------------------------------------------

plot(overallmeanpreds$avgpropmatexpr)

significant_repro <- overallmeanpreds[overallmeanpreds$avgpropmatexpr>0.65 & overallmeanpreds$sig_all_colonies_repr == 1 &
                                        overallmeanpreds$status == "repro",]
significant_sterile <- overallmeanpreds[overallmeanpreds$avgpropmatexpr>0.65 & overallmeanpreds$sig_all_colonies_sterile == 1 &
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

both <- merge(sig_repro, sig_sterile, by="geneID") # 116
unique_sig_sterile <- as.data.frame(sig_sterile[!(sig_sterile$geneID %in% 
                                                    both$geneID ),]) #8
colnames(unique_sig_sterile) <- "geneID"
unique_sig_repro<- as.data.frame(sig_repro[!(sig_repro$geneID %in% 
                                                    both$geneID ),]) #15
colnames(unique_sig_repro) <- "geneID"


head(both)
write.table(both, file="ASE_both_gene_list.txt", col.names = T,
                     row.names = F, quote = F, sep = '\t')
write.table(unique_sig_sterile, file="ASE_sterile_gene_list.txt", col.names = T,
            row.names = F, quote = F, sep = '\t')
write.table(unique_sig_repro, file="ASE_repro_gene_list.txt", col.names = T,
            row.names = F, quote = F, sep = '\t')

all_ASE_genes <- rbind(both, unique_sig_sterile, unique_sig_repro)
write.table(all_ASE_genes, file="ASE_all_gene_list.txt", col.names = T,
            row.names = F, quote = F, sep = '\t')


overallmeanpreds$gene_status <- "non_sig"
overallmeanpreds$gene_status[overallmeanpreds$geneID %in% final_sig_list$geneID] <- "sig"

plot_data <- overallmeanpreds[!duplicated(overallmeanpreds[,c(1,10,16)]),] #2673

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

upset(upset_data, order.by = "freq", text.scale = 4, point.size = 5)

# Significant overlap? 11030 genes in genome, 97 sterile, 93 repro
phyper(116, 124, 10899, 131, lower.tail = F) # p=1.512047e-251 (YES)


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


#### Plot genes unique to either repro or sterile
####
####
####

head(unique_sig_repro)
head(unique_sig_sterile)

repro_unique_data <- overallmeanpreds[overallmeanpreds$geneID %in% unique_sig_repro$geneID,]
sterile_unique_data <- overallmeanpreds[overallmeanpreds$geneID %in% unique_sig_sterile$geneID,]

# Find the top different between status'
genes <- as.vector(unique_sig_repro$geneID)
pdf("repro_unique_genes.pdf")
for (i in genes){
  print(ggplot(subset(overallmeanpreds, geneID == i)) +
          aes(x = colony, color = status, y = propmatexpr, group = status) +
          geom_point(size = 6)+
          geom_line(size = 2)+
          theme_bw()+
          theme(legend.text=element_text(size=18),
                axis.text =element_text(size=16),
                axis.title = element_blank(),
                plot.title = element_text(size=16))+
          scale_color_manual("", breaks=c("non_repro","repro"),
                             values = c("#00CC66","#0066CC"),
                             labels=c("Sterile","Reproductive"))+
          ggtitle(print(i))+
          ylim(0.5,1.0))
}
dev.off()

genes <- as.vector(unique_sig_sterile$geneID)
pdf("sterile_unique_genes.pdf")
for (i in genes){
  print(ggplot(subset(overallmeanpreds, geneID == i)) +
          aes(x = colony, color = status, y = propmatexpr, group = status) +
          geom_point(size = 6)+
          geom_line(size = 2)+
          theme_bw()+
          theme(legend.text=element_text(size=18),
                axis.text =element_text(size=16),
                axis.title = element_blank(),
                plot.title = element_text(size=16))+
          scale_color_manual("", breaks=c("non_repro","repro"),
                             values = c("#00CC66","#0066CC"),
                             labels=c("Sterile","Reproductive"))+
          ggtitle(print(i))+
          ylim(0.5,1.0))
}
dev.off()


# Plot top diff for repro and sterile (rest in supplementary)
# Repro
# LOC100644413 UNC93-like protein MFSD11
# LOC100645303 phosphatidylinositide phosphatase SAC2
# LOC100647110 voltage-dependent calcium channel subunit alpha-2/delta-3


a <- ggplot(subset(overallmeanpreds, geneID == "LOC100644413")) +
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
  ggtitle("LOC100644413")+
  ylim(0.5,1.0)

b <- ggplot(subset(overallmeanpreds, geneID == "LOC100645303")) +
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
  ggtitle("LOC100645303")+
  ylim(0.5,1.0)

c <- ggplot(subset(overallmeanpreds, geneID == "LOC100647110")) +
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
  ggtitle("LOC100647110")+
  ylim(0.5,1.0)


# Sterile
# LOC100647178 venom acid phosphatase Acph-1
# LOC100646178 flavin reductase (NADPH)
# LOC100643395 integral membrane protein GPR180

a1 <- ggplot(subset(overallmeanpreds, geneID == "LOC100647178")) +
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
  ggtitle("LOC100647178")+
  ylim(0.5,1.0)

b1 <- ggplot(subset(overallmeanpreds, geneID == "LOC100646178")) +
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
  ggtitle("LOC100646178")+
  ylim(0.5,1.0)

c1 <- ggplot(subset(overallmeanpreds, geneID == "LOC100643395")) +
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
  ggtitle("LOC100643395")+
  ylim(0.5,1.0)



top_3_diff <- ggarrange(a, b, c, a1, b1, c1, 
                        ncol=3, nrow=2, common.legend = TRUE, legend="right")

annotate_figure(top_3_diff, 
                left = text_grob("Average Proportion of Allelic Expression", 
                                 color = "black", rot = 90, size=20),
                bottom = text_grob("Colony", 
                                 color = "black", size =20,
                                 hjust = 1.5 ))



# Plot general top three showing ASE
head(final_sig_list)

# LOC105666882: uncharacterized (also only in sterile... weird)
# LOC100650626: glucosylceramidase-like
# LOC100642453: heparan sulfate 2-O-sulfotransferase pipe
# LOC100646568: DNA-directed RNA polymerase III subunit RPC7-like


a1 <- ggplot(subset(overallmeanpreds, geneID == "LOC100650626")) +
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

b1 <- ggplot(subset(overallmeanpreds, geneID == "LOC100642453")) +
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
  ggtitle("LOC100642453: heparan sulfate 2-O-sulfotransferase")+
  ylim(0.5,1.0)

c1 <- ggplot(subset(overallmeanpreds, geneID == "LOC100646568")) +
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
  ggtitle("LOC100646568: DNA-directed RNA polymerase II")+
  ylim(0.5,1.0)



top_3_diff <- ggarrange(a1, b1, c1, 
                        ncol=1, nrow=3, common.legend = TRUE, legend="right")

annotate_figure(top_3_diff, 
                left = text_grob("Average Proportion of Allelic Expression", 
                                 color = "black", rot = 90, size=20),
                bottom = text_grob("Colony", 
                                   color = "black", size =20,
                                   hjust = 1.5 ))



## Plot all showing ASE repro against sterile DONE

head(overallmeanpreds)

overallmeanpreds<- as.data.frame(overallmeanpreds)
overallmeanpreds<- overallmeanpreds[overallmeanpreds$propmatexpr >= 0.5,]

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
# S = 1229363078, p-value < 2.2e-16, rho = 0.6137801 

# ------------------ ASM and Paired with ASE -------------------------
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

upset(all_ASM, text.scale = 4, order.by = "freq", point.size = 4)

#Significant? genes in genome 11030, 195 repro, 287 sterile
phyper(97, 195, 10835, 287, lower.tail = F) # p=2.701932e-107 (YES)

# Now have a look at ASM and ASE overlap for repro/sterile
head(all_ASM)
head(upset_data)

colnames(all_ASM) <- c("geneID","Repro AMG", "Sterile AMG")
colnames(upset_data) <- c("geneID", "Repro AEG","Sterile AEG")

ASM_ASE <- merge(all_ASM, upset_data, by="geneID", all=T)
ASM_ASE[is.na(ASM_ASE)] <- 0

ASM_ASE <- ASM_ASE[!duplicated(ASM_ASE),]
upset(ASM_ASE, text.scale = 2.5, order.by = "freq",point.size = 4)

# Take a look at the 6 genes showing ASE and ASM in repro/sterile and the one
# in repro AEG, sterile AEG and repro AMG plus one with sterile and repro AEG
# and sterile AMG

in_all <- ASM_ASE[(ASM_ASE$`Repro AMG` == 1 & ASM_ASE$`Sterile AMG` == 1 &
                    ASM_ASE$`Repro AEG` == 1 & ASM_ASE$`Sterile AEG` ==1),]
# LOC100643777 40S ribosomal protein S6
# LOC100643941 connectin
# LOC100644811 neuroligin-4, Y-linked 
# LOC100644932 AP-1 complex subunit mu-1
# LOC100652132 importin-11
# LOC105665778 regulator of microtubule dynamics protein 1-like


# Significant? 11030 genes in genome, 98 ASE and 385 ASM and 6 overlap
phyper(6, 139, 10932, 385, lower.tail = F) # p=0.2099091 (NO)

# Take a look at one gen overlap
other_one <- ASM_ASE[(ASM_ASE$`Repro AMG` == 1 & ASM_ASE$`Sterile AMG` == 0 &
                        ASM_ASE$`Repro AEG` == 1 & ASM_ASE$`Sterile AEG` ==1),]

# LOC105666711 tyrosine-protein kinase Btk29A (immunity)

other_one_again <- ASM_ASE[(ASM_ASE$`Repro AMG` == 0 & ASM_ASE$`Sterile AMG` == 1 &
                        ASM_ASE$`Repro AEG` == 1 & ASM_ASE$`Sterile AEG` ==1),]

# LOC100643219 putative pre-mRNA-splicing factor ATP-dependent RNA helicase PRP

### Make some kind of plot of the ASM genes ASE levels ... 


head(overallmeanpreds)
head(all_ASM)


ASE_levels_of_ASM <- merge(all_ASM, overallmeanpreds, by="geneID")
ASE_levels_of_ASM <- ASE_levels_of_ASM[,c(1,2,3,4,12,14)]
ASE_levels_of_ASM <- ASE_levels_of_ASM[!duplicated(ASE_levels_of_ASM),]

ASE_levels_of_ASM$`Repro AMG`[ASE_levels_of_ASM$`Repro AMG` == 1] <- ASE_levels_of_ASM$avgpropmatexpr
ASE_levels_of_ASM$`Sterile AMG`[ASE_levels_of_ASM$`Sterile AMG` == 1] <- ASE_levels_of_ASM$avgpropmatexpr

ASE_levels_of_ASM$`Repro AMG`[ASE_levels_of_ASM$`Repro AMG` <0.5 ] <- NA
ASE_levels_of_ASM$`Sterile AMG`[ASE_levels_of_ASM$`Sterile AMG` <0.5 ] <- NA

ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "repro" 
                         & !is.na(ASE_levels_of_ASM$`Repro AMG`)
                         & !is.na(ASE_levels_of_ASM$`Sterile AMG`)] <- "repro_AMG_both"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "repro" 
                              & !is.na(ASE_levels_of_ASM$`Repro AMG`)
                              & is.na(ASE_levels_of_ASM$`Sterile AMG`)] <- "repro_AMG_repro"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "repro" 
                              & is.na(ASE_levels_of_ASM$`Repro AMG`)
                              & !is.na(ASE_levels_of_ASM$`Sterile AMG`)] <- "repro_AMG_sterile"

ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "non_repro" 
                              & !is.na(ASE_levels_of_ASM$`Repro AMG`)
                              & !is.na(ASE_levels_of_ASM$`Sterile AMG`)] <- "sterile_AMG_both"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "non_repro" 
                              & !is.na(ASE_levels_of_ASM$`Repro AMG`)
                              & is.na(ASE_levels_of_ASM$`Sterile AMG`)] <- "sterile_AMG_repro"
ASE_levels_of_ASM$plot_column[ASE_levels_of_ASM$status == "non_repro" 
                              & is.na(ASE_levels_of_ASM$`Repro AMG`)
                              & !is.na(ASE_levels_of_ASM$`Sterile AMG`)] <- "sterile_AMG_sterile"


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
#1    294 0.57671                            
#2    296 0.57713 -2 -0.00042921 0.1094 0.8964

# Can now just do a Kruskal-wallis and dunn test
kruskal.test(model_data$avgpropmatexpr ~ model_data$plot_column) 
# Kruskal-Wallis chi-squared = 28.838, df = 2, p-value = 5.469e-07

library(FSA)
dunnTest(model_data$avgpropmatexpr ~ model_data$plot_column,
         method="bh")
#                Comparison         Z      P.unadj        P.adj
#1    AMG_both - AMG_repro  5.149837 2.607136e-07 7.821408e-07
#2  AMG_both - AMG_sterile  4.147296 3.364246e-05 5.046369e-05
#3 AMG_repro - AMG_sterile -1.851454 6.410429e-02 6.410429e-02

