
###差异OTU分析
setwd("F:\\#R\\Others\\XingXF\\R_markdown_for_published")

library(ggplot2)
library(vegan)
library(reshape2)
library(Rmisc)
#BiocManager::install("edgeR")
library(edgeR)

##Part1, the mfec part
###(1) load the data
###design
design <- read.csv(file = "./data/mfec_metadata.csv", header = T)
names(design)

###all dataset
natural_total <- subset(design, Design == "Natural_total" & Bacteria >= 1000)
natural_endo <- subset(design, Design == "Natural_endo" & Bacteria >= 1000)
syncom_endo <- subset(design, Design == "Syncom_endo" & Bacteria >= 1000)

###natural_total_d5 & natural_endo_d5
natural_total_d5 <- subset(natural_total, Days == "D5" & Experiment != "I1")
natural_endo_d5 <- subset(natural_endo,  Days == "D5")

###syncom_endo_d5
syncom_endo_d5 <- subset(syncom_endo,  Days == "D6" & Experiment != "I2")
syncom_endo_mfec <- subset(syncom_endo_d5, Genotype == "Col0" | Genotype == "mfec")
syncom_endo_fec <- subset(syncom_endo_d5, Genotype == "Col0" | Genotype == "fec")
syncom_endo_min7 <- subset(syncom_endo_d5, Genotype == "Col0" | Genotype == "min7")

###OTU table & taxonomy & fasta
otu_table <- read.csv(file = "./data/asv_table.csv", header = T, row.names = 2) ###this is the rarefied OTU/ASV table

otu_taxonomy <- read.csv(file = "./data/asv_taxonomy.csv", header = T)
otu_taxonomy <- subset(otu_taxonomy, otu_taxonomy$Feature_ID %in% otu_table$Feature_ID)
match(otu_taxonomy$Feature_ID, otu_table$Feature_ID)
row.names(otu_taxonomy) <- row.names(otu_table)

otu_fasta <- read.csv(file = "./data/asv_fasta.csv", header = T)
otu_fasta <- subset(otu_fasta, otu_fasta$Feature_ID %in% otu_table$Feature_ID)
match(otu_fasta$Feature_ID, otu_fasta$Feature_ID)
row.names(otu_fasta) <- row.names(otu_table)

###subset OTUs/ASVs
otu_table_natural_total <- otu_table[, names(otu_table) %in% natural_total_d5$sample_id]
otu_table_natural_endo <- otu_table[, names(otu_table) %in% natural_endo_d5$sample_id]

otu_table_syncom_endo_mfec <- otu_table[, names(otu_table) %in% syncom_endo_mfec$sample_id]
otu_table_syncom_endo_fec <- otu_table[, names(otu_table) %in% syncom_endo_fec$sample_id]
otu_table_syncom_endo_min7 <- otu_table[, names(otu_table) %in% syncom_endo_min7$sample_id]

group_natural_total <- as.character(natural_total_d5$Genotype) 
group_natural_total
group_natural_endo <- as.character(natural_endo_d5$Genotype) 
group_natural_endo

group_syncom_endo_mfec <- as.character(syncom_endo_mfec$Genotype)
group_syncom_endo_mfec
group_syncom_endo_fec <- as.character(syncom_endo_fec$Genotype)
group_syncom_endo_fec
group_syncom_endo_min7 <- as.character(syncom_endo_min7$Genotype)
group_syncom_endo_min7


###(2) Differentially abundant OTUs were detected using EdgeR’s generalized linear model (GLM) approach

###Creates a DGEList object from a table of counts
dgelist_natural_total <- DGEList(counts = otu_table_natural_total, group = group_natural_total)
dgelist_natural_endo <- DGEList(counts = otu_table_natural_endo, group = group_natural_endo)

dgelist_syncom_endo_mfec <- DGEList(counts = otu_table_syncom_endo_mfec, group = group_syncom_endo_mfec)
dgelist_syncom_endo_fec <- DGEList(counts = otu_table_syncom_endo_fec, group = group_syncom_endo_fec)
dgelist_syncom_endo_min7 <- DGEList(counts = otu_table_syncom_endo_min7, group = group_syncom_endo_min7)

###Remove low abundance OTUs
###低表达量的otu受到随机因素影响比较大，故其统计结果也不可靠，还会影响p值校正过程
dgelist_natural_total <- dgelist_natural_total[rowSums(dgelist_natural_total$counts) >= 10, ]
dgelist_natural_endo <- dgelist_natural_endo[rowSums(dgelist_natural_endo$counts) >= 10, ]

dgelist_syncom_endo_mfec <- dgelist_syncom_endo_mfec[rowSums(dgelist_syncom_endo_mfec$counts) >= 10, ]
dgelist_syncom_endo_fec <- dgelist_syncom_endo_fec[rowSums(dgelist_syncom_endo_fec$counts) >= 10, ]
dgelist_syncom_endo_min7 <- dgelist_syncom_endo_min7[rowSums(dgelist_syncom_endo_min7$counts) >= 10, ]

###plot MDS
###This results are same to pcoa plots described before
plotMDS(dgelist_natural_total, col = ifelse(natural_total_d5$Genotype == "Col0", "red", "darkgreen") , dim = c(1, 2))
plotMDS(dgelist_natural_endo, col = ifelse(natural_endo_d5$Genotype == "Col0", "red", "darkgreen") , dim = c(1, 2))

plotMDS(dgelist_syncom_endo_mfec, col = ifelse(syncom_endo_mfec$Genotype == "Col0", "red", "darkgreen"), dim = c(1, 2))
plotMDS(dgelist_syncom_endo_fec, col = ifelse(syncom_endo_fec$Genotype == "Col0", "red", "darkgreen"), dim = c(1, 2))
plotMDS(dgelist_syncom_endo_min7, col = ifelse(syncom_endo_min7$Genotype == "Col0", "red", "darkgreen"), dim = c(1, 2))

###model.matrix
design_natural_total <- model.matrix(~group_natural_total)
design_natural_endo <- model.matrix(~group_natural_endo)

design_syncom_endo_mfec <- model.matrix(~group_syncom_endo_mfec)
design_syncom_endo_fec <- model.matrix(~group_syncom_endo_fec)
design_syncom_endo_min7 <- model.matrix(~group_syncom_endo_min7)

###Estimate dispersions
dge_natural_total <- estimateDisp(dgelist_natural_total, design_natural_total, robust = TRUE)
dge_natural_endo <- estimateDisp(dgelist_natural_endo, design_natural_endo, robust = TRUE)

dge_syncom_endo_mfec <- estimateDisp(dgelist_syncom_endo_mfec, design_syncom_endo_mfec, robust = TRUE) 
dge_syncom_endo_fec <- estimateDisp(dgelist_syncom_endo_fec, design_syncom_endo_fec, robust = TRUE) 
dge_syncom_endo_min7 <- estimateDisp(dgelist_syncom_endo_min7, design_syncom_endo_min7, robust = TRUE) 

###Natural_total GLM: negative binomial generalized log-linear model
fit_natural_total <- glmFit(dge_natural_total, design_natural_total, robust = TRUE)
lrt_natural_total <- glmLRT(fit_natural_total)

topTags(lrt_natural_total, adjust.method = "fdr", sort.by = "PValue")

###Collect all the informations and save the data
###Relative abundance, GLM stats, taxonomy, fasta
natural_total_glmLRT_otu <- cbind(Col0 = rowMeans(dgelist_natural_total$counts[,dgelist_natural_total$samples$group == "Col0"])/1000, 
                                  mfec = rowMeans(dgelist_natural_total$counts[,dgelist_natural_total$samples$group == "mfec"])/1000,
                                  as.data.frame(topTags(lrt_natural_total, n = dim(dgelist_natural_total$counts)[1], adjust.method = "fdr",sort.by = "none")))

natural_total_glmLRT_otu_tax_fasta <- merge.data.frame(natural_total_glmLRT_otu, cbind(otu_taxonomy, Fasta = otu_fasta[,-1]), by="row.names")

###sig used for pchs
natural_total_glmLRT_otu_tax_fasta$sig <- as.factor(ifelse(natural_total_glmLRT_otu_tax_fasta$FDR >= 0.05,"NS",
                                          ifelse(natural_total_glmLRT_otu_tax_fasta$logFC > 0,"Enriched","Depleted")))

###phylum used for colors :"Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others"

natural_total_glmLRT_otu_tax_fasta$color <- as.factor(ifelse(natural_total_glmLRT_otu_tax_fasta$Phylum %in% c("P__Firmicutes", "P__Actinobacteria", "P__Bacteroidetes"), 
                                                             as.character(natural_total_glmLRT_otu_tax_fasta$Phylum),
                                                             ifelse(natural_total_glmLRT_otu_tax_fasta$Class %in% c("C__Alphaproteobacteria", "C__Betaproteobacteria", "C__Gammaproteobacteria"),
                                                             as.character(natural_total_glmLRT_otu_tax_fasta$Class), "P__Others")))

###the lowest level used for sig asv labels
natural_total_glmLRT_otu_tax_fasta$label <- as.factor(ifelse(natural_total_glmLRT_otu_tax_fasta$FDR >= 0.05,"",
                                                             ifelse(natural_total_glmLRT_otu_tax_fasta$Genus != "", as.character(natural_total_glmLRT_otu_tax_fasta$Genus),
                                                                    ifelse(natural_total_glmLRT_otu_tax_fasta$Family != "", as.character(natural_total_glmLRT_otu_tax_fasta$Family), 
                                                                           ifelse(natural_total_glmLRT_otu_tax_fasta$Order != "", as.character(natural_total_glmLRT_otu_tax_fasta$Order),
                                                                                  ifelse(natural_total_glmLRT_otu_tax_fasta$Class != "", as.character(natural_total_glmLRT_otu_tax_fasta$Class), 
                                                                                  as.character(natural_total_glmLRT_otu_tax_fasta$Phylum)))))))

write.csv(natural_total_glmLRT_otu_tax_fasta, './result/natural_total_glmLRT_otu_tax_fasta.csv', quote = FALSE)  

dge_de_natural_total <- decideTestsDGE(lrt_natural_total, adjust.method = 'fdr', p.value = 0.05)
summary(dge_de_natural_total)

###Optimized Visualization using ggplot2
windows(16,9)
ggplot(natural_total_glmLRT_otu_tax_fasta, aes(logCPM, logFC, shape=sig, color=color, fill=color)) +
  geom_point(cex=3+5*(natural_total_glmLRT_otu_tax_fasta$Col0+natural_total_glmLRT_otu_tax_fasta$mfec))+
  scale_shape_manual(values=c(21)) +
  scale_color_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  scale_fill_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="black", linetype="dashed") +
  geom_text(aes(label = label), colour="grey50", vjust = 1.5, cex=3+5*(natural_total_glmLRT_otu_tax_fasta$Col0+natural_total_glmLRT_otu_tax_fasta$mfec)) +
  labs(x=paste("Log10CPM"),y=paste("Log10FC")) +
  ggtitle(paste("mfec vs Col-0")) +
  theme_bw() + 
  theme(title = element_text(size = 10, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))


###Natural_endo GLM: negative binomial generalized log-linear model
fit_natural_endo <- glmFit(dge_natural_endo, design_natural_endo, robust = TRUE)
lrt_natural_endo <- glmLRT(fit_natural_endo)

topTags(lrt_natural_endo, adjust.method = "fdr", sort.by = "PValue")

###Collect all the informations and save the data
###Relative abundance, GLM stats, taxonomy, fasta
natural_endo_glmLRT_otu <- cbind(Col0 = rowMeans(dgelist_natural_endo$counts[,dgelist_natural_endo$samples$group == "Col0"])/1000, 
                                 mfec = rowMeans(dgelist_natural_endo$counts[,dgelist_natural_endo$samples$group == "mfec"])/1000,
                                 as.data.frame(topTags(lrt_natural_endo, n = dim(dgelist_natural_endo$counts)[1], adjust.method = "fdr",sort.by = "none")))

natural_endo_glmLRT_otu_tax_fasta <- merge.data.frame(natural_endo_glmLRT_otu, cbind(otu_taxonomy, Fasta = otu_fasta[,-1]), by="row.names")

###sig used for pchs
natural_endo_glmLRT_otu_tax_fasta$sig <- as.factor(ifelse(natural_endo_glmLRT_otu_tax_fasta$FDR >= 0.05,"NS",
                                                          ifelse(natural_endo_glmLRT_otu_tax_fasta$logFC > 0,"Enriched","Depleted")))

###phylum used for colors :"Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others"

natural_endo_glmLRT_otu_tax_fasta$color <- as.factor(ifelse(natural_endo_glmLRT_otu_tax_fasta$Phylum %in% c("P__Firmicutes", "P__Actinobacteria", "P__Bacteroidetes"), 
                                                            as.character(natural_endo_glmLRT_otu_tax_fasta$Phylum),
                                                            ifelse(natural_endo_glmLRT_otu_tax_fasta$Class %in% c("C__Alphaproteobacteria", "C__Betaproteobacteria", "C__Gammaproteobacteria"),
                                                                   as.character(natural_endo_glmLRT_otu_tax_fasta$Class), "P__Others")))

###the lowest level used for sig asv labels
natural_endo_glmLRT_otu_tax_fasta$label <- as.factor(ifelse(natural_endo_glmLRT_otu_tax_fasta$FDR >= 0.05,"",
                                                            ifelse(natural_endo_glmLRT_otu_tax_fasta$Genus != "", as.character(natural_endo_glmLRT_otu_tax_fasta$Genus),
                                                                   ifelse(natural_endo_glmLRT_otu_tax_fasta$Family != "", as.character(natural_endo_glmLRT_otu_tax_fasta$Family), 
                                                                          ifelse(natural_endo_glmLRT_otu_tax_fasta$Order != "", as.character(natural_endo_glmLRT_otu_tax_fasta$Order),
                                                                                 ifelse(natural_endo_glmLRT_otu_tax_fasta$Class != "", as.character(natural_endo_glmLRT_otu_tax_fasta$Class), 
                                                                                        as.character(natural_endo_glmLRT_otu_tax_fasta$Phylum)))))))

write.csv(natural_endo_glmLRT_otu_tax_fasta, './result/natural_endo_glmLRT_otu_tax_fasta.csv', quote = FALSE)  

dge_de_natural_endo <- decideTestsDGE(lrt_natural_endo, adjust.method = 'fdr', p.value = 0.05)
summary(dge_de_natural_endo)

###Optimized Visualization using ggplot2
windows(16,9)
ggplot(natural_endo_glmLRT_otu_tax_fasta, aes(logCPM, logFC, shape=sig, color=color, fill=color)) +
  geom_point(cex=3+5*(natural_endo_glmLRT_otu_tax_fasta$Col0+natural_endo_glmLRT_otu_tax_fasta$mfec))+
  scale_shape_manual(values=c(25,24,21)) +
  scale_color_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  scale_fill_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="black", linetype="dashed") +
  geom_text(aes(label = label), colour="grey50", vjust = 1.5, cex=3+5*(natural_endo_glmLRT_otu_tax_fasta$Col0+natural_endo_glmLRT_otu_tax_fasta$mfec)) +
  labs(x=paste("Log10CPM"),y=paste("Log10FC")) +
  ggtitle(paste("mfec vs Col-0")) +
  theme_bw() + theme(title = element_text(size = 10, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))


###Syncom_endo_mfec GLM: negative binomial generalized log-linear model
fit_syncom_endo_mfec <- glmFit(dge_syncom_endo_mfec, design_syncom_endo_mfec, robust = TRUE)
lrt_syncom_endo_mfec <- glmLRT(fit_syncom_endo_mfec)

topTags(lrt_syncom_endo_mfec, adjust.method = "fdr", sort.by = "PValue")

###Collect all the informations and save the data
###Relative abundance, GLM stats, taxonomy, fasta
syncom_endo_mfec_glmLRT_otu <- cbind(Col0 = rowMeans(dgelist_syncom_endo_mfec$counts[,dgelist_syncom_endo_mfec$samples$group == "Col0"])/1000, 
                                     mfec = rowMeans(dgelist_syncom_endo_mfec$counts[,dgelist_syncom_endo_mfec$samples$group == "mfec"])/1000,
                                     as.data.frame(topTags(lrt_syncom_endo_mfec, n = dim(dgelist_syncom_endo_mfec$counts)[1], adjust.method = "fdr",sort.by = "none")))

syncom_endo_mfec_glmLRT_otu_tax_fasta <- merge.data.frame(syncom_endo_mfec_glmLRT_otu, cbind(otu_taxonomy, Fasta = otu_fasta[,-1]), by="row.names")

###sig used for pchs
syncom_endo_mfec_glmLRT_otu_tax_fasta$sig <- as.factor(ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$FDR >= 0.05,"NS",
                                                              ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$logFC > 0,"Enriched","Depleted")))

###phylum used for colors :"Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others"

syncom_endo_mfec_glmLRT_otu_tax_fasta$color <- as.factor(ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$Phylum %in% c("P__Firmicutes", "P__Actinobacteria", "P__Bacteroidetes"), 
                                                                as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Phylum),
                                                                ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$Class %in% c("C__Alphaproteobacteria", "C__Betaproteobacteria", "C__Gammaproteobacteria"),
                                                                       as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Class), "P__Others")))

###the lowest level used for sig asv labels
syncom_endo_mfec_glmLRT_otu_tax_fasta$label <- as.factor(ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$FDR >= 0.05,"",
                                                                ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$Genus != "", as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Genus),
                                                                       ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$Family != "", as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Family), 
                                                                              ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$Order != "", as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Order),
                                                                                     ifelse(syncom_endo_mfec_glmLRT_otu_tax_fasta$Class != "", as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Class), 
                                                                                            as.character(syncom_endo_mfec_glmLRT_otu_tax_fasta$Phylum)))))))

write.csv(syncom_endo_mfec_glmLRT_otu_tax_fasta, './result/syncom_endo_mfec_glmLRT_otu_tax_fasta.csv', quote = FALSE)  

dge_de_syncom_endo_mfec <- decideTestsDGE(lrt_syncom_endo_mfec, adjust.method = 'fdr', p.value = 0.05)
summary(dge_de_syncom_endo_mfec)

###Optimized Visualization using ggplot2
windows(16,9)
ggplot(syncom_endo_mfec_glmLRT_otu_tax_fasta, aes(logCPM, logFC, shape=sig, color=color, fill=color)) +
  geom_point(cex=3+5*(syncom_endo_mfec_glmLRT_otu_tax_fasta$Col0+syncom_endo_mfec_glmLRT_otu_tax_fasta$mfec))+
  scale_shape_manual(values=c(25,24,21)) +
  scale_color_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  scale_fill_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="black", linetype="dashed") +
  geom_text(aes(label = label), colour="grey50", vjust = 1.5, cex=3+5*(syncom_endo_mfec_glmLRT_otu_tax_fasta$Col0+syncom_endo_mfec_glmLRT_otu_tax_fasta$mfec)) +
  labs(x=paste("Log10CPM"),y=paste("Log10FC")) +
  ggtitle(paste("mfec vs Col-0")) +
  theme_bw() + theme(title = element_text(size = 10, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))


##########################################################
###Syncom_endo_fec GLM: negative binomial generalized log-linear model
fit_syncom_endo_fec <- glmFit(dge_syncom_endo_fec, design_syncom_endo_fec, robust = TRUE)
lrt_syncom_endo_fec <- glmLRT(fit_syncom_endo_fec)

topTags(lrt_syncom_endo_fec, adjust.method = "fdr", sort.by = "PValue")

###Collect all the informations and save the data
###Relative abundance, GLM stats, taxonomy, fasta
syncom_endo_fec_glmLRT_otu <- cbind(Col0 = rowMeans(dgelist_syncom_endo_fec$counts[,dgelist_syncom_endo_fec$samples$group == "Col0"])/1000, 
                                    fec = rowMeans(dgelist_syncom_endo_fec$counts[,dgelist_syncom_endo_fec$samples$group == "fec"])/1000,
                                    as.data.frame(topTags(lrt_syncom_endo_fec, n = dim(dgelist_syncom_endo_fec$counts)[1], adjust.method = "fdr",sort.by = "none")))

syncom_endo_fec_glmLRT_otu_tax_fasta <- merge.data.frame(syncom_endo_fec_glmLRT_otu, cbind(otu_taxonomy, Fasta = otu_fasta[,-1]), by="row.names")

###sig used for pchs
syncom_endo_fec_glmLRT_otu_tax_fasta$sig <- as.factor(ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$FDR >= 0.05,"NS",
                                                             ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$logFC > 0,"Enriched","Depleted")))

###phylum used for colors :"Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others"

syncom_endo_fec_glmLRT_otu_tax_fasta$color <- as.factor(ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$Phylum %in% c("P__Firmicutes", "P__Actinobacteria", "P__Bacteroidetes"), 
                                                               as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Phylum),
                                                               ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$Class %in% c("C__Alphaproteobacteria", "C__Betaproteobacteria", "C__Gammaproteobacteria"),
                                                                      as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Class), "P__Others")))

###the lowest level used for sig asv labels
syncom_endo_fec_glmLRT_otu_tax_fasta$label <- as.factor(ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$FDR >= 0.05,"",
                                                               ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$Genus != "", as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Genus),
                                                                      ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$Family != "", as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Family), 
                                                                             ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$Order != "", as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Order),
                                                                                    ifelse(syncom_endo_fec_glmLRT_otu_tax_fasta$Class != "", as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Class), 
                                                                                           as.character(syncom_endo_fec_glmLRT_otu_tax_fasta$Phylum)))))))

write.csv(syncom_endo_fec_glmLRT_otu_tax_fasta, './result/syncom_endo_fec_glmLRT_otu_tax_fasta.csv', quote = FALSE)  

dge_de_syncom_endo_fec <- decideTestsDGE(lrt_syncom_endo_fec, adjust.method = 'fdr', p.value = 0.05)
summary(dge_de_syncom_endo_fec)

###Optimized Visualization using ggplot2
windows(16,9)
ggplot(syncom_endo_fec_glmLRT_otu_tax_fasta, aes(logCPM, logFC, shape=sig, color=color, fill=color)) +
  geom_point(cex=3+5*(syncom_endo_fec_glmLRT_otu_tax_fasta$Col0+syncom_endo_fec_glmLRT_otu_tax_fasta$fec))+
  scale_shape_manual(values=c(25,24,21)) +
  scale_color_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  scale_fill_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="black", linetype="dashed") +
  geom_text(aes(label = label), colour="grey50", vjust = 1.5, cex=3+5*(syncom_endo_fec_glmLRT_otu_tax_fasta$Col0+syncom_endo_fec_glmLRT_otu_tax_fasta$fec)) +
  labs(x=paste("Log10CPM"),y=paste("Log10FC")) +
  ggtitle(paste("fec vs Col-0")) +
  theme_bw() + theme(title = element_text(size = 10, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))

######################################################################
###Syncom_endo_min7 GLM: negative binomial generalized log-linear model
fit_syncom_endo_min7 <- glmFit(dge_syncom_endo_min7, design_syncom_endo_min7, robust = TRUE)
lrt_syncom_endo_min7 <- glmLRT(fit_syncom_endo_min7)

topTags(lrt_syncom_endo_min7, adjust.method = "fdr", sort.by = "PValue")

###Collect all the informations and save the data
###Relative abundance, GLM stats, taxonomy, fasta
syncom_endo_min7_glmLRT_otu <- cbind(Col0 = rowMeans(dgelist_syncom_endo_min7$counts[,dgelist_syncom_endo_min7$samples$group == "Col0"])/1000, 
                                     min7 = rowMeans(dgelist_syncom_endo_min7$counts[,dgelist_syncom_endo_min7$samples$group == "min7"])/1000,
                                     as.data.frame(topTags(lrt_syncom_endo_min7, n = dim(dgelist_syncom_endo_min7$counts)[1], adjust.method = "fdr",sort.by = "none")))

syncom_endo_min7_glmLRT_otu_tax_fasta <- merge.data.frame(syncom_endo_min7_glmLRT_otu, cbind(otu_taxonomy, Fasta = otu_fasta[,-1]), by="row.names")

###sig used for pchs
syncom_endo_min7_glmLRT_otu_tax_fasta$sig <- as.factor(ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$FDR >= 0.05,"NS",
                                                              ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$logFC > 0,"Enriched","Depleted")))

###phylum used for colors :"Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others"

syncom_endo_min7_glmLRT_otu_tax_fasta$color <- as.factor(ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$Phylum %in% c("P__Firmicutes", "P__Actinobacteria", "P__Bacteroidetes"), 
                                                                as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Phylum),
                                                                ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$Class %in% c("C__Alphaproteobacteria", "C__Betaproteobacteria", "C__Gammaproteobacteria"),
                                                                       as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Class), "P__Others")))

###the lowest level used for sig asv labels
syncom_endo_min7_glmLRT_otu_tax_fasta$label <- as.factor(ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$FDR >= 0.05,"",
                                                                ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$Genus != "", as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Genus),
                                                                       ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$Family != "", as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Family), 
                                                                              ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$Order != "", as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Order),
                                                                                     ifelse(syncom_endo_min7_glmLRT_otu_tax_fasta$Class != "", as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Class), 
                                                                                            as.character(syncom_endo_min7_glmLRT_otu_tax_fasta$Phylum)))))))

write.csv(syncom_endo_min7_glmLRT_otu_tax_fasta, './result/syncom_endo_min7_glmLRT_otu_tax_fasta.csv', quote = FALSE)  

dge_de_syncom_endo_min7 <- decideTestsDGE(lrt_syncom_endo_min7, adjust.method = 'fdr', p.value = 0.05)
summary(dge_de_syncom_endo_min7)

###Optimized Visualization using ggplot2
windows(16,9)
ggplot(syncom_endo_min7_glmLRT_otu_tax_fasta, aes(logCPM, logFC, shape=sig, color=color, fill=color)) +
  geom_point(cex=3+5*(syncom_endo_min7_glmLRT_otu_tax_fasta$Col0+syncom_endo_min7_glmLRT_otu_tax_fasta$min7))+
  scale_shape_manual(values=c(25,24,21)) +
  scale_color_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  scale_fill_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="black", linetype="dashed") +
  geom_text(aes(label = label), colour="grey50", vjust = 1.5, cex=3+5*(syncom_endo_min7_glmLRT_otu_tax_fasta$Col0+syncom_endo_min7_glmLRT_otu_tax_fasta$min7)) +
  labs(x=paste("Log10CPM"),y=paste("Log10FC")) +
  ggtitle(paste("min7 vs Col-0")) +
  theme_bw() + 
  theme(title = element_text(size = 10, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))




###(3) Optimized Visualization using the EvolView program


##Part2, the cad1 part

###(1) load the data
###design
design <- read.csv(file = "./data/cad1_metadata.csv", header = T)
names(design)

###cad1 syncom endo dataset
cad1_syncom_endo <- subset(design, input_com == "SynCom" & compartment == "endosphere")

###OTU table & taxonomy & fasta
cad1_otu_table <- read.csv(file = "./data/cad1_asv_table.csv", header = T, row.names = 2) ###this is the rarefied OTU/ASV table, rarefied to 3652 reads

cad1_otu_taxonomy <- read.csv(file = "./data/cad1_asv_taxonomy.csv", header = T)
cad1_otu_taxonomy <- subset(cad1_otu_taxonomy, cad1_otu_taxonomy$Feature_ID %in% cad1_otu_table$Feature_ID)
match(cad1_otu_taxonomy$Feature_ID, cad1_otu_table$Feature_ID)
row.names(cad1_otu_taxonomy) <- row.names(cad1_otu_table)

cad1_otu_fasta <- read.csv(file = "./data/cad1_asv_fasta.csv", header = T)
cad1_otu_fasta <- subset(cad1_otu_fasta, cad1_otu_fasta$Feature_ID %in% cad1_otu_table$Feature_ID)
match(cad1_otu_fasta$Feature_ID, cad1_otu_fasta$Feature_ID)
row.names(cad1_otu_fasta) <- row.names(cad1_otu_table)

###subset OTUs/ASVs
otu_table_cad1_syncom_endo <- cad1_otu_table[, names(cad1_otu_table) %in% cad1_syncom_endo$sample_id]

group_cad1_syncom_endo <- c(rep("cad1",20), rep("Col0", 20))
group_cad1_syncom_endo

###(2) Differentially abundant OTUs were detected using EdgeR’s generalized linear model (GLM) approach

###Creates a DGEList object from a table of counts
dgelist_cad1_syncom_endo <- DGEList(counts = otu_table_cad1_syncom_endo, group = group_cad1_syncom_endo)

###Remove low abundance OTUs
###低表达量的otu受到随机因素影响比较大，故其统计结果也不可靠，还会影响p值校正过程
dgelist_cad1_syncom_endo <- dgelist_cad1_syncom_endo[rowSums(dgelist_cad1_syncom_endo$counts) >= 10, ]

###plot MDS
###This results are same to pcoa plots described before
plotMDS(dgelist_cad1_syncom_endo, col = ifelse(cad1_syncom_endo$genotype == "Col0", "red", "darkgreen") , dim = c(1, 2))

###model.matrix
design_cad1_syncom_endo <- model.matrix(~group_cad1_syncom_endo)

###Estimate dispersions
dge_cad1_syncom_endo <- estimateDisp(dgelist_cad1_syncom_endo, design_cad1_syncom_endo, robust = TRUE)

###cad1_syncom_endo GLM: negative binomial generalized log-linear model
fit_cad1_syncom_endo <- glmFit(dge_cad1_syncom_endo, design_cad1_syncom_endo, robust = TRUE)
lrt_cad1_syncom_endo <- glmLRT(fit_cad1_syncom_endo)

topTags(lrt_cad1_syncom_endo, adjust.method = "fdr", sort.by = "PValue")

###Collect all the informations and save the data
###Relative abundance, GLM stats, taxonomy, fasta
cad1_syncom_endo_glmLRT_otu <- cbind(Col0 = rowMeans(dgelist_cad1_syncom_endo$counts[,dgelist_cad1_syncom_endo$samples$group == "Col0"])/3562, 
                                     cad1 = rowMeans(dgelist_cad1_syncom_endo$counts[,dgelist_cad1_syncom_endo$samples$group == "cad1"])/3562,
                                     as.data.frame(topTags(lrt_cad1_syncom_endo, n = dim(dgelist_cad1_syncom_endo$counts)[1], adjust.method = "fdr",sort.by = "none")))

cad1_syncom_endo_glmLRT_otu_tax_fasta <- merge.data.frame(cad1_syncom_endo_glmLRT_otu, cbind(cad1_otu_taxonomy, Fasta = cad1_otu_fasta[,-1]), by="row.names")

###sig used for pchs
cad1_syncom_endo_glmLRT_otu_tax_fasta$sig <- as.factor(ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$FDR >= 0.05,"NS",
                                                              ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$logFC > 0,"Depleted","Enriched")))

###phylum used for colors :"Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others"

cad1_syncom_endo_glmLRT_otu_tax_fasta$color <- as.factor(ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$Phylum %in% c("P__Firmicutes", "P__Actinobacteria", "P__Bacteroidetes"), 
                                                                as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Phylum),
                                                                ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$Class %in% c("C__Alphaproteobacteria", "C__Betaproteobacteria", "C__Gammaproteobacteria"),
                                                                       as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Class), "P__Others")))

###the lowest level used for sig asv labels
cad1_syncom_endo_glmLRT_otu_tax_fasta$label <- as.factor(ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$FDR >= 0.05,"",
                                                                ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$Genus != "", as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Genus),
                                                                       ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$Family != "", as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Family), 
                                                                              ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$Order != "", as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Order),
                                                                                     ifelse(cad1_syncom_endo_glmLRT_otu_tax_fasta$Class != "", as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Class), 
                                                                                            as.character(cad1_syncom_endo_glmLRT_otu_tax_fasta$Phylum)))))))

write.csv(cad1_syncom_endo_glmLRT_otu_tax_fasta, './result/cad1_syncom_endo_glmLRT_otu_tax_fasta.csv', quote = FALSE)  

dge_de_cad1_syncom_endo <- decideTestsDGE(lrt_cad1_syncom_endo, adjust.method = 'fdr', p.value = 0.05)
summary(dge_de_cad1_syncom_endo)

###Optimized Visualization using ggplot2
windows(16,9)
ggplot(cad1_syncom_endo_glmLRT_otu_tax_fasta, aes(logCPM, -logFC, shape=sig, color=color, fill=color)) +
  geom_point(cex=3+5*(cad1_syncom_endo_glmLRT_otu_tax_fasta$Col0+cad1_syncom_endo_glmLRT_otu_tax_fasta$cad1))+
  scale_shape_manual(values=c(25,24,21)) +
  scale_color_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  scale_fill_manual(values=c("red","darkgreen", "green","blue","purple", "cyan", "grey")) +
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed") +
  geom_hline(aes(yintercept=-1), colour="black", linetype="dashed") +
  geom_text(aes(label = label), colour="grey50", vjust = 1.5, cex=3+5*(cad1_syncom_endo_glmLRT_otu_tax_fasta$Col0+cad1_syncom_endo_glmLRT_otu_tax_fasta$cad1)) +
  labs(x=paste("Log10CPM"),y=paste("Log10FC")) +
  ggtitle(paste("cad1 vs Col-0")) +
  theme_bw() + 
  theme(title = element_text(size = 10, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))










