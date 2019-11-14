
setwd("F:\\#R\\Others\\XingXF\\R_markdown_for_published")

library(vegan)
library(ggplot2)
library(agricolae)

###design
design <- read.csv(file = "./data/mfec_metadata.csv", header = T)
names(design)

###all dataset
###remove samples with reads < 1000
natural_total <- subset(design, Design == "Natural_total" & Bacteria >= 1000)
natural_endo <- subset(design, Design == "Natural_endo" & Bacteria >= 1000)
syncom_endo <- subset(design, Design == "Syncom_endo" & Bacteria >= 1000)

###d5/d6 datasets
natural_total_d5 <- subset(natural_total, Days == "D5" & Experiment != "I1")
natural_endo_d5 <- subset(natural_endo,  Days == "D5")
syncom_endo_d5 <- subset(syncom_endo,  Days == "D6" & Experiment != "I2")


###barplot of bacteria in phylum level
phylum_raw <- read.csv(file = "./data/mfec_phyla_raw.csv", header = T, row.names = 1)
phylum_raw <- phylum_raw[,1:21]
phylum_order <- phylum_raw[, order(colSums(phylum_raw),decreasing = T)]
names(phylum_order)
phylum <- cbind.data.frame(phylum_order[,2:7], phylum_order[,1]-rowSums(phylum_order[,2:4]) + rowSums(phylum_order[,8:21]))
names(phylum) <- c("Betaproteobacteria","Gammaproteobacteria", "Alphaproteobacteria","Firmicutes", "Actinobacteria", "Bacteroidetes", "Others")

###
#################################
###genotype average barplots
###natural_total_design
###
phylum_natural_total_d5 <- subset(phylum, row.names(phylum) %in% natural_total_d5$sample_id)
phylum_natural_endo_d5 <- subset(phylum, row.names(phylum) %in% natural_endo_d5$sample_id)
phylum_syncom_endo_d5 <- subset(phylum, row.names(phylum) %in% syncom_endo_d5$sample_id)

###
phylum_natural_total_d5_total <- as.data.frame(decostand(phylum_natural_total_d5,"total"))
rowSums(phylum_natural_total_d5_total)
phylum_natural_endo_d5_total <- as.data.frame(decostand(phylum_natural_endo_d5,"total"))
phylum_syncom_endo_d5_total <- as.data.frame(decostand(phylum_syncom_endo_d5,"total"))

#################################
###genotype average barplots
###natural_total_d5_design
col_phylum <- c("red","darkgreen", "green","blue","purple", "cyan", "grey")

phylum_natural_total_d5_design <- cbind.data.frame(phylum_natural_total_d5_total, natural_total_d5)
phylum_natural_d5_average <- as.data.frame(aggregate(phylum_natural_total_d5_design[,1:7], 
                                                     list(phylum_natural_total_d5_design$Genotype), 
                                                     mean))

row.names(phylum_natural_d5_average) <- phylum_natural_d5_average$Group.1
phylum_natural_d5_average <- phylum_natural_d5_average[,-1]

################Figure 1D
windows()
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_natural_d5_average.bar <- barplot(as.matrix(t(phylum_natural_d5_average)), col = col_phylum,
                                         space = c(1,1), ylab = "Relative Abundance")

legend(x=5, y=1, legend = names(phylum_natural_d5_average), fill = col_phylum, cex = 1)


###natural_endo_d5_design
phylum_natural_endo_d5_design <- cbind.data.frame(phylum_natural_endo_d5_total, natural_endo_d5)
phylum_natural_endo_d5_average <- as.data.frame(aggregate(phylum_natural_endo_d5_design[,1:7], 
                                                          list(phylum_natural_endo_d5_design$Genotype), 
                                                          mean))
row.names(phylum_natural_endo_d5_average) <- phylum_natural_endo_d5_average$Group.1
phylum_natural_endo_d5_average <- phylum_natural_endo_d5_average[,-1]


#######################Figure 1D
windows()
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_natural_endo_d5_average.bar <- barplot(as.matrix(t(phylum_natural_endo_d5_average)), col = col_phylum,
                                              space = c(1,1), ylab = "Relative Abundance")
legend(x=5, y=1, legend = names(phylum_natural_endo_d5_average), fill = col_phylum, cex = 1)



#######################Figure 2B
###syncom_endo_d5_design
phylum_syncom_endo_d5_design <- cbind.data.frame(phylum_syncom_endo_d5_total, syncom_endo_d5)
phylum_syncom_endo_d5_design$Genotype <- factor(phylum_syncom_endo_d5_design$Genotype,levels=c("Col0","fec","min7","mfec"))
phylum_syncom_endo_d5_average <- as.data.frame(aggregate(phylum_syncom_endo_d5_design[,1:7], 
                                                         list(phylum_syncom_endo_d5_design$Genotype), 
                                                         mean))
row.names(phylum_syncom_endo_d5_average) <- phylum_syncom_endo_d5_average$Group.1
phylum_syncom_endo_d5_average <- phylum_syncom_endo_d5_average[,-1]

####

windows(9,6)
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_syncom_endo_d5_average.bar <- barplot(as.matrix(t(phylum_syncom_endo_d5_average)), col = col_phylum,
                                             space = c(1,1), ylab = "Relative Abundance")

legend(x=8.5, y=1, legend = names(phylum_syncom_endo_d5_average), fill = col_phylum, cex = 1)


####################################################################################
###alpha diversity

alpha_div <- read.csv(file = "./data/mfec_alpha_diversity.csv", header = T)
names(alpha_div)

###all dataset
alpha_natural_total <- merge(alpha_div, natural_total, by="sample_id", all.y = T)
alpha_natural_endo <- merge(alpha_div, natural_endo, by="sample_id", all.y = T)
alpha_syncom_endo <- merge(alpha_div, syncom_endo, by="sample_id", all.y = T)

###d5/d6 datasets
alpha_natural_total_d5 <- subset(alpha_natural_total, Days == "D5" & Experiment != "I1" )
alpha_natural_endo_d5 <- subset(alpha_natural_endo,  Days == "D5")
alpha_syncom_endo_d5 <- subset(alpha_syncom_endo,  Days == "D6" & Experiment != "I2")


table(alpha_natural_total_d5$Genotype)
table(alpha_natural_endo_d5$Genotype)
table(alpha_syncom_endo_d5$Genotype)

#################################################################
###natural_total
alpha_natural_total_d5$Group <- substr(alpha_natural_total_d5$Group,10,20)

################################################
windows(3,3)
ggplot(alpha_natural_total_d5, aes(Group, observed_otus)) + ylim(0,100) +
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Observed OTUs") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


windows(3,3)
ggplot(alpha_natural_total_d5, aes(Group, shannon)) + ylim(0,5) +
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Shannon index") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


#########################################################
###显著性检验

sobs.fdr <- kruskal(alpha_natural_total_d5$observed_otus, alpha_natural_total_d5$Group, group=TRUE, p.adj="fdr")
sobs.fdr

shannon.fdr <- kruskal(alpha_natural_total_d5$shannon, alpha_natural_total_d5$Group, group=TRUE, p.adj="fdr")
shannon.fdr

#################################################################
###natural_endo
alpha_natural_endo_d5$Group <- substr(alpha_natural_endo_d5$Group,7,20)
###

windows(3,3)
ggplot(alpha_natural_endo_d5, aes(Group, observed_otus)) + ylim(0,100) +
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Observed OTUs") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


windows(3,3)
ggplot(alpha_natural_endo_d5, aes(Group, shannon)) + ylim(0,5) +
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Shannon index") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


#########################################################
###显著性检验
sobs.fdr <- kruskal(alpha_natural_endo_d5$observed_otus, alpha_natural_endo_d5$Group, group=TRUE, p.adj="fdr")
sobs.fdr

shannon.fdr <- kruskal(alpha_natural_endo_d5$shannon, alpha_natural_endo_d5$Group, group=TRUE, p.adj="fdr")
shannon.fdr


#################################################################
#################################################################
###alpha_syncom_endo
###
alpha_syncom_endo_d5$Group <- substr(alpha_syncom_endo_d5$Group,14,20)
alpha_syncom_endo_d5$Group <- factor(alpha_syncom_endo_d5$Group,levels=c("Col0","fec","min7","mfec"))
###

windows(4,3)
ggplot(alpha_syncom_endo_d5, aes(Group, observed_otus)) + ylim(0,100) + 
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "green", "cyan", "darkgreen")) + 
  labs(x = "", y = "Observed OTUs") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


windows(4,3)
ggplot(alpha_syncom_endo_d5, aes(Group, shannon)) + ylim(0,5) + 
  geom_boxplot(aes(fill=Group)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "green", "cyan", "darkgreen")) + 
  labs(x = "", y = "Shannon index") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


#########################################################
###显著性检验
syncom_endo_sobs.fdr <- kruskal(alpha_syncom_endo_d5$observed_otus, alpha_syncom_endo_d5$Group, group=TRUE, p.adj="fdr")
syncom_endo_sobs.fdr

syncom_endo_shannon.fdr <- kruskal(alpha_syncom_endo_d5$shannon, alpha_syncom_endo_d5$Group, group=TRUE, p.adj="fdr")
syncom_endo_shannon.fdr

##############################################################################
##############################################################################
##############################################################################
###beta diversity
###bray curtis

bray_curtis <- read.csv(file = "./data/mfec_bray_curtis.csv", header = T, row.names = 1)
names(bray_curtis) <- row.names(bray_curtis)

###all dataset
bray_natural_total <- bray_curtis[match(natural_total$sample_id, row.names(bray_curtis)), match(natural_total$sample_id, names(bray_curtis))]
bray_natural_endo <- bray_curtis[match(natural_endo$sample_id, row.names(bray_curtis)), match(natural_endo$sample_id, names(bray_curtis))]
bray_syncom_endo <- bray_curtis[match(syncom_endo$sample_id, row.names(bray_curtis)), match(syncom_endo$sample_id, names(bray_curtis))]

###d5/d6 datasets
bray_natural_total_d5 <- bray_curtis[match(natural_total_d5$sample_id, row.names(bray_curtis)), match(natural_total_d5$sample_id, names(bray_curtis))]
bray_natural_endo_d5 <- bray_curtis[match(natural_endo_d5$sample_id, row.names(bray_curtis)), match(natural_endo_d5$sample_id, names(bray_curtis))]
bray_syncom_endo_d5 <- bray_curtis[match(syncom_endo_d5$sample_id, row.names(bray_curtis)), match(syncom_endo_d5$sample_id, names(bray_curtis))]

#############################################################
######################################################
###bray_natural_total_d5
bray_natural_total_d5.dist <- as.dist(bray_natural_total_d5)

bray_natural_total_d5_pcoa <- cmdscale(bray_natural_total_d5,k=(nrow(bray_natural_total_d5)-1),eig=T)
bray_natural_total_d5_eig <- bray_natural_total_d5_pcoa$eig
bray_natural_total_d5_pcoa123 <- cbind(as.data.frame(scores(bray_natural_total_d5_pcoa)[,c(1,2,3)]), natural_total_d5)
names(bray_natural_total_d5_pcoa123)[1:3] <- c("PCoA1", "PCoA2", "PCoA3")

bray_natural_total_d5_pcoa123$Group <- substr(bray_natural_total_d5_pcoa123$Group,10,20)

windows(4,3)
ggplot(bray_natural_total_d5_pcoa123, aes(PCoA1, PCoA2, colour=Group, shape=Experiment)) +
  geom_point(aes(fill=Group), size=3) +
  scale_colour_manual(values=c("red", "darkgreen")) +
  scale_fill_manual(values=c("red", "darkgreen")) +
  scale_shape_manual(values=c(21:25)) +
  labs(x=paste("PCoA 1 (", format(100 * bray_natural_total_d5_eig[1] / sum(bray_natural_total_d5_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bray_natural_total_d5_eig[2] / sum(bray_natural_total_d5_eig), digits=4), "%)", sep=""))  +
  theme_bw() + 
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text = element_text(size = 10, colour = "black"))

###############


##########################################################################
###bray_natural_endo_d5
bray_natural_endo_d5.dist <- as.dist(bray_natural_endo_d5)

bray_natural_endo_d5_pcoa <- cmdscale(bray_natural_endo_d5,k=(nrow(bray_natural_endo_d5)-1),eig=T)
bray_natural_endo_d5_eig <- bray_natural_endo_d5_pcoa$eig
bray_natural_endo_d5_pcoa123 <- cbind(as.data.frame(scores(bray_natural_endo_d5_pcoa)[,c(1,2,3)]), natural_endo_d5)
names(bray_natural_endo_d5_pcoa123)[1:3] <- c("PCoA1", "PCoA2", "PCoA3")

bray_natural_endo_d5_pcoa123$Group <- substr(bray_natural_endo_d5_pcoa123$Group,7,20)

windows(4,3)
ggplot(bray_natural_endo_d5_pcoa123, aes(PCoA1, PCoA2, colour=Group, shape=Experiment)) +
  geom_point(aes(fill=Group), size=3) +
  scale_colour_manual(values=c("red", "darkgreen")) +
  scale_fill_manual(values=c("red", "darkgreen")) +
  scale_shape_manual(values=c(21:25)) +
  labs(x=paste("PCoA 1 (", format(100 * bray_natural_endo_d5_eig[1] / sum(bray_natural_endo_d5_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bray_natural_endo_d5_eig[2] / sum(bray_natural_endo_d5_eig), digits=4), "%)", sep=""))  +
  theme_bw() + 
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text = element_text(size = 10, colour = "black"))

###############

############################################
############################################
###bray_syncom_endo_d5
bray_syncom_endo_d5.dist <- as.dist(bray_syncom_endo_d5)

bray_syncom_endo_d5_pcoa <- cmdscale(bray_syncom_endo_d5,k=(nrow(bray_syncom_endo_d5)-1),eig=T)
bray_syncom_endo_d5_eig <- bray_syncom_endo_d5_pcoa$eig
bray_syncom_endo_d5_pcoa123 <- cbind(as.data.frame(scores(bray_syncom_endo_d5_pcoa)[,c(1,2,3)]), syncom_endo_d5)
names(bray_syncom_endo_d5_pcoa123)[1:3] <- c("PCoA1", "PCoA2", "PCoA3")

bray_syncom_endo_d5_pcoa123$Group <- substr(bray_syncom_endo_d5_pcoa123$Group,14,20)
bray_syncom_endo_d5_pcoa123$Group <- factor(bray_syncom_endo_d5_pcoa123$Group,levels=c("Col0","fec","min7","mfec"))

windows(4,3)
ggplot(bray_syncom_endo_d5_pcoa123, aes(PCoA1, PCoA2, colour=Group, shape=Experiment)) +
  geom_point(aes(fill=Group), size=3) +
  scale_colour_manual(values=c("red", "green", "cyan", "darkgreen")) +
  scale_fill_manual(values=c("red", "green","cyan", "darkgreen")) +
  scale_shape_manual(values=c(21:25)) +
  labs(x=paste("PCoA 1 (", format(100 * bray_syncom_endo_d5_eig[1] / sum(bray_syncom_endo_d5_eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * bray_syncom_endo_d5_eig[2] / sum(bray_syncom_endo_d5_eig), digits=4), "%)", sep=""))  +
  theme_bw() + 
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        axis.text = element_text(size = 10, colour = "black"))

###############

###################################################################

####anosim
bray_natural_total_anosim <- anosim(bray_natural_total_d5.dist, bray_natural_total_d5_pcoa123$Group,permutations = 999)
summary(bray_natural_total_anosim)

bray_natural_endo_anosim <- anosim(bray_natural_endo_d5.dist, bray_natural_endo_d5_pcoa123$Group,permutations = 999)
summary(bray_natural_endo_anosim)

bray_syncom_endo_anosim <- anosim(bray_syncom_endo_d5.dist, bray_syncom_endo_d5_pcoa123$Group,permutations = 999)
summary(bray_syncom_endo_anosim)

####adonis2
adonis2(bray_natural_total_d5 ~ Experiment, data = natural_total_d5, add=T, sqrt.dist=T, permutations = 999)
adonis2(bray_natural_total_d5 ~ Genotype, data = natural_total_d5, add=T, sqrt.dist=T, permutations = 999)

adonis2(bray_natural_endo_d5 ~ Experiment, data = natural_endo_d5, add=T, sqrt.dist=T, permutations = 999)
adonis2(bray_natural_endo_d5 ~ Genotype, data = natural_endo_d5, add=T, sqrt.dist=T, permutations = 999)

adonis2(bray_syncom_endo_d5 ~ Experiment, data = syncom_endo_d5, add=T, sqrt.dist=T, permutations = 999)
adonis2(bray_syncom_endo_d5 ~ Genotype, data = syncom_endo_d5, add=T, sqrt.dist=T, permutations = 999)


##############################################################################
##############################################################################

######################################
###lefse分析
lefse <- read.csv(file = "./data/mfec_lefse.csv", header = T, row.names = 1)
matrch_id <- match(gsub("-", ".", design$sample_id), names(lefse))

lefse_design <- lefse[,matrch_id]
match(names(lefse_design), gsub("-", ".", design$sample_id))
names(lefse_design) <- design$sample_id

lefse_design_group <- cbind.data.frame(design, t(lefse_design))

lefse_natural_total <- subset(lefse_design_group, lefse_design_group$sample_id %in% natural_total$sample_id)
lefse_natural_total_d5 <- subset(lefse_design_group, lefse_design_group$sample_id %in% natural_total_d5$sample_id)

lefse_natural_endo <- subset(lefse_design_group, lefse_design_group$sample_id %in% natural_endo$sample_id)
lefse_natural_endo_d5 <- subset(lefse_design_group, lefse_design_group$sample_id %in% natural_endo_d5$sample_id)

lefse_syncom_endo <- subset(lefse_design_group, lefse_design_group$sample_id %in% syncom_endo$sample_id)
lefse_syncom_endo_d5 <- subset(lefse_design_group, lefse_design_group$sample_id %in% syncom_endo_d5$sample_id)

write.csv(lefse_natural_total, file = "./result/lefse_natural_total.csv")
write.csv(lefse_natural_total_d5, file = "./result/lefse_natural_total_d5.csv")

write.csv(lefse_natural_endo, file = "./result/lefse_natural_endo.csv")
write.csv(lefse_natural_endo_d5, file = "./result/lefse_natural_endo_d5.csv")

write.csv(lefse_syncom_endo, file = "./result/lefse_syncom_endo.csv")
write.csv(lefse_syncom_endo_d5, file = "./result/lefse_syncom_endo_d5.csv")











