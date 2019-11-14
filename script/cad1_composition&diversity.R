
setwd("F:\\#R\\Others\\XingXF\\R_markdown_for_published")

library(vegan)

###phylum
phylum <- read.csv(file = "./data/cad1_phyla_distrubitions.csv", header = T, row.names = 1)
names(phylum)

#################################
###genotype average barplots
###natural_total_design
###
phylum_natural_total <- subset(phylum, input.com == "AtMix" & compartment == "total")
phylum_natural_endo <- subset(phylum, input.com == "AtMix" & compartment == "endosphere")

phylum_syncom_total <- subset(phylum, input.com == "SynCom" & compartment == "total")
phylum_syncom_endo <- subset(phylum, input.com == "SynCom" & compartment == "endosphere")

table(phylum_syncom_endo$genotype)
table(phylum_natural_total$genotype)
###

phylum_natural_total_total <- as.data.frame(decostand(phylum_natural_total[,1:7],"total"))
rowSums(phylum_natural_total_total)

phylum_natural_endo_total <- as.data.frame(decostand(phylum_natural_endo[,1:7],"total"))

phylum_syncom_total_total <- as.data.frame(decostand(phylum_syncom_total[,1:7],"total"))

phylum_syncom_endo_total <- as.data.frame(decostand(phylum_syncom_endo[,1:7],"total"))

#################################
###
col_phylum <- c("red","darkgreen", "green","blue","purple", "cyan", "grey")


###genotype average barplots
phylum_natural_total_design <- cbind.data.frame(phylum_natural_total_total, phylum_natural_total[,8:12])
phylum_natural_total_design$genotype <- factor(phylum_natural_total_design$genotype,levels=c("Col-0","cad1"))

phylum_natural_average <- as.data.frame(aggregate(phylum_natural_total_design[,1:7], 
                                                  list(phylum_natural_total_design$genotype), 
                                                  mean))
row.names(phylum_natural_average) <- phylum_natural_average$Group.1
phylum_natural_average <- phylum_natural_average[,-1]

##############################

windows()
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_natural_average.bar <- barplot(as.matrix(t(phylum_natural_average)), col = col_phylum,
                                      space = c(1,1), ylab = "Relative Abundance")
legend(x=5, y=1, legend = names(phylum_natural_average), fill = col_phylum, cex = 1)

##############################
############################################################
###natural_endo_design
phylum_natural_endo_design <- cbind.data.frame(phylum_natural_endo_total, phylum_natural_endo)
phylum_natural_endo_design$genotype <- factor(phylum_natural_endo_design$genotype,levels=c("Col-0","cad1"))

phylum_natural_endo_average <- as.data.frame(aggregate(phylum_natural_endo_design[,1:7], 
                                                       list(phylum_natural_endo_design$genotype), 
                                                       mean))

row.names(phylum_natural_endo_average) <- phylum_natural_endo_average$Group.1
phylum_natural_endo_average <- phylum_natural_endo_average[,-1]

###
windows()
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_natural_endo_average.bar <- barplot(as.matrix(t(phylum_natural_endo_average)), col = col_phylum,
                                           space = c(1,1), ylab = "Relative Abundance")

legend(x=5, y=1, legend = names(phylum_natural_endo_average), fill = col_phylum, cex = 1)


###########################################################
###syncom_total_design
phylum_syncom_total_design <- cbind.data.frame(phylum_syncom_total_total, phylum_syncom_total)
phylum_syncom_total_design$genotype <- factor(phylum_syncom_total_design$genotype,levels=c("Col-0","cad1"))

phylum_syncom_total_average <- as.data.frame(aggregate(phylum_syncom_total_design[,1:7], 
                                                       list(phylum_syncom_total_design$genotype), 
                                                       mean))
row.names(phylum_syncom_total_average) <- phylum_syncom_total_average$Group.1
phylum_syncom_total_average <- phylum_syncom_total_average[,-1]

####

windows()
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_syncom_total_average.bar <- barplot(as.matrix(t(phylum_syncom_total_average)), col = col_phylum,
                                           space = c(1,1), ylab = "Relative Abundance")

legend(x=8.5, y=1, legend = names(phylum_syncom_total_average), fill = col_phylum, cex = 1)

###########################################################
###syncom_endo_design
phylum_syncom_endo_design <- cbind.data.frame(phylum_syncom_endo_total, phylum_syncom_endo)
phylum_syncom_endo_design$genotype <- factor(phylum_syncom_endo_design$genotype,levels=c("Col-0","cad1"))

phylum_syncom_endo_average <- as.data.frame(aggregate(phylum_syncom_endo_design[,1:7], 
                                                      list(phylum_syncom_endo_design$genotype), 
                                                      mean))
row.names(phylum_syncom_endo_average) <- phylum_syncom_endo_average$Group.1
phylum_syncom_endo_average <- phylum_syncom_endo_average[,-1]

####

windows()
par(las=2, mar=c(6,4,2,15), xpd=T)
phylum_syncom_endo_average.bar <- barplot(as.matrix(t(phylum_syncom_endo_average)), col = col_phylum,
                                          space = c(1,1), ylab = "Relative Abundance")

legend(x=8.5, y=1, legend = names(phylum_syncom_endo_average), fill = col_phylum, cex = 1)


####################################################################################
####################################################################################
####################################################################################
###alpha diversity

alpha_div <- read.csv(file = "./data/cad1_alpha_diversity.csv", header = T, row.names = 1)
names(alpha_div)
alpha_div$genotype <- factor(alpha_div$genotype,levels=c("Col-0","cad1"))

###all dataset
alpha_natural_total <- subset(alpha_div, input.com == "AtMix" & compartment == "total")
alpha_natural_endo <- subset(alpha_div, input.com == "AtMix" & compartment == "endosphere")

alpha_syncom_total <- subset(alpha_div, input.com == "SynCom" & compartment == "total")
alpha_syncom_endo <- subset(alpha_div, input.com == "SynCom" & compartment == "endosphere")


library(ggplot2)

#################################################################
###alpha-diversity
################################################

windows(3,3)
ggplot(alpha_natural_total, aes(genotype, observed_otus)) + ylim(0,100) +
  geom_boxplot(aes(fill=genotype)) +
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
ggplot(alpha_natural_total, aes(genotype, shannon)) + ylim(0,5) +
  geom_boxplot(aes(fill=genotype)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Shannon index") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


#################################
windows(3,3)
ggplot(alpha_natural_endo, aes(genotype, observed_otus)) + ylim(0,100) +
  geom_boxplot(aes(fill=genotype)) +
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
ggplot(alpha_natural_endo, aes(genotype, shannon)) + ylim(0,5) +
  geom_boxplot(aes(fill=genotype)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Shannon index") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


################################################

windows(3,3)
ggplot(alpha_syncom_total, aes(genotype, observed_otus)) + ylim(0,100) +
  geom_boxplot(aes(fill=genotype)) +
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
ggplot(alpha_syncom_total, aes(genotype, shannon)) + ylim(0,5) +
  geom_boxplot(aes(fill=genotype)) +
  geom_jitter(size=3, bg ="black", position = position_dodge(width = 0.75)) +
  scale_shape_manual(values=c(21:25)) +
  scale_fill_manual(values=c("red", "darkgreen")) + 
  labs(x = "", y = "Shannon index") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle=30, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"))


#################################
windows(3,3)
ggplot(alpha_syncom_endo, aes(genotype, observed_otus)) + ylim(0,100) +
  geom_boxplot(aes(fill=genotype)) +
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
ggplot(alpha_syncom_endo, aes(genotype, shannon)) + ylim(0,5) +
  geom_boxplot(aes(fill=genotype)) +
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
###ÏÔÖøÐÔ¼ìÑé
library(agricolae)
natural_total_sobs.fdr <- kruskal(alpha_natural_total$observed_otus, alpha_natural_total$genotype, group=TRUE, p.adj="fdr")
natural_total_shannon.fdr <- kruskal(alpha_natural_total$shannon, alpha_natural_total$genotype, group=TRUE, p.adj="fdr")

natural_endo_sobs.fdr <- kruskal(alpha_natural_endo$observed_otus, alpha_natural_endo$genotype, group=TRUE, p.adj="fdr")
natural_endo_shannon.fdr <- kruskal(alpha_natural_endo$shannon, alpha_natural_endo$genotype, group=TRUE, p.adj="fdr")

syncom_total_sobs.fdr <- kruskal(alpha_syncom_total$observed_otus, alpha_syncom_total$genotype, group=TRUE, p.adj="fdr")
syncom_total_shannon.fdr <- kruskal(alpha_syncom_total$shannon, alpha_syncom_total$genotype, group=TRUE, p.adj="fdr")

syncom_endo_sobs.fdr <- kruskal(alpha_syncom_endo$observed_otus, alpha_syncom_endo$genotype, group=TRUE, p.adj="fdr")
syncom_endo_shannon.fdr <- kruskal(alpha_syncom_endo$shannon, alpha_syncom_endo$genotype, group=TRUE, p.adj="fdr")


