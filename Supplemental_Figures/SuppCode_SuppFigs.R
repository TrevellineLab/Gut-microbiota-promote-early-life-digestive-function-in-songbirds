################ 
#SUPPLEMENTAL CODE FOR ANALYSES - SUPPLEMENTAL FIGURES AND SUPPORTING ANALYSES #################################
#Gut microbiota promote early-life digestive function in songbirds
#Brian K. Trevelline
#23 MAY 2026
################

#Load packages
library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)
library(Matrix)
library(afex)
library(rcompanion)
library(patchwork)
library(ggdist)
library(gghalves)
library(emmeans)
library(Cairo)
library(outliers)
library(rmarkdown)
library(knitr)


### SUPPLEMENTAL FIGURE S1 - RPCA ORDINATION OF ALL SAMPLE TYPES AND TIMEPOINTS
rpca.df <- read.csv("RPCAordinationFIGS1.csv")

fig.s1.plot <- 
  ggplot(rpca.df, aes(Axis1, Axis2, shape = Treatment, fill = Treatment, col = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_shape_manual(values = c(21,22,23,24,21,22,23,24,24,21,21,21), limits=c("AX-D1-2", "AX-D3", "AX-D5-6", "AX-SICONT", "CV-D1-2","CV-D3", "CV-D5-6", "CV-SICONT", "PR-SICONT", "INOC", "MEC", "BLANK")) +
  scale_color_manual(values = c("#009E73","#009E73", "#009E73","#009E73", "#D55E00", "#D55E00", "#D55E00", "#D55E00", "#CC79A7","#0072B2","#F0E442", "#000000"), limits=c("AX-D1-2", "AX-D3", "AX-D5-6", "AX-SICONT", "CV-D1-2","CV-D3", "CV-D5-6", "CV-SICONT", "PR-SICONT", "INOC", "MEC", "BLANK")) +
  scale_fill_manual(values = c("#009E73","#009E73", "#009E73","#009E73", "#D55E00","#D55E00", "#D55E00", "#D55E00",  "#CC79A7", "#0072B2","#F0E442", "#000000" ), limits=c("AX-D1-2", "AX-D3", "AX-D5-6", "AX-SICONT", "CV-D1-2","CV-D3", "CV-D5-6", "CV-SICONT", "PR-SICONT", "INOC", "MEC", "BLANK")) +
  xlab("RPCA 1 (59.8%)") +
  ylab("RPCA 2 (40.2%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = 'right')

ggsave(filename = "fig.s1.plot.pdf", 
       plot = fig.s1.plot,
       device = cairo_pdf,
       width = 15, 
       height = 12,
       units = "cm",
       dpi = 600)

## FIG S2A - AX/CV RPCA DISTANCE COMPARISION TO BLANKS
rpca.dist <- read.csv("RPCAdistancesFIGS2.csv")

rpca.dist.blank <-
  rpca.dist %>%
  filter(Group1Group2 == "BLANKAX" |
           Group1Group2 == "BLANKCV")

pairwisePermutationTest(
  formula = Distance ~ Group1Group2,
  data = rpca.dist.blank,
  method = "fdr")

fig.s2.a.plot <- 
  ggplot(rpca.dist.blank, aes(Group1Group2, Distance, fill = Group1Group2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("BLANKAX", "BLANKCV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,4.25), breaks = c(0,1,2,3,4)) +
  theme_bw() +
  annotate("text", label = "***", x = 1.5, y = 4.2, size = 5) +
  ylab("Pairwise RPCA distances to BLANK") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s2.a.plot.pdf", 
       plot = fig.s2.a.plot,
       device = cairo_pdf,
       width = 10, 
       height = 10,
       units = "cm",
       dpi = 600)

## RPCA DISTANCE COMPARISION 
## FIG S2B - AX/CV RPCA DISTANCE COMPARISION TO INOC
rpca.dist.inoc <-
  rpca.dist %>%
  filter(Group1Group2 == "INOCAX" |
           Group1Group2 == "INOCCV")

pairwisePermutationTest(
  formula = Distance ~ Group1Group2,
  data = rpca.dist.inoc,
  method = "fdr")

fig.s2.b.plot <- 
  ggplot(rpca.dist.inoc, aes(Group1Group2, Distance, fill = Group1Group2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("INOCAX", "INOCCV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,4.25), breaks = c(0,1,2,3,4)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 4.2, yend = 4.2) +
  annotate("text", label = "***", x = 1.5, y = 4.2, size = 5) +
  ylab("Pairwise RPCA distances to INOC") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s2.b.plot.pdf", 
       plot = fig.s2.b.plot,
       device = cairo_pdf,
       width = 10, 
       height = 10,
       units = "cm",
       dpi = 600)

## FIG S3 - CV RECIPIENT VS NON-RECIPIENT RPCA DISTANCE COMPARISION
rpca.dist.inoc.conv.paired <- read.csv("RPCAdistancesFIGS3.csv")

#TEST EFFECT OF GROUP ON DISTANCE CONTROLLING FOR RANDOM EFFECTS OF INDIVIDUAL
recipient.test <- lmer(Distance ~ InoculaRecipient + (1|SubjectID1), data = rpca.dist.inoc.conv.paired)
summary(recipient.test)

fig.s3.plot <- 
  ggplot(rpca.dist.inoc.conv.paired, aes(InoculaRecipient, Distance, fill = InoculaRecipient)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("#999999", "#999999")) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
  theme_bw() +
  annotate("text", label = "*", x = 1.5, y = 4, size = 5) +
  ylab("Pairwise RPCA distances to INOC") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s3.plot.pdf", 
       plot = fig.s3.plot,
       device = cairo_pdf,
       width = 10, 
       height = 10,
       units = "cm",
       dpi = 600)

## FIG S4A - COMPARISON OF MICROBIOTA COMMUNITY COMPOSITION USING GENUS-LEVEL MICROBIAL ABUNDANCES
rel.abund <- read.csv("ASVRelativeAbundances_Genus.csv")

fig.s4.a.plot <- 
  ggplot(rel.abund, aes(y = RelativeAbundance, x = PlottingTreatment, fill = PlottingGenus)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_x_discrete(labels = c("MEC","AX","CV","INOC","PR")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = c("#009E73","#D55E00","#F0E442","#0072B2","#CC7987","#E69F00" ,"#56B4E9","dark gray"), limits = c( "A", "B","C","D","E","F","G","H"), labels = c("Catellicoccus","Enterococcus","Klebsiella","Lactobacillus","Lactococcus","Serratia","Staphylococcus", "Other"), name = "Bacterial genus") +
  theme_bw() +
  ylab("Relative Abundance") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.justification = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"))

ggsave(filename = "fig.s4.a.plot.pdf", 
       plot = fig.s4.a.plot,
       device = cairo_pdf,
       width = 16, 
       height = 14,
       units = "cm",
       dpi = 600)

## FIG S4B - PERCENT INOC ASVS SUCCESFULLY TRANSFERRED TO CV NESTLINGS
### MATCH ASVS IN INOC AND CONV PAIRS
asv.df <- read.csv("INOC-CV_ASVs.csv")
rownames(asv.df) <- asv.df$FeatureID

#INOCULATIONS
inoc.2301.total.asvs <- count(asv.df, INOC2301 >= 1)
inoc.2302.total.asvs <- count(asv.df, INOC2302 >= 1)
inoc.2304.total.asvs <- count(asv.df, INOC2304 >= 1)
inoc.2305.total.asvs <- count(asv.df, INOC2305 >= 1)
inoc.2306.total.asvs <- count(asv.df, INOC2306 >= 1)
inoc.2311.total.asvs <- count(asv.df, INOC2311 >= 1)
inoc.2317.total.asvs <- count(asv.df, INOC2317 >= 1)
inoc.2318.total.asvs <- count(asv.df, INOC2318 >= 1)

#INOC2301
#NL093 = 51.3% TRANSFERRED
non.trans.asvs.NL093 <- filter(asv.df, INOC2301 >= 1  & NL093 == 0)
trans.asvs.NL093 <- filter(asv.df, INOC2301 >= 1  & NL093 >= 1)
trans.eff.NL093 <- nrow(trans.asvs.NL093)/inoc.2301.total.asvs[2,2]

#INOC2301
#NL094 = 30.7%
non.trans.asvs.NL094 <- filter(asv.df, INOC2301 >= 1  & NL094 == 0)
trans.asvs.NL094 <- filter(asv.df, INOC2301 >= 1  & NL094 >= 1)
trans.eff.NL094 <- nrow(trans.asvs.NL094)/inoc.2301.total.asvs[2,2]

#INOC2311
#NL107 = 28.8%
non.trans.asvs.NL107 <- filter(asv.df, INOC2311 >= 1  & NL107 == 0)
trans.asvs.NL107 <- filter(asv.df, INOC2311 >= 1  & NL107 >= 1)
trans.eff.NL107 <- nrow(trans.asvs.NL107)/inoc.2311.total.asvs[2,2]

#INOC2317
#NL118 = 33.3%
non.trans.asvs.NL118 <- filter(asv.df, INOC2317 >= 1  & NL118 == 0)
trans.asvs.NL118 <- filter(asv.df, INOC2317 >= 1  & NL118 >= 1)
trans.eff.NL118 <- nrow(trans.asvs.NL118)/inoc.2317.total.asvs[2,2]

#INOC2302
#NL123 = 34.2%
non.trans.asvs.NL123 <- filter(asv.df, INOC2302 >= 1  & NL123 == 0)
trans.asvs.NL123 <- filter(asv.df, INOC2302 >= 1  & NL123 >= 1)
trans.eff.NL123 <- nrow(trans.asvs.NL123)/inoc.2302.total.asvs[2,2]

#INOC2302
#NL124 = 50.0%
non.trans.asvs.NL124 <- filter(asv.df, INOC2302 >= 1  & NL124 == 0)
trans.asvs.NL124 <- filter(asv.df, INOC2302 >= 1  & NL124 >= 1)
trans.eff.NL124 <- nrow(trans.asvs.NL124)/inoc.2302.total.asvs[2,2]

#INOC2305
#NL141 = 81.8%
non.trans.asvs.NL141 <- filter(asv.df, INOC2305 >= 1  & NL141 == 0)
trans.asvs.NL141 <- filter(asv.df, INOC2305 >= 1  & NL141 >= 1)
trans.eff.NL141 <- nrow(trans.asvs.NL141)/inoc.2305.total.asvs[2,2]

#INOC2305
#NL143 = 59.1%
non.trans.asvs.NL143 <- filter(asv.df, INOC2305 >= 1  & NL143 == 0)
trans.asvs.NL143 <- filter(asv.df, INOC2305 >= 1  & NL143 >= 1)
trans.eff.NL143 <- nrow(trans.asvs.NL143)/inoc.2305.total.asvs[2,2]

#INOC2305
#NL146 = 59.1%
non.trans.asvs.NL146 <- filter(asv.df, INOC2305 >= 1  & NL146 == 0)
trans.asvs.NL146 <- filter(asv.df, INOC2305 >= 1  & NL146 >= 1)
trans.eff.NL146 <- nrow(trans.asvs.NL146)/inoc.2305.total.asvs[2,2]

#INOC2306
#NL151 = 77.8%
non.trans.asvs.NL151 <- filter(asv.df, INOC2306 >= 1  & NL151 == 0)
trans.asvs.NL151 <- filter(asv.df, INOC2306 >= 1  & NL151 >= 1)
trans.eff.NL151 <- nrow(trans.asvs.NL151)/inoc.2306.total.asvs[2,2]

#INOC2306
#NL152 = 77.8%
non.trans.asvs.NL152 <- filter(asv.df, INOC2306 >= 1  & NL152 == 0)
trans.asvs.NL152 <- filter(asv.df, INOC2306 >= 1  & NL152 >= 1)
trans.eff.NL152 <- nrow(trans.asvs.NL152)/inoc.2306.total.asvs[2,2]

#INOC2318
#NL158 = 68.8%
non.trans.asvs.NL158 <- filter(asv.df, INOC2318 >= 1  & NL158 == 0)
trans.asvs.NL158 <- filter(asv.df, INOC2318 >= 1  & NL158 >= 1)
trans.eff.NL158 <- nrow(trans.asvs.NL158)/inoc.2318.total.asvs[2,2]

#INOC2304
#NL165 = 62.5%
non.trans.asvs.NL165 <- filter(asv.df, INOC2304 >= 1  & NL165 == 0)
trans.asvs.NL165 <- filter(asv.df, INOC2304 >= 1  & NL165 >= 1)
trans.eff.NL165 <- nrow(trans.asvs.NL165)/inoc.2304.total.asvs[2,2]

trans.eff <- c(trans.eff.NL093, 
               trans.eff.NL094,
               trans.eff.NL107,
               trans.eff.NL118,
               trans.eff.NL123,
               trans.eff.NL124,
               trans.eff.NL141,
               trans.eff.NL143,
               trans.eff.NL146,
               trans.eff.NL151,
               trans.eff.NL152,
               trans.eff.NL158,
               trans.eff.NL165)

trans.eff.df <- data.frame(
  AnimalID = c("NL093", 
               "NL094",
               "NL107",
               "NL118",
               "NL123",
               "NL124",
               "NL141",
               "NL143",
               "NL146",
               "NL151",
               "NL152",
               "NL158",
               "NL165"),
  TransferEfficiency = trans.eff
)

#EXTRACT TAXONOMY FOR EACH MATCHING ASV BY NESTLING
trans.asvs <- c(rownames(trans.asvs.NL093), 
                rownames(trans.asvs.NL094),
                rownames(trans.asvs.NL107),
                rownames(trans.asvs.NL118),
                rownames(trans.asvs.NL123),
                rownames(trans.asvs.NL124),
                rownames(trans.asvs.NL141),
                rownames(trans.asvs.NL143),
                rownames(trans.asvs.NL146),
                rownames(trans.asvs.NL151),
                rownames(trans.asvs.NL152),
                rownames(trans.asvs.NL158),
                rownames(trans.asvs.NL165))

trans.asvs.unique <- unique(trans.asvs)
#filter(asv.df, FeatureID == trans.asvs.unique)
transferred.asvs <- merge(data.frame(FeatureID=trans.asvs.unique), asv.df, all.x=TRUE)
write.csv(transferred.asvs, "TransferredASVs.csv")

avg.trans.eff <- mean(trans.eff)
range(trans.eff)
sd(trans.eff)

#NON-TRANSFERRED ASVS
non.trans.asvs <- c(rownames(non.trans.asvs.NL093), 
                rownames(non.trans.asvs.NL094),
                rownames(non.trans.asvs.NL107),
                rownames(non.trans.asvs.NL118),
                rownames(non.trans.asvs.NL123),
                rownames(non.trans.asvs.NL124),
                rownames(non.trans.asvs.NL141),
                rownames(non.trans.asvs.NL143),
                rownames(non.trans.asvs.NL146),
                rownames(non.trans.asvs.NL151),
                rownames(non.trans.asvs.NL152),
                rownames(non.trans.asvs.NL158),
                rownames(non.trans.asvs.NL165))

non.trans.asvs.unique <- unique(non.trans.asvs)
#filter(asv.df, FeatureID == non.trans.asvs.unique)
non.transferred.asvs <- merge(data.frame(FeatureID=non.trans.asvs.unique), asv.df, all.x=TRUE)
write.csv(non.transferred.asvs, "NonTransferredASVs.csv")

fig.s4.b.plot <- 
  ggplot(trans.eff.df, aes(1, TransferEfficiency)) +
  geom_boxplot(outlier.shape = NA, fill = "#999999") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_x_discrete(labels = c("X")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.50,0.75,1)) +
  theme_bw() +
  ylab("Proportion of INOC ASVs transfered to CV") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.position = 'none')

ggsave(filename = "fig.s4.b.plot.pdf", 
       plot = fig.s4.b.plot,
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

## FIG S6 - Mean Villus Width
#IMPORT HISTOLOGY DATA
hist.df <- read.csv("VilliMeasurements.csv")

#VILLI WIDTH
width.test <- lmer(MeanVillusWidth ~ Treatment + Sex + (1|Parentage) + (1|IndividualID), data = hist.df)
summary(width.test)

med.vil.wid <- aggregate(MeanVillusWidth ~ IndividualID + Treatment, data = hist.df, FUN = median, na.rm = TRUE)

fig.s6.plot <- 
  ggplot(hist.df, aes(Treatment, MeanVillusWidth, fill = Treatment)) +
  ggdist::stat_halfeye(width = 0.7) + 
  geom_half_point(side = "l", alpha = 0.3) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,225)) +
  annotate("segment", x = 1.0, xend = 1.3, y = 67.096, yend = 67.096) +
  annotate("segment", x = 1.0, xend = 1.3, y = 59.944, yend = 59.944) +
  annotate("segment", x = 1.0, xend = 1.3, y = 80.888, yend = 80.888) +
  annotate("segment", x = 1.0, xend = 1.3, y = 95.378, yend = 95.378) +
  annotate("segment", x = 1.0, xend = 1.3, y = 58.418, yend = 58.418) +
  annotate("segment", x = 1.0, xend = 1.3, y = 78.484, yend = 78.484) +
  annotate("segment", x = 2.0, xend = 2.3, y = 92.835, yend = 92.835) +
  annotate("segment", x = 2.0, xend = 2.3, y = 77.861, yend = 77.861) +
  annotate("segment", x = 2.0, xend = 2.3, y = 76.038, yend = 76.038) +
  annotate("segment", x = 2.0, xend = 2.3, y = 79.202, yend = 79.202) +
  annotate("segment", x = 2.0, xend = 2.3, y = 95.664, yend = 95.664) +
  annotate("segment", x = 2.0, xend = 2.3, y = 95.970, yend = 95.970) +
  annotate("text", label = "NS", x = 1.5, y = 225, size = 3.4) +
  theme_bw() +
  ylab(expression(paste('Mean villi width (', "\u03bc", 'm)'))) +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s6.plot.pdf", 
       plot = fig.s6.plot,
       device = cairo_pdf,
       width = 5, 
       height = 7.5,
       units = "cm",
       dpi = 600)

### FIG S7 A - Gizzard Mass
morph.df <- read.csv("MorphologyData.csv")

gizzard.test <- lmer(GizzardMassMg ~ Treatment * Sex + TarsusLength + (1|Parentage), data = morph.df)
summary(gizzard.test)

fig.s7.a.plot <- 
  ggplot(morph.df, aes(x = Treatment, y = GizzardMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(200,600), breaks = c(200, 400, 600)) +
  theme_bw() +
  annotate("text", label = "NS", x = 1.5, y = 600, size = 4) +
  ylab("Gizzard mass (mg)") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s7.a.plot.pdf", 
       plot = fig.s7.a.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

#FIG S7B - Heart
heart.test <- lmer(HeartMassMg ~ Treatment + Sex + TarsusLength + (1|Parentage), data = morph.df)
summary(heart.test)

fig.s7.b.plot <- 
  ggplot(morph.df, aes(x = Treatment, y = HeartMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,110), breaks = c(0,50,100)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 105, yend = 105) +
  annotate("text", label = "NS", x = 1.5, y = 110, size = 4) +
  ylab("Heart mass (mg)") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s7.b.plot.pdf", 
       plot = fig.s7.b.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)


###FIG S8 - NESTLINGS BODY CONDITION MORPHOMETRICS
#FIG S8A - BODY MASS
mass.test <- lmer(FastedBodyMass ~ Treatment + Sex + (1|Parentage), data = morph.df)
summary(mass.test)

fig.s8.a.plot <- 
  ggplot(morph.df, aes(x = Treatment, y = FastedBodyMass, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(5,10.2), breaks = c(6,8, 10)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 10, yend = 10) +
  annotate("text", label = "NS", x = 1.5, y = 10.2, size = 4) +
  ylab("Body mass (g)") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s8.a.plot.pdf", 
       plot = fig.s8.a.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)


#FIG S8B - TARSUS LENGTH
tarsus.test <- lmer(TarsusLength ~ Treatment + Sex + (1|Parentage), data = morph.df)
summary(tarsus.test)

fig.s8.b.plot <- 
  ggplot(morph.df, aes(x = Treatment, y = TarsusLength, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(10,14.3), breaks = c(10,12,14)) +
  theme_bw() +
  annotate("text", label = "NS", x = 1.5, y = 14.3, size = 4) +
  ylab("Tarsus length (mm)") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s8.b.plot.pdf", 
       plot = fig.s8.b.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)


# FIG S8C - WING CHORD
wing.test <- lmer(WingChord ~ Treatment + Sex + (1|Parentage), data = morph.df)
summary(wing.test)

fig.s8.c.plot <- 
  ggplot(morph.df, aes(x = Treatment, y = WingChord, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(13,23.2), breaks = c(14,16,18,20,22)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 22.8, yend = 22.8) +
  annotate("text", label = "NS", x = 1.5, y = 23.2, size = 4) +
  ylab("Wing chord (mm)") +
  xlab("") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s8.c.plot.pdf", 
       plot = fig.s8.c.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

# FIG S9 - GROWTH AND FEEDING RATES OF AX AND CV NESTLINGS
# LOAD SUMMARY SE HELPER FUNCTION
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


## IMPORT RAW DATA GROWTH AND FEEDING DATA 
feeding.growth.df <- read.csv("FeedingGrowthData.csv")

# FILTER TO REMOVE DAY 7 MEASUREMENTS (not a full day of feedings)
feeding.growth.df.filtered <- 
  feeding.growth.df %>%
  filter(Day != 7)

#AX NESTLINGS BY INDIVIDUALS OVER TIME
ax.data.filtered <- filter(feeding.growth.df.filtered, Treatment == "AX")

ax.data.filtered.bm.summary <- summarySE(ax.data.filtered, measurevar="AvgDailyBodyMass", groupvars= c("Day", "AnimalID"))
ax.data.filtered.vol.summary <- summarySE(ax.data.filtered, measurevar="TotalDailyFormulaVolume", groupvars= c("Day", "AnimalID"))

ax.summary <- ax.data.filtered.bm.summary
ax.summary$TotalDailyFormulaVolume <- ax.data.filtered.vol.summary$TotalDailyFormulaVolume
ax.summary$Volume_se <- ax.data.filtered.vol.summary$se
ax.summary <- ax.summary %>% mutate(Treatment = "AX")

#CV NESTLINGS BY INDIVIDUALS OVER TIME
cv.data.filtered <- filter(feeding.growth.df.filtered, Treatment == "CV")

cv.data.filtered.bm.summary <- summarySE(cv.data.filtered, measurevar="AvgDailyBodyMass", groupvars= c("Day", "AnimalID"))
cv.data.filtered.vol.summary <- summarySE(cv.data.filtered, measurevar="TotalDailyFormulaVolume", groupvars= c("Day", "AnimalID"))

cv.summary <- cv.data.filtered.bm.summary
cv.summary$TotalDailyFormulaVolume <- cv.data.filtered.vol.summary$TotalDailyFormulaVolume
cv.summary$Volume_se <- cv.data.filtered.vol.summary$se
cv.summary <- cv.summary %>% mutate(Treatment = "CV")

#COMBINE AX AD CV SUMMARY TABLES
ax.cv.summary.df <- bind_rows(ax.summary, cv.summary)

#FIG S9A - BODY MASS BY TREATMENT GROUP OVER TIME
daily.mass.test <- lmer(AvgDailyBodyMass ~ Treatment * Day + (1|AnimalID), data = ax.cv.summary.df)
summary(daily.mass.test)

fig.s9.a.plot <- 
  ggplot(ax.summary, aes(Day, AvgDailyBodyMass, fill = AnimalID, color = AnimalID)) +
  geom_point(data = ax.summary, shape =21, alpha = 0.7, size = 1.5, color = "#009E73",  fill = "#009E73") +
  geom_path(data = ax.summary, color = "#009E73", aes(Day, AvgDailyBodyMass)) +
  geom_point(data = cv.summary, shape =21, alpha = 0.7, size = 1.5, color = "#D55E00", fill = "#D55E00") +
  geom_path(data = cv.summary, color = "#D55E00", aes(Day, AvgDailyBodyMass)) +
  xlab("Days post-hatch") +
  ylab("Mean daily body mass (g)") +
  annotate("text", x = 3, y = 10, label = "NS", size = 5) +
  scale_x_continuous(limits = c(0, 6), breaks = c(0,1,2,3,4,5,6)) + 
  theme_bw() +
  theme(legend.position = 'none', 
        plot.margin = margin(20,40,5,5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(filename = "fig.s9.a.plot.pdf", 
       plot = fig.s9.a.plot,
       width = 10, 
       height = 10,
       device = cairo_pdf,
       units = "cm",
       dpi = 600)

#FORMULA VOL BY TREATMENT GROUP
daily.intake.test <- lmer(TotalDailyFormulaVolume ~ Treatment * Day + (1|AnimalID), data = ax.cv.summary.df)
summary(daily.intake.test)

fig.s9.b.plot <- 
  ggplot(ax.summary, aes(Day, TotalDailyFormulaVolume, fill = AnimalID, color = AnimalID)) +
  geom_point(data = ax.summary, shape =21, alpha = 0.7, size = 1.5, color = "#009E73",  fill = "#009E73") +
  geom_path(data = ax.summary, color = "#009E73", aes(Day, TotalDailyFormulaVolume)) +
  geom_point(data = cv.summary, shape =21, alpha = 0.7, size = 1.5, color = "#D55E00", fill = "#D55E00") +
  geom_path(data = cv.summary, color = "#D55E00", aes(Day, TotalDailyFormulaVolume)) +
  xlab("Days post-hatch") +
  ylab(expression(paste('Total daily formula intake (', "\u03bc", 'L)'))) +
  annotate("text", x = 3, y = 6500, label = "NS", size = 5) +
  scale_x_continuous(limits = c(0, 6), breaks = c(0,1,2,3,4,5,6)) + 
  theme_bw() +
  theme(legend.position = 'none', 
        plot.margin = margin(20,40,5,5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(filename = "fig.s9.b.plot.pdf", 
       plot = fig.s9.b.plot,
       width = 10, 
       height = 10,
       device = cairo_pdf,
       units = "cm",
       dpi = 600)

# FIG S10 - METABOLOMICS ORDINATIONS
plsda.df <-  read.csv("PLSDA.csv")

#FIG.S10A - PLSDA ORDINATION C18 NEG
fig.s10.a.plot <- 
  ggplot(plsda.df%>%
           filter(Column == "C18NEG"), 
         aes(Comp1, Comp2, fill = Treatment, col = Treatment)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  scale_color_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  stat_ellipse(level = 0.95) +
  xlab("Component 1 (17.5%)") +
  ylab("Component 2 (7.1%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = 'none')

ggsave(filename = "fig.s10.a.plot.pdf", 
       plot = fig.s10.a.plot,
       width = 8, 
       height = 8,
       units = "cm",
       dpi = 600)

#FIG.S10B - PLSDA ORDINATION C18 POS
fig.s10.b.plot <- 
  ggplot(plsda.df%>%
           filter(Column == "C18POS"), aes(Comp1, Comp2, fill = Treatment, col = Treatment)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  scale_color_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  stat_ellipse(level = 0.95) +
  xlab("Component 1 (25.2%)") +
  ylab("Component 2 (6.8%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = 'none')

ggsave(filename = "fig.s10.b.plot.pdf", 
       plot = fig.s10.b.plot,
       width = 8, 
       height = 8,
       units = "cm",
       dpi = 600)

#FIG.S10C - PLSDA ORDINATION HLIC NEG
fig.s10.c.plot <- 
  ggplot(plsda.df%>%
           filter(Column == "HILICNEG"), aes(Comp1, Comp2, fill = Treatment, col = Treatment)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  scale_color_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  stat_ellipse(level = 0.95) +
  xlab("Component 1 (23.0%)") +
  ylab("Component 2 (11.9%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = 'none')

ggsave(filename = "fig.s10.c.plot.pdf", 
       plot = fig.s10.c.plot,
       width = 8, 
       height = 8,
       units = "cm",
       dpi = 600)

#FIG.S10D - PLSDA ORDINATION HLIC POS
fig.s10.d.plot <- 
  ggplot(plsda.df%>%
           filter(Column == "HILICPOS"), aes(Comp1, Comp2, fill = Treatment, col = Treatment)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  scale_color_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  stat_ellipse(level = 0.95) +
  xlab("Component 1 (49.2%)") +
  ylab("Component 2 (5.2%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = 'none')

ggsave(filename = "fig.s10.d.plot.pdf", 
       plot = fig.s10.d.plot,
       width = 8, 
       height = 8,
       units = "cm",
       dpi = 600)




### FIG S12 A - Eicosapentaenoic acid
df <- read.csv("Metabolomics.csv")

fig.s.12.a.plot <- 
  ggplot(df, aes(x = Treatment, y = EicosapentaenoicAcid, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  #geom_hline(aes(yintercept = 0.680488996), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab("Eicosapentaenoic acid") +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 7, yend = 7) +
  annotate("text", label = "*", x = 1.5, y = 7, size = 5) +
  #annotate("segment", x = 1, xend = 2, y = 5.4, yend = 5.4) +
  #annotate("text", label = "**", x = 1.5, y = 5.5, size = 5) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s.12.a.plot.pdf", 
       plot = fig.s.12.a.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 300)

### FIG S12 A - Glyceric acid
fig.s.12.b.plot <- 
  ggplot(df, aes(x = Treatment, y = GlycericAcid, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  #geom_hline(aes(yintercept = 0.680488996), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab("Glyceric acid") +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 7, yend = 7) +
  annotate("text", label = "*", x = 1.5, y = 7, size = 5) +
  #annotate("segment", x = 1, xend = 2, y = 5.4, yend = 5.4) +
  #annotate("text", label = "**", x = 1.5, y = 5.5, size = 5) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.s.12.b.plot.pdf", 
       plot = fig.s.12.b.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 300)



