################ 
#SUPPLEMENTAL CODE FOR ANALYSES - FIGURE 1
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

#MAIN FIGURES#
#FIGURE 1 - The gut microbiota promotes early-life digestive development in songbirds

#Import QIIME 2 data output
rpca.df <- read.csv("RPCA.csv")

#FIG.1A - RPCA ORDINATION
fig.1.a.plot <- 
  ggplot(rpca.df, aes(RPCA_Axis1, RPCA_Axis2, fill = Treatment, col = Treatment)) +
  geom_point(size = 3, alpha = 0.7, shape = 21) +
  scale_color_manual(values = c("#009E73", "#D55E00", "#0072B2","#000000"), limits=c("AX", "CV","INOC", "BLANK")) +
  scale_fill_manual(values = c("#009E73", "#D55E00", "#0072B2","#A9A9A9"), limits=c("AX", "CV","INOC","BLANK")) +
  stat_ellipse(level = 0.95) +
  xlab("RPCA 1 (72.1%)") +
  ylab("RPCA 2 (27.9%)") +
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
        legend.position = c(0.87, 0.2),
        legend.box.background = element_rect(color = "black", linetype = "solid"),
        legend.box.margin = margin(-5, -2, -5, -5))

ggsave(filename = "fig.1.a.plot.pdf", 
       plot = fig.1.a.plot,
       width = 9.5, 
       height = 8,
       units = "cm",
       dpi = 600)

#FIG.1B - CHAO 1 ASV RICHNESS
#Import CHAO1 DATA
chao.df <- read.csv("Chao1.csv")

#TEST CHAO1 ACROSS AX, CV, INOC
chao.df.filtered <-
  chao.df %>%
  filter(Treatment != "CTRL")

#TEST EFFECT OF TREATMENT ON CHAO1 IN AX, CV, INOC
chao1.test <- lmer(Chao1 ~ Treatment + (1|IndividualID), data = chao.df.filtered)
emmeans(chao1.test, list(pairwise ~ Treatment), adjust = "fdr")

chao.df.mean <- aggregate(Chao1 ~ IndividualID + Treatment, chao.df.filtered, mean)

#Import NESTLING DATA
nestling.df <- read.csv("NestlingData.csv")

#TEST EFFECT OF TREATMENT ON MEAN CHAO1 IN AX AND CV NESTLINGS WITH FIXED EFFECT OF SEX AND RANDOM EFFECT OF PARENTAGE
chao1.test <- lmer(MeanChao1 ~ Treatment + Sex + (1|Parentage), data = nestling.df)
summary(chao1.test)

fig.1.b.plot <- 
  ggplot(chao.df.mean, aes(x = Treatment, y = Chao1, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  geom_hline(aes(yintercept = 7), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00", "#0072B2"), limits = c("AX", "CV", "INOC")) +
  scale_x_discrete(labels = c("AX","CV","INOC")) +
  scale_y_continuous(limits = c(0,51)) +
  theme_bw() +
  annotate("segment", x = 2, xend = 3, y = 49.5, yend = 49.5) +
  annotate("text", label = "NS", x = 2.5, y = 50.9, size = 3.4) +
  annotate("segment", x = 1, xend = 2, y = 39.5, yend = 39.5) +
  annotate("text", label = "*", x = 1.5, y = 39.5, size = 5) +
  ylab("Mean Chao1 ASV richness") +
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

ggsave(filename = "fig.1.b.plot.pdf", 
       plot = fig.1.b.plot,
       width = 6, 
       height = 8,
       units = "cm",
       dpi = 600)

#FIG 1C - Bacterial Abundance (16S copy number/g intestinal contents)
#TEST EFFECT OF TREATMENT ON BACTERIAL ABUNDNANCE IN AX, CV, INOC
qpcr.df <- read.csv("SuppTableS4.csv")

qpcr.df.filtered <-
  qpcr.df %>% 
  filter(Treatment != "BLANK")

abund.test <- lm(QuantityMeanNGCopyNumberCorrectedPerGramLog10 ~ Treatment, data = qpcr.df.filtered)
emmeans(abund.test, list(pairwise ~ Treatment), adjust = "fdr")

#RUN LMER FOR AX AND CV INTESTINAL CONTENTS ONLY WITH SEX AS FIXED EFFECT AND PARENTAGE AS RANDOM EFFECT
abund.test <- lmer(QuantityMeanNGCopyNumberCorrectedPerGramLog10 ~ Treatment + Sex + (1|Parentage), data = nestling.df)
summary(abund.test)

fig.1.c.plot <- 
  ggplot(qpcr.df.filtered, aes(x = Treatment, y = as.numeric(QuantityMeanNGCopyNumberCorrectedPerGramLog10), fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  geom_hline(aes(yintercept = 0.985331395), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00", "#0072B2"), limits = c("AX", "CV", "INOC")) +
  scale_x_discrete(labels = c("AX","CV","INOC")) +
  scale_y_continuous(limits = c(0,6.1)) +
  ylab(expression(paste("log"["10"], "(mean 16S rRNA copies/g)"))) +
  theme_bw() +
  annotate("segment", x = 2, xend = 3, y = 5.9, yend = 5.9) +
  annotate("text", label = "**", x = 2.5, y = 5.9, size = 5) +
  annotate("segment", x = 1, xend = 2, y = 4.7, yend = 4.7) +
  annotate("text", label = "**", x = 1.5, y = 4.7, size = 5) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.1.c.plot.pdf", 
       plot = fig.1.c.plot,
       width = 6, 
       height = 8,
       units = "cm",
       dpi = 600)

#FIG 1E-H - HISTOLOGY DATA
hist.df <- read.csv("VilliMeasurements.csv")

#FIG 1E - VILLUS HEIGHT
villi.height.test <- lmer(VillusHeight ~ Treatment * Sex + (1|Parentage) + (1|IndividualID), data = hist.df)
summary(villi.height.test)

#MEDIAN VILLUS HEIGHT
med.vil.len <- aggregate(VillusHeight ~ IndividualID + Treatment, data = hist.df, FUN = median, na.rm = TRUE)

fig.1.e.plot <- 
  ggplot(hist.df, aes(Treatment, VillusHeight, fill = Treatment)) +
  ggdist::stat_halfeye(width = 0.7) + 
  #geom_boxplot() +
  geom_half_point(side = "l", alpha = 0.1, range_scale = 0.35, color = "#000000") +
  #geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,1000), breaks = c(0,200,400,600,800)) +
  annotate("segment", x = 1.0, xend = 1.3, y = 427.4999, yend = 427.4999) +
  annotate("segment", x = 1.0, xend = 1.3, y = 253.1266, yend = 253.1266) +
  annotate("segment", x = 1.0, xend = 1.3, y = 375.5115, yend = 375.5115) +
  annotate("segment", x = 1.0, xend = 1.3, y = 391.3025, yend = 391.3025) +
  annotate("segment", x = 1.0, xend = 1.3, y = 369.8048, yend = 369.8048) +
  annotate("segment", x = 1.0, xend = 1.3, y = 413.0883, yend = 413.0883) +
  annotate("segment", x = 2.0, xend = 2.3, y = 387.5180, yend = 387.5180) +
  annotate("segment", x = 2.0, xend = 2.3, y = 458.2572, yend = 458.2572) +
  annotate("segment", x = 2.0, xend = 2.3, y = 475.5933, yend = 475.5933) +
  annotate("segment", x = 2.0, xend = 2.3, y = 411.1535, yend = 411.1535) +
  annotate("segment", x = 2.0, xend = 2.3, y = 446.9223, yend = 446.9223) +
  annotate("segment", x = 2.0, xend = 2.3, y = 446.2550, yend = 446.2550) +
  #annotate("segment", x = 1, xend = 2, y = 990, yend = 990) +
  annotate("text", label = "**", x = 1.5, y = 1000, size = 5) +
  theme_bw() +
  ylab(expression(paste('Villi height (', "\u03bc", 'm)'))) +
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

ggsave(filename = "fig.1.e.plot.pdf", 
       plot = fig.1.e.plot,
       device = cairo_pdf,
       width = 5, 
       height = 7.5,
       units = "cm",
       dpi = 600)

#FIG. 1F - CRYPT DEPTH
crypt.test <- lmer(CryptDepth ~ Treatment * Sex + (1|Parentage) + (1|IndividualID), data = hist.df)
summary(crypt.test)

med.crypt.dep <- aggregate(CryptDepth ~ IndividualID + Treatment, data = hist.df, FUN = median, na.rm = TRUE)

fig.1.f.plot <- 
  ggplot(hist.df, aes(Treatment, CryptDepth, fill = Treatment)) +
  ggdist::stat_halfeye(width = 0.7) + 
  #geom_boxplot() +
  geom_half_point(side = "l", alpha = 0.1, range_scale = 0.35, color = "#000000") +
  #geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  #scale_y_continuous(limits = c(0,225)) +
  annotate("segment", x = 1.0, xend = 1.3, y = 144.260, yend = 144.260) +
  annotate("segment", x = 1.0, xend = 1.3, y = 123.550, yend = 123.550) +
  annotate("segment", x = 1.0, xend = 1.3, y = 144.220, yend = 144.220) +
  annotate("segment", x = 1.0, xend = 1.3, y = 185.540, yend = 185.540) +
  annotate("segment", x = 1.0, xend = 1.3, y = 159.370, yend = 159.370) +
  annotate("segment", x = 1.0, xend = 1.3, y = 139.775, yend = 139.775) +
  annotate("segment", x = 2.0, xend = 2.3, y = 188.555, yend = 188.555) +
  annotate("segment", x = 2.0, xend = 2.3, y = 193.370, yend = 193.370) +
  annotate("segment", x = 2.0, xend = 2.3, y = 184.500, yend = 184.500) +
  annotate("segment", x = 2.0, xend = 2.3, y = 165.660, yend = 165.660) +
  annotate("segment", x = 2.0, xend = 2.3, y = 174.900, yend = 174.900) +
  annotate("segment", x = 2.0, xend = 2.3, y = 192.955, yend = 192.955) +
  #annotate("segment", x = 1, xend = 2, y = 520, yend = 520) +
  annotate("text", label = "**", x = 1.5, y = 525, size = 5) +
  theme_bw() +
  ylab(expression(paste('Crypt depth (', "\u03bc", 'm)'))) +
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

ggsave(filename = "fig.1.f.plot.pdf", 
       plot = fig.1.f.plot,
       device = cairo_pdf,
       width = 5, 
       height = 7.5,
       units = "cm",
       dpi = 600)

#FIG. 1G - GOBLET CELLS
muc.df <- read.csv("MucosaMeasurements.csv")

muc.df.filtered <-
  muc.df %>%
  filter(mucosal_perimeter_mm > 5) 

#goblet_like_per_mm
goblet.test <- lmer(goblet_like_per_mm ~ Treatment + Sex + (1|Parentage) + (1|IndividualID), data = muc.df.filtered)
summary(goblet.test)

#MEDIAN BY INDIVIDUAL
med.muc.gob <- aggregate(goblet_like_per_mm ~ IndividualID + Treatment, data = muc.df.filtered, FUN = median, na.rm = TRUE)

fig.1.g.plot <- 
  ggplot(muc.df.filtered, aes(Treatment, goblet_like_per_mm, fill = Treatment)) +
  ggdist::stat_halfeye(width = 0.7) + 
  #geom_boxplot() +
  geom_half_point(side = "l", alpha = 0.5, range_scale = 0.35, color = "#000000") +
  #geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,250)) +
  annotate("segment", x = 1.0, xend = 1.3, y = 45.71718, yend = 45.71718) +
  annotate("segment", x = 1.0, xend = 1.3, y = 38.49662, yend = 38.49662) +
  annotate("segment", x = 1.0, xend = 1.3, y = 53.32439, yend = 53.32439) +
  annotate("segment", x = 1.0, xend = 1.3, y = 76.52062, yend = 76.52062) +
  annotate("segment", x = 1.0, xend = 1.3, y = 101.48446, yend = 101.48446) +
  annotate("segment", x = 1.0, xend = 1.3, y = 30.40627, yend = 30.40627) +
  annotate("segment", x = 2.0, xend = 2.3, y = 59.43893, yend = 59.43893) +
  annotate("segment", x = 2.0, xend = 2.3, y = 46.54737, yend = 46.54737) +
  annotate("segment", x = 2.0, xend = 2.3, y = 121.43877, yend = 121.43877) +
  annotate("segment", x = 2.0, xend = 2.3, y = 51.60546, yend = 51.60546) +
  annotate("segment", x = 2.0, xend = 2.3, y = 55.54546, yend = 55.54546) +
  annotate("segment", x = 2.0, xend = 2.3, y = 93.87913, yend = 93.87913) +
  #annotate("segment", x = 1, xend = 2, y = 248, yend = 248) +
  annotate("text", label = "*", x = 1.5, y = 250, size = 5) +
  theme_bw() +
  ylab("Mean goblet cells/mm") +
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

ggsave(filename = "fig.1.g.plot.pdf", 
       plot = fig.1.g.plot,
       device = cairo_pdf,
       width = 5, 
       height = 7.5,
       units = "cm",
       dpi = 600)

#FIG. 1H - MUCOSAL THICKNESS
thickness.test <- lmer(thickness_mean_um ~ Treatment * Sex + (1|Parentage) + (1|IndividualID), data = muc.df.filtered)
summary(thickness.test)

#MEDIAN BY INDIVIDUAL
med.muc.thick <- aggregate(thickness_mean_um ~ IndividualID + Treatment, data = muc.df.filtered, FUN = median, na.rm = TRUE)

fig.1.h.plot <- 
  ggplot(muc.df.filtered, aes(Treatment, thickness_mean_um, fill = Treatment)) +
  ggdist::stat_halfeye(width = 0.7) + 
  #geom_boxplot() +
  geom_half_point(side = "l", alpha = 0.5, range_scale = 0.35, color = "#000000") +
  #geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("#009E73", "#D55E00"), limits=c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  #scale_y_continuous(limits = c(1.5,3.0), breaks = c(2.0,3.0)) +
  annotate("segment", x = 1.0, xend = 1.3, y = 1.962736, yend = 1.962736) +
  annotate("segment", x = 1.0, xend = 1.3, y = 1.707603, yend = 1.707603) +
  annotate("segment", x = 1.0, xend = 1.3, y = 2.087693, yend = 2.087693) +
  annotate("segment", x = 1.0, xend = 1.3, y = 1.958111, yend = 1.958111) +
  annotate("segment", x = 1.0, xend = 1.3, y = 2.055509, yend = 2.055509) +
  annotate("segment", x = 1.0, xend = 1.3, y = 1.693997, yend = 1.693997) +
  annotate("segment", x = 2.0, xend = 2.3, y = 1.952202, yend = 1.952202) +
  annotate("segment", x = 2.0, xend = 2.3, y = 2.173275, yend = 2.173275) +
  annotate("segment", x = 2.0, xend = 2.3, y = 2.310892, yend = 2.310892) +
  annotate("segment", x = 2.0, xend = 2.3, y = 1.849653, yend = 1.849653) +
  annotate("segment", x = 2.0, xend = 2.3, y = 2.163801, yend = 2.163801) +
  annotate("segment", x = 2.0, xend = 2.3, y = 2.403334, yend = 2.403334) +
  #annotate("segment", x = 1, xend = 2, y = 2.99, yend = 2.99) +
  annotate("text", label = "**", x = 1.5, y = 3, size = 5) +
  theme_bw() +
  ylab(expression(paste('Mean mucosa thickness (', "\u03bc", 'm)'))) +
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

ggsave(filename = "fig.1.h.plot.pdf", 
       plot = fig.1.h.plot,
       device = cairo_pdf,
       width = 5, 
       height = 7.5,
       units = "cm",
       dpi = 600)

