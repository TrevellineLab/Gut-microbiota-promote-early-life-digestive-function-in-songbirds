################ 
#SUPPLEMENTAL CODE FOR ANALYSES - FIGURE 2  #################################
#Gut microbiota promote early-life digestive function in songbirds
#Brian K. Trevelline
#23 MAY 2026

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

### FIGURE 2

### LOAD MORPHOLOGY DATA
df <- read.csv("MorphologyData.csv")

### FIG 2A - LIVER MASS

#REMOVE OUTLIER NL113
grubbs.test(df$LiverMassMg)
df.filtered <- 
  df %>%
  filter(AnimalID != "NL113")

fig.2.a.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = LiverMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,450), breaks = c(0,200,400)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 435, yend = 435) +
  annotate("text", label = "*", x = 1.5, y = 435, size = 5) +
  ylab("Liver mass (mg)") +
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

ggsave(filename = "fig.2.a.plot.pdf", 
       plot = fig.2.a.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

#FIG 2B - LARGE INTESTINE MASS

fig.2.b.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = LIMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,250), breaks = c(0,100, 200)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 240, yend = 240) +
  annotate("text", label = "*", x = 1.5, y = 240, size = 5) +
  ylab("Large intestine mass (mg)") +
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

ggsave(filename = "fig.2.b.plot.pdf", 
       plot = fig.2.b.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

#FIG 2C - SMALL INTESTINE MASS

fig.2.c.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = SIMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(300,750), breaks = c(300, 500, 700)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 730, yend = 730) +
  annotate("text", label = "NS", x = 1.5, y = 745, size = 4) +
  ylab("Small intestine mass (mg)") +
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

ggsave(filename = "fig.2.c.plot.pdf", 
       plot = fig.2.c.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)


### FIG 2D - Pancreas Mass

fig.2.d.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = PancreasMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,112), breaks = c(0,50,100)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 106, yend = 106) +
  annotate("text", label = "NS", x = 1.5, y = 111.5, size = 4) +
  ylab("Pancreas mass (mg)") +
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

ggsave(filename = "fig.2.d.plot.pdf", 
       plot = fig.2.d.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)


#FIG 2X - INTESTINAL MASS

fig.2.x.plot <- 
  ggplot(df, aes(x = Treatment, y = IntestinalMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  #scale_y_continuous(limits = c(0,250), breaks = c(0,100, 200)) +
  theme_bw() +
  #annotate("segment", x = 1, xend = 2, y = 240, yend = 240) +
  #annotate("text", label = "*", x = 1.5, y = 240, size = 5) +
  ylab("Intestinal mass (mg)") +
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

ggsave(filename = "fig.2.x.plot.pdf", 
       plot = fig.2.x.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

### FIG 2E - METABOLOMICS VOLCANO PLOT
df <- read.csv("MetaboliteVolcano.csv")

fig.2.e.plot <-
  ggplot(df, aes(x = mean_difference_CV_minus_AX, y = -log10(p_value_group), fill = SigDirectionFC, label = Label)) +
  geom_vline(xintercept = 0, color = "BLACK", linetype = 1) +
  geom_hline(yintercept = 1.30102999566, color = "#808080", linetype = 2) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_text(hjust=-0.2, vjust=0.5, size =2) +
  scale_fill_manual(values = c("#009E73","#D55E00","#808080"), limits = c("SigDown", "SigUp", "N")) +
  theme_bw() +
  xlab("mean difference CV minus AX") +
  ylab("-log10(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = "fig.2.e.plot.pdf", 
       plot = fig.2.e.plot,
       width = 14, 
       height = 10,
       units = "cm",
       dpi = 600)

### FIGURE 1E INSET
#Load data frame from Python stats output
#Define the paired vectors

df <- read.csv("Fig2E_INSET.csv")

#Reshape the data to "long" format for ggplot2
df.pivot <- pivot_longer(df, 
                         c(AX, CV),
                         names_to = "Variable", 
                         values_to = "Value")

fig.2.e.inset <- 
  ggplot(df.pivot, aes(x = Variable, y = Value, fill = Variable)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  geom_line(aes(group = ID), color = "BLACK", linetype = "solid") +
  scale_y_continuous(limits = c(-0.4,0.45), breaks = c(-0.4,-0.2,0,0.2,0.4)) +
  scale_x_discrete(labels = c("AX","CV")) +
  ylab("Mean relative\npeak intensity") +
  theme_classic() +
  annotate("text", label = "**", x = 1.5, y = 0.44, size = 5) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        plot.background = element_rect(fill = "transparent"),
        legend.position = 'none')

ggsave(filename = "fig.2.e.inset.pdf", 
       plot = fig.2.e.inset,
       width = 4, 
       height = 5,
       units = "cm",
       dpi = 600)

# RELATIVE PEAK INTENSITY BOXPLOTS
df <- read.csv("Metabolomics.csv")

### FIG 2F - Cholic acid
fig.2.f.plot <- 
  ggplot(df, aes(x = Treatment, y = CholicAcid, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  #geom_hline(aes(yintercept = 0.680488996), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab("Cholic acid") +
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

ggsave(filename = "fig.2.f.plot.pdf", 
       plot = fig.2.f.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 300)

### FIG 2G - Beta-muricholic acid
fig.2.g.plot <- 
  ggplot(df, aes(x = Treatment, y = BetaMuricholicAcid01, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  #geom_hline(aes(yintercept = 0.680488996), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab("Beta-muricholic acid") +
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

ggsave(filename = "fig.2.g.plot.pdf", 
       plot = fig.2.g.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 300)

### FIG 2H - Taurochenodeoxycholic acid
fig.2.h.plot <- 
  ggplot(df, aes(x = Treatment, y = TaurochenodeoxycholicAcidSalt, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  #geom_hline(aes(yintercept = 0.680488996), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab("Taurochenodeoxycholic acid") +
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
        axis.title.y = element_text(size = 12, hjust = 0.2),
        legend.position = 'none')

ggsave(filename = "fig.2.h.plot.pdf", 
       plot = fig.2.h.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 300)

### FIG 2I - Glycholic acid
fig.2.i.plot <- 
  ggplot(df, aes(x = Treatment, y = GlycocholicAcid, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  #geom_hline(aes(yintercept = 0.680488996), col = "#000000", linetype = 2) + #Y intercept value represents max abundance in controls
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab("Glycholic acid") +
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

ggsave(filename = "fig.2.i.plot.pdf", 
       plot = fig.2.i.plot,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 300)








