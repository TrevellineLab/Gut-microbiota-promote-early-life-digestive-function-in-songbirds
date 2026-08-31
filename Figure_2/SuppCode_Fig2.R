################ 
#SUPPLEMENTAL CODE FOR ANALYSES - FIGURE 2  #################################
#Gut microbiota promote early-life digestive function in songbirds
#Brian K. Trevelline
#6 JULY 2026

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

liver.test <- lmer(log10(LiverMassMg) ~ Treatment + log10(TarsusLength) * Sex + (1|Parentage), data = df.filtered)
summary(liver.test)

fig.2.a.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = LiverMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,450), breaks = c(0,200,400)) +
  theme_bw() +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

#FIG 2B - LARGE INTESTINE MASS
li.mass.test <- lmer(log10(LIMassMg) ~ Treatment + log10(TarsusLength) * Sex + (1|Parentage), data = df)
summary(li.mass.test)

fig.2.b.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = LIMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,250), breaks = c(0,100, 200)) +
  theme_bw() +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

#FIG 2C - SMALL INTESTINE MASS
si.mass.test <- lmer(log10(SIMassMg) ~ Treatment + log10(TarsusLength) * Sex + (1|Parentage), data = df)
summary(si.mass.test)

fig.2.c.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = SIMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(300,750), breaks = c(300, 500, 700)) +
  theme_bw() +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)


### FIG 2D - Pancreas Mass
panc.mass.test <- lmer(log10(PancreasMassMg) ~ Treatment + log10(TarsusLength) * Sex + (1|Parentage), data = df)
summary(panc.mass.test)

fig.2.d.plot <- 
  ggplot(df.filtered, aes(x = Treatment, y = PancreasMassMg, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(0,112), breaks = c(0,50,100)) +
  theme_bw() +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

### FIGURE 1E - METABOLOMIC VOLCANO PLOT
volcano.df <- read.csv("lme_results_all_group_plus_sex_fixed_with_mean_difference.csv")

fig.2.e.plot <- 
  ggplot(volcano.df, aes(x=mean_difference_CV_minus_AX, y = -log10(p_value_group), fill = p_value_group_significance_direction, label = Label)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.7, color = "BLACK") +
  geom_text(hjust=0, vjust=1.7, size =2.2) +
  scale_fill_manual(values = c("#009E73","#D55E00","#808080"), limits = c("Negative","Positive", "NS")) +
  geom_vline(xintercept = 0, color = "BLACK", linetype = 1) +
  geom_hline(yintercept = 1.30102999566, color = "#808080", linetype = 2) +
  xlab("Mean difference (CV-AX)") +
  ylab(expression(-log[10](p-value))) +
  scale_y_continuous(limits = c(0.0, 3.05), breaks = c(0.0,0.5,1.0,1.5,2.0,2.5,3.0))+
  scale_x_continuous(limits = c(-0.875,0.875), breaks = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8)) +
  theme_bw() +
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
       device = cairo_pdf,
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

fig.2.e.inset.plot <- 
  ggplot(df.pivot, aes(x = Variable, y = Value, fill = Variable)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.7, size = 1) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  geom_line(aes(group = ID), color = "black", linetype = "solid") +
  scale_y_continuous(limits = c(-0.4,0.45), breaks = c(-0.4,-0.2,0,0.2,0.4)) +
  scale_x_discrete(labels = c("AX","CV")) +
  ylab("Mean diff. bile acids") +
  theme_classic() +
  annotate("text", label = "**", x = 1.5, y = 0.445, size = 4) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.background = element_rect(fill = "white"),
        legend.position = 'none')

ggsave(filename = "fig.2.e.inset.plot.pdf", 
       plot = fig.2.e.inset.plot,
       device = cairo_pdf,
       width = 3.1, 
       height = 3.5,
       units = "cm",
       dpi = 600)

# RELATIVE PEAK INTENSITY BOXPLOTS
metabolomics.df <- read.csv("Metabolomics.csv")

### FIG 2F - Cholic acid
fig.2.f.plot <- 
  ggplot(metabolomics.df <- read.csv("Metabolomics.csv"), aes(x = Treatment, y = CholicAcid, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7.0), breaks = c(5.0,6.0,7.0)) +
  ylab(expression(-log[10]("CA rel. intensity"))) +
  theme_bw() +
  annotate("text", label = "**", x = 1.5, y = 7, size = 5) +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

### FIG 2G - Glycholic acid
fig.2.g.plot <- 
  ggplot(metabolomics.df <- read.csv("Metabolomics.csv"), aes(x = Treatment, y = GlycocholicAcid, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab(expression(-log[10]("GCA rel. intensity"))) +
  theme_bw() +
  annotate("text", label = "*", x = 1.5, y = 7, size = 5) +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

### FIG 2H - Taurochenodeoxycholic acid
fig.2.h.plot <- 
  ggplot(metabolomics.df <- read.csv("Metabolomics.csv"), aes(x = Treatment, y = TaurochenodeoxycholicAcidSalt, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab(expression(-log[10]("TCDCA rel. intensity"))) +
  theme_bw() +
  annotate("text", label = "*", x = 1.5, y = 7, size = 5) +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

### FIG 2G - Beta-muricholic acid
fig.2.i.plot <- 
  ggplot(metabolomics.df <- read.csv("Metabolomics.csv"), aes(x = Treatment, y = BetaMuricholicAcid01, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  scale_fill_manual(values = c("#009E73","#D55E00"), limits = c("AX", "CV")) +
  scale_x_discrete(labels = c("AX","CV")) +
  scale_y_continuous(limits = c(4.5,7), breaks = c(5,6,7)) +
  ylab(expression(-log[10]("BMA rel. intensity"))) +
  theme_bw() +
  annotate("text", label = "*", x = 1.5, y = 7, size = 5) +
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
       device = cairo_pdf,
       width = 3.5, 
       height = 7,
       units = "cm",
       dpi = 600)

