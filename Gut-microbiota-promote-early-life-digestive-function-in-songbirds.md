Gut microbiota promote early-life digestive function in songbirds
================
Brian K. Trevelline
March 21, 2026

``` r
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
```

Figure 1 - The gut microbiota promotes early-life digestive development
in songbirds

``` r
#FIG.1A - RPCA ORDINATION
#Import QIIME 2 data output
rpca.df <- read.csv("GFHOSP_QIIME2_RPCA_FIG1A.csv")

#Plot data
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
```

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## ℹ Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](Gut-microbiota-promote-early-life-digestive-function-in-songbirds_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#FIG.1B - CHAO 1 ASV RICHNESS
#Import CHAO1 AND NESTLING DATA
chao.df <- read.csv("GFHOSP_QIIME2_Chao1.csv")
nestling.df <- read.csv("GFHOSP_NESTLING_DATA.csv")

#TEST CHAO1 ACROSS AX, CV, INOC
chao.df.filtered <-
  chao.df %>%
  filter(Treatment != "CTRL")

chao1.test <- lmer(Chao1 ~ Treatment + (1|IndividualID), data = chao.df.filtered)
emmeans(chao1.test, list(pairwise ~ Treatment), adjust = "fdr")
```

    ## $`emmeans of Treatment`
    ##  Treatment emmean   SE    df lower.CL upper.CL
    ##  AX          17.3 1.69  26.7     13.9     20.8
    ##  CV          23.2 1.80  26.0     19.5     27.0
    ##  INOC        30.2 4.55 115.8     21.2     39.3
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of Treatment`
    ##  1         estimate   SE    df t.ratio p.value
    ##  AX - CV      -5.91 2.47  26.3  -2.393  0.0362
    ##  AX - INOC   -12.91 4.86 108.8  -2.657  0.0272
    ##  CV - INOC    -7.00 4.90 107.2  -1.429  0.1559
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## P value adjustment: fdr method for 3 tests

``` r
chao.df.mean <- aggregate(Chao1 ~ IndividualID + Treatment, chao.df.filtered, mean)

#TEST EFFECT OF TREATMENT ON MEAN CHAO1 IN AX AND CV NESTLINGS WITH RANDOM EFFECTS FOR PARENTAGE AND SEX
chao1.test <- lmer(MeanChao1 ~ Treatment + (1|ParentagePair) + (1|Sex), data = nestling.df)
chao1.null <-lmer(MeanChao1 ~ 1 + (1|ParentagePair) + (1|Sex), data = nestling.df)
anova(chao1.test, chao1.null)
```

    ## Data: nestling.df
    ## Models:
    ## chao1.null: MeanChao1 ~ 1 + (1 | ParentagePair) + (1 | Sex)
    ## chao1.test: MeanChao1 ~ Treatment + (1 | ParentagePair) + (1 | Sex)
    ##            npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
    ## chao1.null    4 195.93 201.26 -93.966   187.93                       
    ## chao1.test    5 192.42 199.08 -91.209   182.42 5.5151  1    0.01885 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
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
```

![](Gut-microbiota-promote-early-life-digestive-function-in-songbirds_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#FIG 1C
#Bacterial Abundance (16S copy number/g intestinal contents)

#Import data
qpcr.df <- read.csv("SuppTableS4_031326.csv")

#TEST EFFECT OF TREATMENT ON BACTERIAL ABUNDNANCE IN AX, CV, INOC
qpcr.df.filtered <-
  qpcr.df %>% 
  filter(Treatment != "BLANK")

abund.test <- lm(QuantityMeanNGCopyNumberCorrectedPerGramLog10 ~ Treatment, data = qpcr.df.filtered)
emmeans(abund.test, list(pairwise ~ Treatment), adjust = "fdr")
```

    ## $`emmeans of Treatment`
    ##  Treatment emmean    SE df lower.CL upper.CL
    ##  AX          1.82 0.245 33     1.32     2.32
    ##  CV          2.87 0.264 33     2.34     3.41
    ##  INOC        4.37 0.336 33     3.69     5.05
    ## 
    ## Confidence level used: 0.95 
    ## 
    ## $`pairwise differences of Treatment`
    ##  1         estimate    SE df t.ratio p.value
    ##  AX - CV      -1.05 0.360 33  -2.918  0.0063
    ##  AX - INOC    -2.55 0.416 33  -6.127  <.0001
    ##  CV - INOC    -1.50 0.427 33  -3.509  0.0020
    ## 
    ## P value adjustment: fdr method for 3 tests

``` r
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
```

![](Gut-microbiota-promote-early-life-digestive-function-in-songbirds_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
