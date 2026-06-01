# Gut microbiota promote early-life digestive function in songbirds

This repository contains the supplemental R code and input data used to reproduce Figure 1 and Figure 2 from the manuscript:

**Gut microbiota promote early-life digestive function in songbirds**  
Brian K. Trevelline et al.

This repo also includes all data files used to generate supplemental figures. 

## Repository structure

```text
.
├── README.md
├── Figure_1/
│   ├── SuppCode_Fig1.R
│   ├── data/
│   │   ├── Chao1.csv
│   │   ├── MucosaMeasurements.csv
│   │   ├── NestlingData.csv
│   │   ├── RPCA.csv
│   │   ├── SuppTableS4.csv
│   │   └── VilliMeasurements.csv
│   └── output/
│       ├── fig.1.a.plot.pdf
│       ├── fig.1.b.plot.pdf
│       ├── fig.1.c.plot.pdf
│       ├── fig.1.e.plot.pdf
│       ├── fig.1.f.plot.pdf
│       ├── fig.1.g.plot.pdf
│       └── fig.1.h.plot.pdf
└── Figure_2/
    ├── SuppCode_Fig2.R
    ├── data/
    │   └── [Figure 2 input files]
    └── output/
        └── [Figure 2 output files]

## Figure 1 supplemental code

This folder contains the R code, input data, and output files used to reproduce Figure 1.

The main script is:

```text
SuppCode_Fig1.R
```

Dowload all data files and the R script to your working R directory. 

Running the script imports the Figure 1 input data, performs the associated statistical analyses, and writes the Figure 1 panel PDFs.

## Figure 1
## What the code does

The Figure 1 script analyzes how gut microbiota treatment is associated with microbiota composition, bacterial richness, bacterial abundance, and intestinal morphology in nestling songbirds.

The script generates these figure panels:

| Panel | Output file | Analysis |
|---|---|---|
| Figure 1A | `fig.1.a.plot.pdf` | RPCA ordination of microbiota composition |
| Figure 1B | `fig.1.b.plot.pdf` | Chao1 ASV richness |
| Figure 1C | `fig.1.c.plot.pdf` | 16S rRNA copy number per gram intestinal contents |
| Figure 1E | `fig.1.e.plot.pdf` | Villus height |
| Figure 1F | `fig.1.f.plot.pdf` | Crypt depth |
| Figure 1G | `fig.1.g.plot.pdf` | Goblet cells per mm |
| Figure 1H | `fig.1.h.plot.pdf` | Mucosal thickness |

## Input files

### `RPCA.csv`

Used for Figure 1A. This file contains RPCA ordination coordinates and sample metadata. The script plots `RPCA_Axis1` against `RPCA_Axis2`, colors points by treatment, and adds 95% confidence ellipses.

### `Chao1.csv`

Used for Figure 1B. This file contains Chao1 ASV richness estimates. The script removes control samples, tests treatment effects on Chao1 richness, calculates mean richness by individual, and generates a treatment-level richness plot.

### `NestlingData.csv`

Used in the Figure 1B and Figure 1C statistical models. This file contains nestling-level metadata and summary measurements, including treatment, sex, parentage, mean Chao1 richness, and bacterial-abundance estimates.

### `SuppTableS4.csv`

Used for Figure 1C. This file contains qPCR-derived bacterial abundance estimates. The script analyzes log10-transformed 16S rRNA copy number per gram intestinal contents.

### `VilliMeasurements.csv`

Used for Figure 1E and Figure 1F. This file contains intestinal histology measurements, including villus height and crypt depth.

### `MucosaMeasurements.csv`

Used for Figure 1G and Figure 1H. This file contains mucosal measurements, including goblet-cell counts, mucosal perimeter, and mucosal thickness. The script filters this dataset to retain measurements with `mucosal_perimeter_mm > 5`.

## Software requirements

The code was written in R and loads the following packages:

```r
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

To install these packages, run:

```r
install.packages(c(
  "tidyverse",
  "cowplot",
  "lme4",
  "lmerTest",
  "Matrix",
  "afex",
  "rcompanion",
  "patchwork",
  "ggdist",
  "gghalves",
  "emmeans",
  "Cairo",
  "outliers",
  "rmarkdown",
  "knitr"
))
```

## How to run the code

Open R or RStudio and set the working directory to the `Figure_1/` folder.

Then run:

```r
source("SuppCode_Fig1.R")
```

The script assumes that all input CSV files are in the same folder as `SuppCode_Fig1.R`. The generated PDF files are written to the same folder.

## Statistical analyses

The script fits linear models and linear mixed-effects models to test treatment effects. Pairwise treatment comparisons are estimated with `emmeans`, with FDR correction where specified in the script.

### Chao1 ASV richness

The script first tests treatment effects across AX, CV, and INOC samples, excluding CTRL samples:

```r
chao1.test <- lmer(Chao1 ~ Treatment + (1 | IndividualID), data = chao.df.filtered)
emmeans(chao1.test, list(pairwise ~ Treatment), adjust = "fdr")
```

It then tests mean Chao1 richness in AX and CV nestlings while including sex as a fixed effect and parentage as a random effect:

```r
chao1.test <- lmer(MeanChao1 ~ Treatment + Sex + (1 | Parentage), data = nestling.df)
summary(chao1.test)
```

### Bacterial abundance

The script tests treatment effects on qPCR-derived bacterial abundance, excluding BLANK samples:

```r
abund.test <- lm(
  QuantityMeanNGCopyNumberCorrectedPerGramLog10 ~ Treatment,
  data = qpcr.df.filtered
)
emmeans(abund.test, list(pairwise ~ Treatment), adjust = "fdr")
```

It also tests abundance in AX and CV intestinal contents while including sex and parentage:

```r
abund.test <- lmer(
  QuantityMeanNGCopyNumberCorrectedPerGramLog10 ~ Treatment + Sex + (1 | Parentage),
  data = nestling.df
)
summary(abund.test)
```

### Villus height

The script tests effects of treatment, sex, and their interaction on villus height, with parentage and individual identity as random effects:

```r
villi.height.test <- lmer(
  VillusHeight ~ Treatment * Sex + (1 | Parentage) + (1 | IndividualID),
  data = hist.df
)
summary(villi.height.test)
```

### Crypt depth

The script tests effects of treatment, sex, and their interaction on crypt depth, with parentage and individual identity as random effects:

```r
crypt.test <- lmer(
  CryptDepth ~ Treatment * Sex + (1 | Parentage) + (1 | IndividualID),
  data = hist.df
)
summary(crypt.test)
```

### Goblet cells per mm

The script filters mucosal measurements to samples with `mucosal_perimeter_mm > 5`, then tests effects of treatment and sex on goblet cells per mm, with parentage and individual identity as random effects:

```r
goblet.test <- lmer(
  goblet_like_per_mm ~ Treatment + Sex + (1 | Parentage) + (1 | IndividualID),
  data = muc.df.filtered
)
summary(goblet.test)
```

### Mucosal thickness

The script tests effects of treatment, sex, and their interaction on mucosal thickness, with parentage and individual identity as random effects:

```r
thickness.test <- lmer(
  thickness_mean_um ~ Treatment * Sex + (1 | Parentage) + (1 | IndividualID),
  data = muc.df.filtered
)
summary(thickness.test)
```

## Figure panel details

### Figure 1A: RPCA ordination

The script imports `RPCA.csv`, plots `RPCA_Axis1` and `RPCA_Axis2`, colors samples by treatment, and adds 95% confidence ellipses. Axis labels report the percent variation shown by each RPCA axis.

Output:

```text
fig.1.a.plot.pdf
```

### Figure 1B: Chao1 ASV richness

The script imports `Chao1.csv`, removes `CTRL` samples, calculates individual-level mean Chao1 richness, and plots Chao1 richness across AX, CV, and INOC treatments. A dashed horizontal line marks the maximum richness observed in controls.

Output:

```text
fig.1.b.plot.pdf
```

### Figure 1C: bacterial abundance

The script imports `SuppTableS4.csv`, removes `BLANK` samples, and plots log10-transformed mean 16S rRNA copy number per gram intestinal contents across AX, CV, and INOC treatments. A dashed horizontal line marks the maximum blank/control value.

Output:

```text
fig.1.c.plot.pdf
```

### Figure 1E: villus height

The script imports `VilliMeasurements.csv`, fits a mixed-effects model for villus height, calculates median villus height by individual, and plots villus-height distributions for AX and CV nestlings.

Output:

```text
fig.1.e.plot.pdf
```

### Figure 1F: crypt depth

Using `VilliMeasurements.csv`, the script fits a mixed-effects model for crypt depth, calculates median crypt depth by individual, and plots crypt-depth distributions for AX and CV nestlings.

Output:

```text
fig.1.f.plot.pdf
```

### Figure 1G: goblet cells per mm

The script imports `MucosaMeasurements.csv`, filters to measurements with `mucosal_perimeter_mm > 5`, fits a mixed-effects model for goblet cells per mm, and plots goblet-cell density for AX and CV nestlings.

Output:

```text
fig.1.g.plot.pdf
```

### Figure 1H: mucosal thickness

Using the filtered mucosal-measurement data, the script fits a mixed-effects model for mucosal thickness and plots mucosal-thickness distributions for AX and CV nestlings.

Output:

```text
fig.1.h.plot.pdf
```

## Treatment labels

The script uses the following treatment labels:

| Label | Description |
|---|---|
| `AX` | Antibiotic-treated or microbiota-reduced nestlings |
| `CV` | Conventional nestlings |
| `INOC` | Inoculum sample |
| `BLANK` | Blank control |
| `CTRL` | Control sample |

Refer to the manuscript methods for the full experimental design and treatment descriptions.



