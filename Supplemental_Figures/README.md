# Supplemental Figures Code

This folder contains the code used to generate supplemental figures and supporting analyses for the manuscript:

**Gut microbiota promote early-life digestive function in songbirds**

The main script is:

```text
SuppCode_SuppFigs.R
```

This script performs statistical analyses and generates publication-ready PDF figure panels for the supplemental results. The analyses include microbiome ordinations, RPCA distance comparisons, ASV transfer-efficiency calculations, histological measurements, organ and body-size morphometrics, growth and feeding trajectories, and metabolomics visualizations.

Fig. S11 is generated using the python script Code_For_Plots_FigS11.py.

## Overview

The script is organized by supplemental figure. Each section imports the relevant input data, performs any required statistical test or data transformation, generates a figure, and saves the result as a PDF.

The workflow is mostly linear: running the script from top to bottom will reproduce the supplemental figure panels, assuming all required input files are present in the working directory.

## Required R packages

The script uses the following R packages:

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

Most figures are generated with `ggplot2`, with additional support from `ggdist`, `gghalves`, and `Cairo` for specialized plotting and PDF output. Mixed-effects models are fit with `lme4` and `lmerTest`.

## Input files

The script expects the following input files to be available in the working directory:

```text
RPCAordinationFIGS1.csv
RPCAdistancesFIGS2.csv
RPCAdistancesFIGS3.csv
ASVRelativeAbundances_Genus.csv
INOC-CV_ASVs.csv
VilliMeasurements.csv
MorphologyData.csv
FeedingGrowthData.csv
PLSDA.csv
lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
Metabolomics.csv
```

## Output files

The script writes figure panels as PDF files, including:

```text
fig.s1.plot.pdf
fig.s2.a.plot.pdf
fig.s2.b.plot.pdf
fig.s3.plot.pdf
fig.s4.a.plot.pdf
fig.s4.b.plot.pdf
fig.s6.plot.pdf
fig.s7.a.plot.pdf
fig.s7.b.plot.pdf
fig.s8.a.plot.pdf
fig.s8.b.plot.pdf
fig.s8.c.plot.pdf
fig.s9.a.plot.pdf
fig.s9.b.plot.pdf
fig.s10.a.plot.pdf
fig.s10.b.plot.pdf
fig.s10.c.plot.pdf
fig.s10.d.plot.pdf
fig.s11.plot.pdf
fig.s.12.a.plot.pdf
fig.s.12.b.plot.pdf
```

The ASV transfer-efficiency section also writes two CSV files:

```text
TransferredASVs.csv
NonTransferredASVs.csv
```

These files contain ASVs classified as successfully transferred or not transferred from inocula to conventionalized nestlings.

## Supplemental Figure S1: RPCA ordination

This section imports RPCA ordination coordinates and plots microbial community composition across all sample types and time points. Samples are plotted by the first two RPCA axes, with point shape and color indicating treatment or sample category.

Input:

```text
RPCAordinationFIGS1.csv
```

Output:

```text
fig.s1.plot.pdf
```

## Supplemental Figure S2: RPCA distance comparisons

This section compares pairwise RPCA distances between treatment groups and reference sample types.

Panel A compares AX and CV samples to blanks.

Panel B compares AX and CV samples to inocula.

The script uses permutation-based pairwise tests with FDR correction, then generates boxplots with jittered individual points.

Input:

```text
RPCAdistancesFIGS2.csv
```

Outputs:

```text
fig.s2.a.plot.pdf
fig.s2.b.plot.pdf
```

## Supplemental Figure S3: Recipient versus non-recipient RPCA distances

This section tests whether inoculum-recipient status predicts RPCA distance to inoculum samples while accounting for repeated measures or paired structure using a mixed-effects model.

Input:

```text
RPCAdistancesFIGS3.csv
```

Output:

```text
fig.s3.plot.pdf
```

## Supplemental Figure S4: Microbiota composition and ASV transfer efficiency

Panel A plots genus-level microbial relative abundances across sample groups.

Panel B calculates the proportion of inoculum ASVs detected in each corresponding CV nestling, summarizes transfer efficiency across individuals, and plots the distribution of transfer efficiencies.

This section also extracts and writes tables of transferred and non-transferred ASVs.

Inputs:

```text
ASVRelativeAbundances_Genus.csv
INOC-CV_ASVs.csv
```

Outputs:

```text
fig.s4.a.plot.pdf
fig.s4.b.plot.pdf
TransferredASVs.csv
NonTransferredASVs.csv
```

## Supplemental Figure S6: Villus width

This section analyzes intestinal villus width using a mixed-effects model with treatment and sex as fixed effects and parentage and individual identity as random effects.

Model:

```r
MeanVillusWidth ~ Treatment + Sex + (1 | Parentage) + (1 | IndividualID)
```

The figure shows villus-width distributions for AX and CV nestlings.

Input:

```text
VilliMeasurements.csv
```

Output:

```text
fig.s6.plot.pdf
```

## Supplemental Figure S7: Organ mass

This section analyzes gizzard and heart mass using mixed-effects models. Gizzard mass is modeled with a treatment-by-sex interaction and tarsus length as a body-size covariate. Heart mass is modeled with treatment, sex, and tarsus length as predictors.

Models:

```r
GizzardMassMg ~ Treatment + Sex + TarsusLength + (1 | Parentage)

HeartMassMg ~ Treatment + Sex + TarsusLength + (1 | Parentage)
```

Input:

```text
MorphologyData.csv
```

Outputs:

```text
fig.s7.a.plot.pdf
fig.s7.b.plot.pdf
```

## Supplemental Figure S8: Body-condition morphometrics

This section analyzes nestling body mass, tarsus length, and wing chord. Each trait is modeled as a function of treatment and sex, with parentage included as a random effect.

Models:

```r
FastedBodyMass ~ Treatment + Sex + (1 | Parentage)

TarsusLength ~ Treatment + Sex + (1 | Parentage)

WingChord ~ Treatment + Sex + (1 | Parentage)
```

Input:

```text
MorphologyData.csv
```

Outputs:

```text
fig.s8.a.plot.pdf
fig.s8.b.plot.pdf
fig.s8.c.plot.pdf
```

## Supplemental Figure S9: Growth and feeding rates

This section analyzes daily growth and feeding trajectories for AX and CV nestlings.

The script first removes day 7 observations because they do not represent a full day of feeding data. It then summarizes body mass and formula intake by day and animal ID. Mixed-effects models test whether treatment, day, or their interaction predicts body mass or formula intake over time.

Input:

```text
FeedingGrowthData.csv
```

Outputs:

```text
fig.s9.a.plot.pdf
fig.s9.b.plot.pdf
```

## Supplemental Figure S10: Metabolomics ordinations

This section plots PLS-DA ordinations for four metabolomics modes:

```text
C18 negative mode
C18 positive mode
HILIC negative mode
HILIC positive mode
```

Each panel shows AX and CV samples in the first two PLS-DA components, with 95% ellipses.

Input:

```text
PLSDA.csv
```

Outputs:

```text
fig.s10.a.plot.pdf
fig.s10.b.plot.pdf
fig.s10.c.plot.pdf
fig.s10.d.plot.pdf
```

## Supplemental Figure S11: Metabolomics p-value distributions

This section uses Python code to generate histograms of metabolite-level p-values from linear mixed-effects model results. Metabolites are grouped by analysis mode based on their name prefixes.

Input:

```text
lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
```

Output:

```text
fig.s11.plot.pdf
```

Note: this section is written in Python rather than R and requires:

```python
import re
import pandas as pd
import matplotlib.pyplot as plt
```

If running the full workflow from R, this section should either be moved to a separate Python script or run in an environment that supports both R and Python execution.

## Supplemental Figure S12: Selected metabolite abundances

This section plots treatment differences for selected metabolites, including eicosapentaenoic acid and glyceric acid.

Input:

```text
Metabolomics.csv
```

Outputs:

```text
fig.s.12.a.plot.pdf
fig.s.12.b.plot.pdf
```

## Running the code

Place all input CSV files in the same directory as the script, then run:

```r
source("SuppCode_SuppFigs.R")
```

Because the script uses relative file paths, the working directory should be set to the folder containing the script and input files.

Example:

```r
setwd("path/to/Supplemental_Figures")
source("SuppCode_SuppFigs.R")
```

