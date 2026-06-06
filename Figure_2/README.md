# Figure 2 Code

This folder contains the code used to generate Figure 2 and supporting statistical analyses for the manuscript:

**Gut microbiota promote early-life digestive function in songbirds**

The scripts in this folder generate the Figure 2 organ-mass panels, metabolomics volcano plot, metabolite inset, bile-acid abundance panels, and supporting mixed-effects model results.

## Scripts

```text
SuppCode_Fig2.R
Figure2_log10_ANCOVA_Fisher_within_ParentagePair_100perm.R
Code_For_Metabolite_Differential_Abundance.R
Code_For_VolcanoPlot_Fig2.py
Code_For_Inset.py
script_for_FDR_correction_for_bile_acids.py
```

## Overview

The Figure 2 workflow has two main components:

1. **Digestive organ morphology analyses**
   - Compares AX and CV nestlings for liver, intestine, pancreas, and related organ-mass traits.
   - Generates the main organ-mass panels for Figure 2.
   - Fits mixed-effects models that account for body size, sex, and parentage structure.

2. **Metabolomics analyses**
   - Fits metabolite-level mixed-effects models comparing AX and CV nestlings while adjusting for sex and parentage.
   - Generates a metabolomics volcano plot.
   - Performs Benjamini-Hochberg FDR correction for bile-acid p-values within annotation category and analysis mode.
   - Produces selected metabolite abundance plots and a paired inset plot.

The scripts are written primarily in R, with Python used for the volcano plot and inset paired t-test.

## Recommended run order

Run the scripts in this order if reproducing the full Figure 2 workflow from raw or intermediate input files:

```text
1. Code_For_Metabolite_Differential_Abundance.R
2. Code_For_VolcanoPlot_Fig2.py
3. script_for_FDR_correction_for_bile_acids.py
4. Code_For_Inset.py
5. Figure2_log10_ANCOVA_Fisher_within_ParentagePair_100perm.R
6. SuppCode_Fig2.R
```

`SuppCode_Fig2.R` is the main figure-generation script. Some of its inputs are outputs or manually prepared summaries from the other scripts.

## Required R packages

The R scripts use the following packages:

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
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(nlme)
library(purrr)
library(tibble)
library(ggplot2)
```

Some scripts load overlapping sets of packages. The metabolite differential abundance script will attempt to install missing R packages from CRAN.

## Required Python packages

The Python scripts use:

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
```

## Input files

The scripts expect the following input files, depending on which part of the workflow is being run:

```text
MorphologyData.csv
DigestiveOrgans.csv
MetaboliteVolcano.csv
Fig2E_INSET.csv
Metabolomics.csv
merged_metabolomics.csv
SuppTableS5.xlsx
```

The metabolite differential abundance script expects metabolite abundances in `merged_metabolomics.csv` and sample metadata in `SuppTableS5.xlsx`.

## Output files

The main figure-generation script writes the following Figure 2 panel PDFs:

```text
fig.2.a.plot.pdf
fig.2.b.plot.pdf
fig.2.c.plot.pdf
fig.2.d.plot.pdf
fig.2.x.plot.pdf
fig.2.e.plot.pdf
fig.2.e.inset.pdf
fig.2.f.plot.pdf
fig.2.g.plot.pdf
fig.2.h.plot.pdf
fig.2.i.plot.pdf
```

The organ-mass ANCOVA and Fisher-combined-test script writes model summaries, permutation results, plotting data, and PDFs with filenames beginning:

```text
Figure2_organ_mass_LOG10_ANCOVA_lmer_Treatment_by_Sex_with_SIplusLI_no_clutch_color_no_lines
```

The metabolite differential abundance script writes:

```text
lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
```

The volcano-plot script writes:

```text
volcano_plot_all_group_plus_sex_fixed_mean_difference.pdf
```

## Script descriptions

## `SuppCode_Fig2.R`

This is the main Figure 2 plotting script. It imports organ morphology data, metabolomics summary data, and selected metabolite abundance data, then generates the primary Figure 2 panels.

The first section uses `MorphologyData.csv` to plot organ-mass comparisons between AX and CV nestlings. The script removes the outlier individual `NL113` before generating the liver, large intestine, small intestine, and pancreas panels. It also includes an additional intestinal-mass plot.

Generated panels include:

```text
fig.2.a.plot.pdf    Liver mass
fig.2.b.plot.pdf    Large intestine mass
fig.2.c.plot.pdf    Small intestine mass
fig.2.d.plot.pdf    Pancreas mass
fig.2.x.plot.pdf    Total intestinal mass
```

The second section uses `MetaboliteVolcano.csv` to generate the Figure 2 metabolomics volcano plot. The volcano plot shows log2 fold change on the x-axis and -log10 p-value on the y-axis, with significant metabolites colored according to direction of change.

Output:

```text
fig.2.e.plot.pdf
```

The third section uses `Fig2E_INSET.csv` to generate a paired inset plot for Figure 2E. The input is reshaped from wide to long format, plotted as paired AX and CV values, and saved as a PDF.

Output:

```text
fig.2.e.inset.pdf
```

The final section uses `Metabolomics.csv` to plot selected bile-acid or metabolite abundances, including cholic acid, beta-muricholic acid, taurochenodeoxycholic acid, and glycocholic acid.

Outputs:

```text
fig.2.f.plot.pdf
fig.2.g.plot.pdf
fig.2.h.plot.pdf
fig.2.i.plot.pdf
```

## `Figure2_log10_ANCOVA_Fisher_within_ParentagePair_100perm.R`

This script performs the main organ-mass statistical analysis for Figure 2 using mixed-effects ANCOVA models.

For each organ, the script fits the model:

```r
log10_OrganMass ~ log10_TarsusLength + Treatment * Sex + (1 | ParentagePair)
```

The model tests whether treatment predicts organ mass while accounting for body size, sex, treatment-by-sex interaction, and parentage-pair structure.

The script analyzes the following organ traits:

```text
LiverMassMg
EmptySIMassMg
EmptyLIMassMg
EmptySIplusLIMassMg
PancreasMassMg
EmptyGizzardMassMg
HeartMassMg
```

The script also calculates a combined small-plus-large intestine mass variable:

```text
EmptySIplusLIMassMg
```

Individual `NL113` is removed before model fitting.

The script writes Type III ANOVA results for each organ and generates a raw organ-mass violin plot with model-derived treatment significance annotations.

It also performs a permutation-calibrated Fisher combined test across five digestive organ traits:

```text
LiverMassMg
EmptySIMassMg
EmptyLIMassMg
PancreasMassMg
EmptyGizzardMassMg
```

For this aggregate test, treatment labels are permuted within `ParentagePair`, preserving parentage structure and within-pair treatment counts. The observed Fisher statistic is compared against a null distribution generated from 100 permutations.

Key outputs include:

```text
*_type3_anova.csv
*_observed_Treatment_pvalues_5_digestive_organs_for_permutation_Fisher.csv
*_observed_Fisher_statistic_Treatment_5_digestive_organs.csv
*_Permutation_null_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.csv
*_Permutation_summary_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.csv
*_Permutation_null_Fisher_statistic_Treatment_5_digestive_organs_within_ParentagePair.pdf
*_raw_organ_mass_plot_data.csv
*_raw_organ_mass_violinplot_custom_two_row_layout_TREATMENT_BRACKET_NO_TOPRIGHT_stats_from_model.pdf
```

## `Code_For_Metabolite_Differential_Abundance.R`

This script fits metabolite-level linear mixed-effects models.

Inputs:

```text
merged_metabolomics.csv
SuppTableS5.xlsx
```

The metabolomics file should contain one metabolite per row, with the first column giving the metabolite ID or name and the remaining columns giving sample-level abundance values.

The metadata file should contain sample IDs in the first column and include columns identifying parentage, sex, and treatment group. The script attempts to detect columns corresponding to group, sex, and parentage.

For each metabolite, the script fits:

```r
value ~ group + sex
random = ~ 1 | parentage
```

Here, `AX` is used as the reference group, so the group coefficient estimates CV relative to AX after adjustment for sex and parentage.

The script also reports unadjusted group means and the mean difference:

```text
mean CV - mean AX
```

Output:

```text
lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
```

This output is used by the volcano-plot script.

## `Code_For_VolcanoPlot_Fig2.py`

This Python script generates a volcano plot from the metabolite-level mixed-effects model results.

Input:

```text
lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
```

The x-axis is the unadjusted mean difference:

```text
mean_difference_CV_minus_AX
```

The y-axis is:

```text
-log10(p_value_group)
```

Metabolites with group p-values below 0.05 are plotted as significant. The script draws a horizontal threshold line at p = 0.05 and a vertical line at zero mean difference.

Output:

```text
volcano_plot_all_group_plus_sex_fixed_mean_difference.pdf
```

## `script_for_FDR_correction_for_bile_acids.py`

This Python script performs Benjamini-Hochberg false-discovery-rate correction for bile-acid metabolite p-values from the linear mixed-effects models comparing CV and AX birds.

The script applies FDR correction separately within bile-acid annotation sets and analysis modes. In the provided version, two bile-acid p-value vectors are defined directly in the script:

```text
C18_NEG_BileAcids
HILIC_NEG_BileAcids
```

The script uses:

```python
multipletests(pvals, method="fdr_bh")
```

to return BH-adjusted q-values in the original p-value order.

This correction is intended to correspond to the bile-acid subset of the metabolomics results, with uncorrected p-values reported in Table S7.

Inputs:

```text
Hard-coded bile-acid p-value vectors from Table S7
```

Output:

```text
Printed BH-adjusted q-values for each bile-acid set
```

Note: in the uploaded version, the p-value vectors are named `C18_NEG_BileAcids` and `HILIC_NEG_BileAcids`, but the correction calls use `set1` and `set2`. To run the script as written, define:

```python
set1 = C18_NEG_BileAcids
set2 = HILIC_NEG_BileAcids
```

before calling `bh_correct()`, or replace `set1` and `set2` directly with the named arrays.

## `Code_For_Inset.py`

This Python script performs a paired t-test for the inset associated with Figure 2E.

The script defines paired AX and CV vectors directly in the code, then runs:

```python
stats.ttest_rel(x, y)
```

The included result is:

```text
Mean difference (x - y) = 0.5153748694
SD of differences = 0.1726849
n = 5
df = 4
paired t-test t-statistic = 6.6709
two-sided p-value = 0.00263
```

This statistical result corresponds to the paired comparison shown in the Figure 2E inset.

## Running the code

Set the working directory to the `Figure_2` folder before running the scripts.

Example R workflow:

```r
setwd("path/to/Figure_2")

source("Code_For_Metabolite_Differential_Abundance.R")
source("Figure2_log10_ANCOVA_Fisher_within_ParentagePair_100perm.R")
source("SuppCode_Fig2.R")
```

Example Python workflow:

```bash
python Code_For_VolcanoPlot_Fig2.py
python script_for_FDR_correction_for_bile_acids.py
python Code_For_Inset.py
```

## Notes

The Figure 2 scripts are intended for figure reproduction and supporting analysis rather than for use as a generalized analysis package. Several plotting parameters, axis limits, colors, significance labels, and output filenames are manually specified to match the final manuscript figure layout.

The organ-mass plotting and organ-mass statistical scripts use related but not identical input files. `SuppCode_Fig2.R` uses `MorphologyData.csv` for the main panel plotting, while `Figure2_log10_ANCOVA_Fisher_within_ParentagePair_100perm.R` uses `DigestiveOrgans.csv` for the log10 ANCOVA and permutation-based aggregate test.

The Python volcano-plot script depends on the CSV produced by the metabolite differential abundance R script. The bile-acid FDR correction script currently uses hard-coded p-value vectors and should be checked against the final Table S7 values before release. The Figure 2E inset plotting in `SuppCode_Fig2.R` depends on `Fig2E_INSET.csv`, while `Code_For_Inset.py` documents and reproduces the paired t-test for the inset comparison.
