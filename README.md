# Gut microbiota promote early-life digestive function in songbirds

This repository contains the supplemental R code and input data used to reproduce Figure 1 and Figure 2 from the manuscript:

**Gut microbiota promote early-life digestive function in songbirds**  
Brian K. Trevelline et al.

This repo also includes all data files used to generate supplemental figures. 

## Repository structure

```text
.
├── README.md
├── QIIME_code/
├── 16S_qPCR_data/
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
|   ├── Code_For_Inset.py
|   ├── Code_For_Metabolite_Differential_Abundance.R
|   ├── Code_For_VolcanoPlot_Fig2.py
|   ├── script_for_FDR_correction_for_bile_acids.R
|   ├── SuppCode_Fig2.R
|   ├── data/
│   │   ├── data_normalized_C18_NEG.csv
│   │   ├── data_normalized_C18_POS.csv
│   │   ├── data_normalized_HILIC_NEG.csv
│   │   ├── data_normalized_HILIC_POS.csv
│   │   ├── SuppTableS4.csv
│   │   ├── DigestiveOrgans.csv
│   │   ├── Fig2E_INSET.csv
│   │   ├── MetaboliteVolcano.csv
│   │   ├── Metabolomics.csv
│   │   ├── MorphologyData.csv
│   │   ├── lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
│   │   └── merged_metabolomics.csv
|   └── output/
|       ├── Figure2_organ_mass_LOG10_ANCOVA_lmer_Treatment_by_Sex_with_SIplusLI_no_clutch_color_no_lines_type3_anova.csv
|       ├── fig.2.a.plot.pdf
|       ├── fig.2.b.plot.pdf
|       ├── fig.2.c.plot.pdf
|       ├── fig.2.d.plot.pdf
|       ├── fig.2.e.inset.pdf
|       ├── fig.2.f.plot.pdf
|       ├── fig.2.g.plot.pdf
|       ├── fig.2.h.plot.pdf
|       ├── fig.2.i.plot.pdf
|       ├── fig.2.x.plot.pdf
|       └── fig.2.e.plot
├── Supplemental_Figures/
    ├── SuppCode_SuppFigs.R
    ├── data/
    │   ├── ASVRelativeAbundances_Genus.csv
    │   ├── FeedingGrowthData.csv
    │   ├── INOC-CV_ASVs.csv
    │   ├── Metabolomics.csv
    │   ├── MorphologyData.csv
    │   ├── NonTransferredASVs.csv
    │   ├── PLSDA.csv
    │   ├── RPCAdistancesFIGS2.csv
    │   ├── RPCAdistancesFIGS3.csv
    │   ├── RPCAordinationFIGS1.csv
    |   └── lme_results_all_group_plus_sex_fixed_with_mean_difference.csv
    └── output/
        ├── fig.s.12.a.plot.pdf
        ├── fig.s.12.b.plot.pdf
        ├── fig.s1.plot.pdf
        ├── fig.s2.a.plot.pdf
        ├── fig.s2.b.plot.pdf
        ├── fig.s3.plot.pdf
        ├── fig.s4.a.plot.pdf
        ├── fig.s4.b.plot.pdf
        ├── fig.s6.plot.pdf
        ├── fig.s7.a.plot.pdf
        ├── fig.s7.b.plot.pdf
        ├── fig.s8.a.plot.pdf
        ├── fig.s8.b.plot.pdf
        ├── fig.s8.c.plot.pdf
        ├── fig.s9.a.plot.pdf
        ├── fig.s9.b.plot.pdf
        ├── fig.s10.a.plot.pdf
        ├── fig.s10.b.plot.pdf
        ├── fig.s10.c.plot.pdf
        ├── fig.s10.d.plot.pdf
        └── fig.s11.plot.pdf
```


Refer to the manuscript methods for the full experimental design and treatment descriptions.



