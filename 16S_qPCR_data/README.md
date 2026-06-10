# 16S qPCR Bacterial DNA Quantification

This folder contains raw and exported results from bacterial 16S qPCR assays run with the Zymo Research Femto™ Bacterial DNA Quantification Kit and analyzed in Applied Biosystems Design & Analysis Software. The assay target in the exported results is `Bacteria`, and the primary analysis output is a well-level CSV table containing Cq values, standard-curve quantities, amplification status, and plate-level standard-curve metrics.

## Files

| File | Description |
|---|---|
| `_e2006_femto_bacterial_dna_quantification_kit_ver.pdf` | Zymo Research Femto™ Bacterial DNA Quantification Kit manual. The kit is designed to quantify bacterial DNA by real-time PCR using bacterial DNA standards and a no-template control. |
| `16S_qPCR_PL1_20230920_115804.eds` through `16S_qPCR_PL7_20230921_100914.eds` | Native Applied Biosystems QuantStudio/Design & Analysis Software experiment files for plates PL1–PL7. |
| `Results_Report_2026-06-09-*.pdf` | Exported Applied Biosystems plate results reports corresponding to the individual `.eds` files. These reports include run metadata, well tables, amplification plots, multicomponent plots, and standard-curve information. |
| `16S_qPCR_Results.csv` | Combined well-level qPCR results table for all plates. This is the main machine-readable results file. |

## Assay and instrument summary

The Femto™ kit is designed for bacterial DNA quantification from purified mixed DNA samples and supports detection across a 20 fg to 20 ng range from as little as 1 µl of sample. The kit uses a real-time qPCR premix and a standard dilution series with a no-template control.

The qPCR runs were performed on an Applied Biosystems QuantStudio™ 7 Pro System using a 96-well 0.2-ml block. The exported reports were generated with Applied Biosystems Design & Analysis Software v2.7.0 using the Primary Analysis and Standard Curve plugins. Quantification cycle calling used the Baseline Threshold method.

## Dataset structure

`16S_qPCR_Results.csv` contains 324 well-level records from seven plates:

| Plate | Number of wells in exported CSV |
|---|---:|
| PL1 | 72 |
| PL2 | 24 |
| PL3 | 69 |
| PL4 | 24 |
| PL5 | 63 |
| PL6 | 24 |
| PL7 | 48 |

The CSV contains the following fields:

```text
Plate, Well, Well Position, Sample, Quantity, Target, Dye, Task, Reporter,
Quencher, Amp Status, Cq, Cq Mean, Cq Confidence, Cq SD, Quantity Mean,
Quantity SD, Auto Threshold, Threshold, Auto Baseline, Baseline Start,
Baseline End, Y-Intercept, R2, Slope, Efficiency, Standard Deviation,
Standard Error
```

Key fields:

- `Plate`: qPCR plate identifier (`PL1`–`PL7`).
- `Well Position`: physical well location on the 96-well plate.
- `Sample`: standard, no-template control, or biological/sample identifier.
- `Task`: qPCR role assigned by the analysis software, including `STANDARD`, `UNKNOWN`, and `NTC`.
- `Cq`: quantification cycle for the individual well.
- `Cq Mean` and `Cq SD`: replicate-level Cq summary statistics where available.
- `Quantity`, `Quantity Mean`, and `Quantity SD`: quantity estimates from the plate-specific standard curve.
- `Amp Status`: amplification call from the analysis software (`AMP`, `INCONCLUSIVE`, or `NO_AMP`).
- `R2`, `Slope`, `Efficiency`, `Y-Intercept`: plate-level standard-curve metrics repeated on each row for that plate.

## Standard-curve performance

Plate-level standard-curve metrics in `16S_qPCR_Results.csv` are:

| Plate | R² | Slope | Efficiency (%) | Y-intercept | SD | SE |
|---|---:|---:|---:|---:|---:|---:|
| PL1 | 0.9977 | -3.5879 | 89.98 | 39.3624 | 0.3654 | 0.0399 |
| PL2 | 0.9992 | -3.6273 | 88.66 | 38.8908 | 0.2102 | 0.0229 |
| PL3 | 0.9972 | -3.7303 | 85.38 | 40.6161 | 0.4163 | 0.0454 |
| PL4 | 0.9985 | -3.6225 | 88.82 | 38.6799 | 0.2936 | 0.0320 |
| PL5 | 0.9991 | -3.5860 | 90.05 | 38.5319 | 0.2216 | 0.0242 |
| PL6 | 0.9990 | -3.6926 | 86.56 | 40.7955 | 0.2438 | 0.0266 |
| PL7 | 0.9986 | -3.7073 | 86.10 | 41.0091 | 0.2939 | 0.0321 |

All seven plates have high linearity by R². Amplification calls in the combined CSV are distributed as 302 `AMP`, 15 `INCONCLUSIVE`, and 7 `NO_AMP` wells.

## Recommended use

Use `16S_qPCR_Results.csv` as the primary input for downstream analysis. The `.eds` files are retained as raw instrument files, and the PDF reports are retained as human-readable audit records for each plate.

A typical analysis workflow is:

1. Load `16S_qPCR_Results.csv`.
2. Filter or flag wells with `Amp Status` equal to `INCONCLUSIVE` or `NO_AMP`.
3. Separate `STANDARD`, `NTC`, and `UNKNOWN` rows using the `Task` column.
4. Summarize unknown samples by `Plate` and `Sample`, using `Cq Mean`, `Quantity Mean`, and `Quantity SD` where appropriate.
5. Confirm that no-template controls behave as expected before interpreting low-quantity samples.
6. Interpret quantities using the units and dilution assumptions defined by the Femto™ kit standard curve and the experimental setup.

## Minimal analysis example

```python
import pandas as pd

qpcr = pd.read_csv("16S_qPCR_Results.csv")

# Plate-level standard curve summary
std_curve = (
    qpcr.groupby("Plate")[["R2", "Slope", "Efficiency", "Y-Intercept"]]
        .first()
        .reset_index()
)

# Unknown-sample summary, excluding non-amplifying and inconclusive wells
unknown_summary = (
    qpcr[(qpcr["Task"] == "UNKNOWN") & (qpcr["Amp Status"] == "AMP")]
        .groupby(["Plate", "Sample"], as_index=False)
        .agg(
            n_wells=("Well", "count"),
            mean_cq=("Cq", "mean"),
            sd_cq=("Cq", "std"),
            mean_quantity=("Quantity", "mean"),
            sd_quantity=("Quantity", "std"),
        )
)
```

## Notes and cautions

- Quantity values are standard-curve-derived estimates. Confirm whether they represent DNA mass, copy number, or another exported unit before comparing them across projects.
- Cq and quantity summaries should be interpreted only after evaluating replicate consistency and NTC behavior.
- Plate-specific standard curves should be retained in downstream analyses because the standard-curve parameters differ slightly among plates.
- The PDF reports and `.eds` files should be archived with the CSV to preserve provenance and enable reanalysis.

