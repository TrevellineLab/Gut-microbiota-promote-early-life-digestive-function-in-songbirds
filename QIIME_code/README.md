# QIIME2 microbiome analysis inputs and exported results

This directory contains the QIIME2 supplemental analysis code, sample metadata, contaminant/feature exclusion lists, and exported tables used to generate microbiome analyses for **“Gut microbiota promote early-life digestive function in songbirds.”** The workflow processes paired-end 16S rRNA amplicon reads, denoises reads with DADA2, assigns taxonomy with a SILVA 138 classifier trained for the 515F-Y/926R primer region, removes contaminant feature sets, and generates alpha- and beta-diversity outputs for Figure 1A/1B and Figure S1.
The classifier file is provided in the Zenodo link referenced in the manuscript. 

## File overview

### Main workflow

| File | Description |
|---|---|
| `SuppCode_QIIME2` | Supplemental QIIME2 command script. It imports paired-end reads, performs DADA2 denoising, filters ASVs by length and taxonomy, removes contaminant ASVs, runs Gemelli RPCA, performs PERMANOVA, exports ordinations and distance matrices, and exports selected feature tables with taxonomy. |
| `sample-manifest.csv` | QIIME2 paired-end manifest. Each sample has forward and reverse FASTQ paths. |
| `sample-metadata.txt` | Full sample metadata table with `SampleID`, `IndividualID`, `Treatment`, and `TreatmentDay`. |
| `sample-metadata-collapsed.txt` | Metadata after grouping samples by individual ID for Figure 1A/Figure 1B analyses. |
| `sample-metadata-collapsed-TreatmentDay.txt` | Minimal metadata for treatment-day-collapsed summaries and taxa barplots. |

### ASV exclusion lists

These files contain feature IDs used to remove select ASVs (e.g., contaminants) from the feature table.

| File | Contents | Number of ASVs |
|---|---:|---:|
| `ctrl-asvs(2).txt` | ASVs detected in environmental/control samples. | 72 |
| `food-asvs(2).txt` | ASVs detected in food/formula samples. | 33 |
| `mec-asvs(2).txt` | ASVs detected in meconium samples. | 91 |
| `blank-asvs(2).txt` | ASVs detected in blank samples. | 36 |

### Exported Figure 1 files

| File | Description |
|---|---|
| `Fig1A_distance-matrix.tsv` | Exported RPCA distance matrix for the Figure 1A beta-diversity analysis after filtering control, food, and meconium ASVs; grouping by individual; excluding PR samples; and running Gemelli RPCA. The matrix contains 50 samples. |
| `Fig1A_ordination.txt` | Exported RPCA ordination for Figure 1A. |
| `filtered-table-no-ctrls-food-blank-CV-1-2-INOC-transposed-with-taxonomy.tsv` | Transposed feature table with taxonomy for the comparison of `CV-1-2` and `INOC` samples. Contains 244 feature rows and sample abundance columns plus taxonomy/confidence fields. |
| `filtered-table-no-ctrls-food-blank-AX-1-2-MEC-transposed-with-taxonomy.tsv` | Transposed feature table with taxonomy for the comparison of `AX-1-2` and `MEC` samples. Contains 276 feature rows and sample abundance columns plus taxonomy/confidence fields. |

### Exported Figure S1 files

| File | Description |
|---|---|
| `FigS1_distance-matrix.tsv` | Exported RPCA distance matrix for the Figure S1 beta-diversity analysis using the filtered table after removal of control and food ASVs. The matrix contains 159 samples. |
| `FigS1_ordination.txt` | Exported RPCA ordination for Figure S1. |

## Analysis workflow

The workflow in `SuppCode_QIIME2` proceeds through the following major steps.

1. **Import paired-end reads**
   - Uses `sample-manifest.csv` as a QIIME2 paired-end manifest.
   - Imports reads as `SampleData[PairedEndSequencesWithQuality]`.
   - Summarizes demultiplexed reads with `qiime demux summarize`.

2. **Denoise reads with DADA2**
   - Runs `qiime dada2 denoise-paired`.
   - Forward reads are truncated at 230 bp.
   - Reverse reads are truncated at 160 bp.
   - Outputs the ASV table, representative sequences, and denoising statistics.

3. **Filter representative sequences by expected amplicon length**
   - Retains sequences between 368 and 376 nucleotides, corresponding to the expected 372-nt amplicon length ±4 nt.
   - Separately exports the deviating representative sequences so they can be removed from the feature table.

4. **Build phylogenetic tree**
   - Aligns filtered representative sequences with MAFFT.
   - Masks the alignment.
   - Builds an unrooted FastTree tree.
   - Midpoint-roots the tree.

5. **Train and apply taxonomy classifier**
   - Extracts SILVA 138 reference reads using primers 515F-Y and 926R.
   - Trains a naive Bayes classifier.
   - Classifies filtered representative sequences with `classify-sklearn`.

6. **Filter non-target and contaminant features**
   - Removes length-deviating ASVs.
   - Removes archaeal, chloroplast, mitochondrial, and eukaryotic features.
   - Removes ASVs identified in control samples using `ctrl-asvs.txt`.
   - Removes ASVs identified in food/formula samples using `food-asvs.txt`.
   - For specific downstream analyses, additionally removes ASVs identified in meconium or blank samples using `mec-asvs.txt` or `blank-asvs.txt`.

7. **Generate Figure S1 beta-diversity output**
   - Runs Gemelli RPCA on `filtered-table-no-ctrls-food.qza`.
   - Generates an Emperor biplot.
   - Runs PERMANOVA by `TreatmentDay` with pairwise tests.
   - Exports the RPCA ordination and distance matrix.

8. **Generate Figure 1A/Figure 1B outputs**
   - Removes meconium ASVs.
   - Calculates Chao1 alpha diversity.
   - Groups samples by `IndividualID` using `mean-ceiling`.
   - Excludes PR samples for the grouped beta-diversity analysis.
   - Runs Gemelli RPCA.
   - Runs PERMANOVA by `Treatment` with pairwise tests.
   - Exports the RPCA ordination and distance matrix.

9. **Generate grouped taxa and feature-abundance tables**
   - Removes blank ASVs.
   - Filters samples to `MEC`, `AX-1-2`, `CV-1-2`, `INOC`, and `PR-SICONT`.
   - Groups by `TreatmentDay` for taxa barplots.
   - Exports transposed feature tables with taxonomy for:
     - `CV-1-2` versus `INOC`
     - `AX-1-2` versus `MEC`

## Software and reference data

The supplemental code activates the following QIIME2 environment:

```bash
source activate qiime2-2023.5
```

The workflow also assumes that the following reference files are available in the working directory or in paths visible to QIIME2:

```text
silva-138-99-seqs.qza
silva-138-99-tax.qza
```

The trained classifier is generated within the workflow as:

```text
silva-138-99-515FY-926R-nb-classifier.qza
```

## Expected input structure

The manifest uses `$PWD/MC-WY-336/...fastq.gz` paths. To rerun the workflow without editing the manifest, place the raw FASTQ files in a directory named `MC-WY-336` inside the working directory, or update `sample-manifest.csv` to point to the correct absolute or relative paths.

Expected input files before running the QIIME2 workflow:

```text
sample-manifest.csv
sample-metadata.txt
sample-metadata-collapsed.txt
sample-metadata-collapsed-TreatmentDay.txt
ctrl-asvs.txt
food-asvs.txt
mec-asvs.txt
blank-asvs.txt
silva-138-99-seqs.qza
silva-138-99-tax.qza
MC-WY-336/*.fastq.gz
```

## Reproducing the analysis

From a shell with QIIME2 installed, run the supplemental code file from the directory containing the manifest, metadata, ASV lists, reference files, and raw FASTQ directory:

```bash
bash SuppCode_QIIME2
```

## Outputs

The workflow produces QIIME2 artifacts (`.qza`), visualizations (`.qzv`), exported ordinations, exported distance matrices, and exported feature tables. The uploaded exported results correspond to:

```text
Fig1A_distance-matrix.tsv
Fig1A_ordination.txt
FigS1_distance-matrix.tsv
FigS1_ordination.txt
filtered-table-no-ctrls-food-blank-CV-1-2-INOC-transposed-with-taxonomy.tsv
filtered-table-no-ctrls-food-blank-AX-1-2-MEC-transposed-with-taxonomy.tsv
```

