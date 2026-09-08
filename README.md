# Single-Cell Multiome (RNA + ATAC) Analysis Template

A commented R template for integrated single-cell RNA-seq + ATAC-seq
("multiome") analysis using **Seurat** and **Signac**, covering the full
workflow from raw 10x Genomics Multiome data to transcription-factor
regulatory analysis:

- Quality control (RNA + ATAC metrics)
- SCTransform (RNA) and TF-IDF/LSI (ATAC) dimensionality reduction
- Weighted Nearest Neighbor (WNN) integration of both modalities
- Clustering and cluster annotation
- Marker gene identification (RNA and SCT assays)
- Gene activity scores from chromatin accessibility
- Transcription factor motif enrichment with chromVAR
- Candidate regulator identification (combining RNA and motif evidence)
- Variant (e.g. GWAS SNP) to regulatory element to transcription factor
  integration, including coverage plots around variants of interest

## Acknowledgements

This template follows the general structure and best practices described in
the official **[Seurat](https://satijalab.org/seurat/)** and
**[Signac](https://stuartlab.org/signac/)** vignettes, in particular the 10x
Genomics Multiome (RNA + ATAC) WNN analysis vignette. It is meant as a
reusable starting point / tutorial, not as a description of any specific
published study — all group names, gene lists, and variant categories in
this script are placeholders to be replaced with your own data.

## Requirements

R (>= 4.2 recommended) with the following packages:

```r
install.packages(c("readxl", "writexl", "ggplot2", "dplyr", "patchwork"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Seurat",
  "Signac",
  "EnsDb.Hsapiens.v86",                 # swap for your species/genome build
  "BSgenome.Hsapiens.UCSC.hg38",        # swap for your species/genome build
  "TxDb.Hsapiens.UCSC.hg38.knownGene",  # swap for your species/genome build
  "Rsamtools",
  "AnnotationHub",
  "chromVAR",
  "JASPAR2020",
  "TFBSTools",
  "motifmatchr",
  "presto"
))

# Optional, used only for one alternative plotting style:
install.packages("SCpubr")
```

## Input data

You need, from 10x Genomics Cell Ranger ARC (or equivalent) output:

- `filtered_feature_bc_matrix.h5` — combined RNA + ATAC filtered count matrix
- `atac_fragments.tsv.gz` (+ its `.tbi` index, created automatically by the
  script if missing) — ATAC fragments file

Optionally, for the downstream variant-integration section:

- One or more variant tables (e.g. exported from a GWAS regulome browser)
  with at minimum: chromosome, start, end, and a variant ID column
- A "regulome" annotation table to merge with the final SNP/TF/accessibility
  summary
- A supplementary table (e.g. from a published paper) listing reference
  marker genes per cell type, if you want to score your clusters against an
  external reference signature

## Configuration

All paths and key parameters are set in a single block at the top of
`multiome_analysis.R` — no need to edit the rest of the script:

```r
BASE_DIR <- "path/to/your/project"
RAW_DATA_DIR <- file.path(BASE_DIR, "raw_data")
...
GENOME_BUILD <- "hg38"
SPECIES_TAXID <- 9606
QC_MIN_RNA_FEATURES <- 200
...
```

Update the QC thresholds after inspecting the violin plots the script
produces — the default values are a starting point, not a universal rule.

## Sections marked "PROJECT-SPECIFIC"

Throughout the script, comments flag sections that encode choices specific
to a particular biological system or study, for example:

- Cluster-to-cell-type name mapping (`RenameIdents`)
- Reference marker gene panels (`reference_markers_group*`, `canonical_markers`)
- Cell type / group labels used throughout (`GroupA`, `GroupB`, ...)
- The two variant tables used in the SNP/TF integration section

Replace these with the labels, genes, and tables relevant to your own
analysis before running the script.

## Usage

The script is meant to be run interactively (e.g. in RStudio), section by
section, since several steps benefit from inspecting intermediate plots
(QC violin plots, elbow plots, UMAPs) before deciding on parameters or
filtering thresholds. The Seurat object is saved to disk (`RDS_FILE`) after
each major step, so you can reload it (`readRDS(RDS_FILE)`) and resume from
any point without recomputing earlier steps.

## License

Distributed under the MIT License — see [LICENSE](LICENSE).
