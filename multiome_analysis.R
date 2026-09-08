# =============================================================================
# INTEGRATED scRNA-seq + scATAC-seq ANALYSIS WITH WEIGHTED NEAREST NEIGHBOR (WNN)
# =============================================================================
#
# NOTE ON DATA / PUBLICATION STATUS:
# This script was developed on an unpublished dataset. Before making this
# repository (or parts of it) public, review the sections marked
# "PROJECT-SPECIFIC" below: cluster labels, marker gene panels, and the
# GWAS/variant table names may reveal details of an unpublished study.
# Consider keeping the repository private until the associated work is
# published, or replacing project-specific identifiers with placeholders.
#
# =============================================================================
# USER CONFIGURATION — edit these paths and parameters for your project
# =============================================================================

# Root working directory for this analysis (all other paths are derived from it)
BASE_DIR <- "path/to/your/project"

# Folder containing the raw 10x multiome data (filtered_feature_bc_matrix.h5
# and the ATAC fragments file)
RAW_DATA_DIR <- file.path(BASE_DIR, "raw_data")
H5_FILE <- file.path(RAW_DATA_DIR, "filtered_feature_bc_matrix.h5")
FRAGMENTS_FILE <- file.path(RAW_DATA_DIR, "atac_fragments.tsv.gz")

# Folder containing GWAS/regulome variant tables used in the downstream
# SNP-TF-accessibility integration (PROJECT-SPECIFIC section, see below)
REGULOME_DIR <- file.path(BASE_DIR, "regulome")
SNP_TABLE_GROUP_A_FILE <- file.path(REGULOME_DIR, "groupA_regulome_hg19.xlsx")
SNP_TABLE_GROUP_B_FILE <- file.path(REGULOME_DIR, "groupB_regulome_hg19.xlsx")
REGULOME_COMBINED_FILE <- file.path(REGULOME_DIR, "regulome_combined.tsv")

# Supplementary table with published marker genes used to build reference
# signatures (e.g. a supplementary table from a reference paper)
SUPPLEMENTARY_TABLE_FILE <- file.path(BASE_DIR, "supplementary_table.xlsx")

# Output folder (Seurat object, marker tables, plots)
OUTPUT_DIR <- BASE_DIR
RDS_FILE <- file.path(OUTPUT_DIR, "multiome_obj.rds")
COVERAGE_GROUP_A_DIR <- file.path(OUTPUT_DIR, "coverage_groupA")
COVERAGE_GROUP_B_DIR <- file.path(OUTPUT_DIR, "coverage_groupB")

# Genome build and species (used for genome annotation, motif matching, chromVAR)
GENOME_BUILD <- "hg38"
SPECIES_TAXID <- 9606   # NCBI taxonomy ID, 9606 = Homo sapiens

# ENCODE blacklist AnnotationHub record ID for your genome build.
# Run `query(AnnotationHub(), c("blacklist", GENOME_BUILD))` interactively to
# find the correct ID for your genome/AnnotationHub version.
BLACKLIST_AH_ID <- "AH107307"

# QC thresholds (adjust after inspecting the violin plots produced below)
QC_MIN_RNA_FEATURES <- 200
QC_MAX_RNA_FEATURES <- 4500
QC_MAX_PERCENT_MT <- 20
QC_MIN_ATAC_COUNTS <- 1000
QC_MAX_ATAC_COUNTS <- 25000
QC_MIN_TSS_ENRICHMENT <- 3
QC_MAX_NUCLEOSOME_SIGNAL <- 2

# =============================================================================
# LIBRARIES
# =============================================================================
library(readxl)
library(writexl)
library(Seurat)        # Multi-omic analysis and WNN
library(Signac)        # ATAC data handling
library(ggplot2)       # Visualization
library(dplyr)         # Data manipulation
library(EnsDb.Hsapiens.v86)              # Genome annotations (adjust for your species/build)
library(BSgenome.Hsapiens.UCSC.hg38)     # Reference genome (adjust for your species/build)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Transcriptome, used for TSS enrichment
library(Rsamtools)     # To index the fragments file
library(patchwork)     # To combine plots
library(AnnotationHub) # Used to fetch the blacklist

# For post-accessibility analysis
# If you don't have BiocManager, install it first
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Install all required packages (some from CRAN, some from Bioconductor) — run once
#  BiocManager::install(c(
#    "chromVAR",           # motif accessibility analysis
#    "JASPAR2020",         # JASPAR motif database
#    "TFBSTools",          # TFBS analysis tools
#    "motifmatchr",        # motif matching against the genome
#    "BSgenome.Hsapiens.UCSC.hg38",  # hg38 reference genome
#    "presto"              # fast differential testing (wilcoxauc)
#    ))

library(chromVAR)     # scATAC-seq motif accessibility analysis
library(JASPAR2020)   # JASPAR motif models
library(TFBSTools)    # TFBS analysis
library(motifmatchr)  # motif matching
library(BSgenome.Hsapiens.UCSC.hg38)   # required by chromVAR
library(presto)       # fast differential expression testing

# =============================================================================
# DATA LOADING AND OBJECT CREATION
# =============================================================================

# Load counts from a 10x file (HDF5 format)
h5_data <- Read10X_h5(H5_FILE)

# Extract RNA and ATAC count matrices
rna_counts <- h5_data$`Gene Expression`  # RNA matrix
atac_counts <- h5_data$Peaks             # ATAC matrix (peaks)

# Keep peaks on standard chromosomes only
standard_chrs <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)   # standard chromosomes for the reference genome (chr1-22, X, Y, M)
peaks_gr <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))  # convert row names (strings) to GRanges so they can be filtered
keep_peaks <- as.vector(seqnames(peaks_gr) %in% standard_chrs)   # logical vector to keep only peaks on standard chromosomes
atac_counts <- atac_counts[keep_peaks, ]   # filter the ATAC count matrix

# Create the multi-omic Seurat object (RNA + ATAC)
multiome_obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = "multiome_project")
# NOTE: no extra filtering parameters are passed here, otherwise cell-count mismatches can occur

# Index the fragments file (must be in the same folder as the .tsv.gz)
if (!file.exists(paste0(FRAGMENTS_FILE, ".tbi"))) {
  indexTabix(FRAGMENTS_FILE, format = "bed")
}

chrom_assay <- CreateChromatinAssay(   # create the chromatin assay
  counts = atac_counts,       # peak count matrix (rows = peaks, columns = cells)
  sep = c(":", "-"),          # separator used in peak names
  genome = GENOME_BUILD,      # reference genome
  fragments = FRAGMENTS_FILE, # link the (indexed) fragments file to the assay, enabling QC metrics
  min.cells = 1,               # no peak/cell filtering at this stage
  min.features = 0
)
multiome_obj[["ATAC"]] <- chrom_assay  # add the ChromatinAssay (Signac) as "ATAC"

multiome_obj  # check that Assays: RNA, ATAC are present

# Add genome annotations to the ATAC assay
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)  # extract gene annotations as GRanges
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- GENOME_BUILD  # explicitly set the genome metadata
Annotation(multiome_obj[["ATAC"]]) <- annotations  # assign the annotations to the ChromatinAssay

# =============================================================================
# QC (RNA + ATAC)
# =============================================================================
multiome_obj[["percent.mt"]] <- PercentageFeatureSet(multiome_obj, assay = "RNA", pattern = "^MT-")   # percentage of mitochondrial genes
multiome_obj <- NucleosomeSignal(multiome_obj, assay = "ATAC")  # nucleosome signal, based on fragment periodicity. High values (> 2-4) indicate abnormal fragmentation; low values are desirable.
multiome_obj <- TSSEnrichment(multiome_obj, assay = "ATAC")     # high TSS enrichment (> 2-4) indicates good sample prep and sequencing (fragments cut at the right nucleosome positions, good ATAC library quality).
                                                                 # Low TSS enrichment (< 2) signals problems: e.g. excess background noise, non-specific fragmentation, or DNA degradation.
ah <- AnnotationHub()  # search for the blacklist matching your genome build
blacklist_query <- query(ah, c("blacklist", GENOME_BUILD))  # the record ID can change over time, it's safer to use query()
print(blacklist_query)  # inspect the results to pick the correct record
blacklist_gr <- blacklist_query[[BLACKLIST_AH_ID]]  # double-check the ID before assigning it to a variable
multiome_obj$blacklist_ratio <- FractionCountsInRegion(object = multiome_obj, assay = 'ATAC', regions = blacklist_gr)

# =============================================================================
# UNIFIED QC FILTERING
# =============================================================================
# First inspect the QC metric distributions
VlnPlot(multiome_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)   # RNA QC before filtering
VlnPlot(multiome_obj, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "blacklist_ratio"), ncol = 4, pt.size = 0)  # ATAC QC before filtering

multiome_obj <- subset(x = multiome_obj,  # thresholds set based on the violin plots above (PROJECT-SPECIFIC: adjust to your own data)
  subset = nFeature_RNA > QC_MIN_RNA_FEATURES & nFeature_RNA < QC_MAX_RNA_FEATURES & percent.mt < QC_MAX_PERCENT_MT &
           nCount_ATAC > QC_MIN_ATAC_COUNTS & nCount_ATAC < QC_MAX_ATAC_COUNTS &
           TSS.enrichment > QC_MIN_TSS_ENRICHMENT & nucleosome_signal < QC_MAX_NUCLEOSOME_SIGNAL
)

VlnPlot(multiome_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)   # RNA QC after filtering
VlnPlot(multiome_obj, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "blacklist_ratio"), ncol = 4, pt.size = 0)  # ATAC QC after filtering

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# RNA ASSAY ANALYSIS
# =============================================================================
DefaultAssay(multiome_obj) <- "RNA"
multiome_obj <- SCTransform(multiome_obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:10, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# WNN workflows use SCTransform, which replaces the NormalizeData -> FindVariableFeatures -> ScaleData
# workflow with a negative-binomial regularization model that accounts for technical noise
# (heteroscedasticity). It is recommended for multi-omic data because it makes RNA signal sharper
# and more comparable across experiments. We then continue with PCA and UMAP as usual.

VizDimLoadings(multiome_obj, dims = 1:2, reduction = "pca")  # inspect component loadings
ElbowPlot(multiome_obj, ndims = 10)   # inspect variance explained by each component
DimPlot(multiome_obj, reduction = "umap.rna", label = TRUE, repel = TRUE, label.box = TRUE) + ggtitle("RNA clusters")

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# ATAC ASSAY ANALYSIS
# =============================================================================
DefaultAssay(multiome_obj) <- "ATAC"
multiome_obj <- RunTFIDF(multiome_obj)  # TF-IDF (Term Frequency - Inverse Document Frequency) normalization of peak counts. Balances a peak's frequency within a cell against its rarity across the whole dataset, giving more weight to informative peaks.
multiome_obj <- FindTopFeatures(multiome_obj, min.cutoff = "q0")    # select the most variable peaks. 'q0' means use all peaks with variance > 0 (i.e. not constant)
multiome_obj <- RunSVD(multiome_obj)  # Singular Value Decomposition (SVD) on the TF-IDF matrix; the resulting components are called LSI (Latent Semantic Indexing), the ATAC analog of PCA. The first components capture the main sources of biological and technical variability.
multiome_obj <- RunUMAP(multiome_obj, reduction = 'lsi', dims = 2:10, reduction.name = "umap.atac", reduction.key = "atacUMAP_")   # UMAP using LSI components 2-10 (component 1 is excluded because it often correlates with sequencing depth)
DimPlot(multiome_obj, reduction = "umap.atac", label = TRUE, repel = TRUE, label.box = TRUE) + ggtitle("ATAC clusters")

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# MULTI-MODAL INTEGRATION (Weighted Nearest Neighbor Analysis)
# =============================================================================
multiome_obj <- FindMultiModalNeighbors(multiome_obj, reduction.list = list("pca", "lsi"), dims.list = list(1:10, 2:10)) # build the WNN graph using both reductions, so each cell's nearest neighbors combine the weight of both modalities. LSI component 1 is excluded for ATAC.
multiome_obj <- RunUMAP(multiome_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")   # build a UMAP based on the WNN graph

# WNN clustering
multiome_obj <- FindClusters(multiome_obj, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE) # identify cell clusters via shared nearest neighbor (SNN); first computes k-nearest neighbors, then builds the SNN graph
DimPlot(multiome_obj, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.box = TRUE) + ggtitle("WNN clusters")

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# PROJECT-SPECIFIC: comparison against reference signatures (AddModuleScore)
# =============================================================================
# Compares gene sets from a supplementary table (e.g. from a reference paper)
# against your own clusters, using AddModuleScore. Adjust `cell_type_column`
# and the group names below to match your own reference table and cell types.
supplementary <- read_excel(SUPPLEMENTARY_TABLE_FILE, sheet = "Supplementary Table 1", skip = 2)

cell_type_column <- "Cell type"   # PROJECT-SPECIFIC: column name in the supplementary table identifying the cell type/group
reference_cell_types <- c("GroupA", "GroupB", "GroupC", "GroupD")  # PROJECT-SPECIFIC: replace with your own group labels
top_n_genes <- 100

top_genes_by_group <- list()
for (grp in reference_cell_types) {
  tab <- subset(supplementary, .data[[cell_type_column]] == grp)
  tab <- tab[order(tab$avg_log2FC, decreasing = TRUE), ]
  top_genes_by_group[[grp]] <- head(tab$gene, top_n_genes)
}

# AddModuleScore requires a minimum number of genes actually present in the dataset
top_genes_filtered <- lapply(top_genes_by_group, function(g) intersect(g, rownames(multiome_obj[["RNA"]])))

multiome_obj <- NormalizeData(multiome_obj, assay = "RNA")
for (grp in names(top_genes_filtered)) {
  multiome_obj <- AddModuleScore(multiome_obj, features = list(top_genes_filtered[[grp]]), name = paste0("Top", top_n_genes, "_", grp), assay = "RNA")
}

vln_plots <- lapply(names(top_genes_filtered), function(grp) {
  VlnPlot(multiome_obj, features = paste0("Top", top_n_genes, "_", grp, "1"), group.by = "seurat_clusters", pt.size = 0) + ggtitle(paste(grp, "signature"))
})
wrap_plots(vln_plots)  # combine all signature violin plots into a single figure

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# MARKER IDENTIFICATION ON THE WNN CLUSTERS (independent of the reference table)
# =============================================================================
DefaultAssay(multiome_obj) <- "SCT"
multiome_obj_markers_SCT <- FindAllMarkers(multiome_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
multiome_obj_markers_SCT %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)    # sort by avg_log2FC
multiome_obj_markers_SCT <- subset(multiome_obj_markers_SCT, subset = multiome_obj_markers_SCT$p_val_adj < 0.05)

saveRDS(multiome_obj, file = RDS_FILE)
write.csv(multiome_obj_markers_SCT, file.path(OUTPUT_DIR, "multiome_obj_markers_SCT.csv"), row.names = FALSE)
gc()

# =============================================================================
# PROJECT-SPECIFIC: cluster annotation
# =============================================================================
# Renames numeric cluster IDs to biologically meaningful cell type labels,
# based on marker genes / the reference signatures above. Replace the mapping
# below with the cell types relevant to your own biological system.
table(multiome_obj$seurat_clusters)    # inspect cluster sizes
Idents(multiome_obj) <- "seurat_clusters"
multiome_obj <- RenameIdents(multiome_obj, c(
  "0" = "GroupA",
  "1" = "GroupB",
  "2" = "GroupA",
  "3" = "GroupC",
  "4" = "GroupA",
  "5" = "GroupD",
  "6" = "unassigned",
  "7" = "unassigned"
))  # PROJECT-SPECIFIC: replace with your own cluster -> cell type mapping

multiome_obj[["Cell_Type"]] <- multiome_obj@active.ident
DimPlot(multiome_obj, reduction = "wnn.umap", label = TRUE, group.by = "Cell_Type", label.box = TRUE)

multiome_obj <- subset(multiome_obj, idents = c("GroupD", "unassigned"), invert = TRUE)  # drop low-quality/unwanted clusters (invert = TRUE means "exclude these, keep everything else")
DimPlot(multiome_obj, reduction = "wnn.umap", label = TRUE, group.by = "Cell_Type", label.box = TRUE)  # visualize after cleanup. The Cell_Type factor still contains the removed levels; use multiome_obj$Cell_Type <- droplevels(multiome_obj$Cell_Type) to drop them completely if needed.

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# GENE ACTIVITY (chromatin-based proxy for gene expression)
# =============================================================================
# gene_activity is a (genes x cells) count matrix (number of fragments)
gene_activity <- GeneActivity(
  object = multiome_obj,        # Seurat object
  assay = "ATAC",                # assay containing the peaks
  features = NULL,               # if NULL, uses all genes from the annotation
  extend.upstream = 2000,        # include 2 kb upstream of the TSS
  extend.downstream = 0          # gene body + upstream only
)
# NOTE: with ATAC-only data (no RNA), gene activity lets you assign an identity to
# cells (e.g. if TCF4 chromatin is "active", that cell is likely a certain cell type).
# In a multiomic setting it works as a cross-check: you can see whether marker
# chromatin accessibility matches RNA expression in the same clusters. Gene activity
# is therefore a bridge between the chromatin world and the gene-centric world,
# helping interpret ATAC clusters using the same genes used for RNA.

multiome_obj[["gene_activity"]] <- CreateAssayObject(counts = gene_activity)   # new assay with raw counts
multiome_obj <- NormalizeData(multiome_obj, assay = "gene_activity", normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(multiome_obj) <- "gene_activity"
# PROJECT-SPECIFIC: marker gene panels for each cell type, taken from a reference
# paper. Replace with the marker genes relevant to your own biological system.
reference_markers_groupA <- c("GENE1", "GENE2", "GENE3", "GENE4")
reference_markers_groupB <- c("GENE1", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9")
reference_markers_groupC <- c("GENE1", "GENE5", "GENE6", "GENE7", "GENE10", "GENE11")

# Feature plots on the integrated (WNN) UMAP - reflects the final populations
fp_groupA <- FeaturePlot(multiome_obj, features = reference_markers_groupA, reduction = "wnn.umap", pt.size = 0.1, max.cutoff = "q90", ncol = 2)
fp_groupB <- FeaturePlot(multiome_obj, features = reference_markers_groupB, reduction = "wnn.umap", pt.size = 0.1, max.cutoff = "q90", ncol = 3)
fp_groupC <- FeaturePlot(multiome_obj, features = reference_markers_groupC, reduction = "wnn.umap", pt.size = 0.1, max.cutoff = "q90", ncol = 3)
fp_groupA
fp_groupB
fp_groupC

# =============================================================================
# ADDITIONAL RNA MARKER EXPLORATION
# =============================================================================
# When a dataset is noisy, it can help to inspect additional RNA-level markers
# to get a clearer picture of the clusters.
DefaultAssay(object = multiome_obj) <- "RNA"
Idents(multiome_obj) <- "Cell_Type"
multiome_obj <- NormalizeData(multiome_obj, assay = "RNA")
multiome_obj_markers_RNA <- FindAllMarkers(multiome_obj, only.pos = TRUE, assay = "RNA", min.pct = 0.2)
multiome_obj_markers_RNA %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
multiome_obj_markers_RNA <- subset(multiome_obj_markers_RNA, subset = multiome_obj_markers_RNA$p_val_adj < 0.05)
# PROJECT-SPECIFIC: canonical marker genes for your system, e.g. from a reference
# publication describing the cell populations of interest.
canonical_markers <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10")

VlnPlot(multiome_obj, features = c("nFeature_RNA"), pt.size = 0)
DotPlot(multiome_obj, features = canonical_markers) + RotatedAxis()

FeaturePlot(multiome_obj, features = canonical_markers, reduction = "umap.rna", pt.size = 0.1, max.cutoff = "q90", ncol = 3)

write.csv(multiome_obj_markers_RNA, file.path(OUTPUT_DIR, "multiome_obj_markers_RNA.csv"), row.names = FALSE)
saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# PROJECT-SPECIFIC: SNP/variant focus — regulatory element overlap
# =============================================================================
# The following sections integrate genetic variants (e.g. GWAS hits) with
# chromatin accessibility and transcription-factor motifs. Replace the two
# input variant tables and downstream group labels with those relevant to
# your own study; consider whether the group names themselves reveal
# unpublished results before making this section public.
DefaultAssay(multiome_obj) <- "ATAC"
Idents(multiome_obj) <- multiome_obj$Cell_Type
table(Idents(multiome_obj))

df_groupA <- read_excel(SNP_TABLE_GROUP_A_FILE, sheet = "GroupA")   # read the variant tables (first sheet only, containing all variants)
df_groupB <- read_excel(SNP_TABLE_GROUP_B_FILE, sheet = "GroupB")

# Window size around each variant for coverage plots
window_bp <- 2000   # number of bases to show left and right of the variant

# Create output folders for the coverage plots
dir.create(COVERAGE_GROUP_A_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COVERAGE_GROUP_B_DIR, showWarnings = FALSE, recursive = TRUE)

# Function to draw and save the CoveragePlot for a single variant
plot_snp_coverage <- function(chrom, start, end, snp_id, out_dir) {    # chrom, start, end: genome coordinates of the variant, snp_id: rs identifier, out_dir: where to save the PNG
  region_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start - window_bp, end = end + window_bp))  # build the GRanges for the region to display

  p <- CoveragePlot(
    object   = multiome_obj,
    region   = region_gr,
    group.by = NULL,          # NULL uses the active identity (cell types)
    extend.upstream   = 0,    # don't extend beyond the region limits
    extend.downstream = 0,
    annotation = TRUE          # show gene names
  )

  snp_center <- floor((start + end) / 2)  # add a dashed red line at the center of the variant
  p <- p & geom_vline(xintercept = snp_center, linetype = "dashed", color = "red", linewidth = 0.6)
  p <- p + plot_annotation(title = paste0(snp_id, " (", chrom, ":", start, "-", end, ")"), theme = theme(plot.title = element_text(hjust = 0.5)))  # add a single title for the whole figure

  file_name <- file.path(out_dir, paste0(snp_id, ".png"))   # save as PNG named after the variant ID
  ggsave(file_name, plot = p, width = 10, height = 6)
  cat("  Saved", file_name, "\n")
}

# Loop over variant group A
cat("===== Processing variants from df_groupA =====\n")
for (i in seq_len(nrow(df_groupA))) {
  plot_snp_coverage(chrom = df_groupA$chrom[i],
                     start = df_groupA$start[i],
                     end = df_groupA$end[i],
                     snp_id = df_groupA$rsids[i],
                     out_dir = COVERAGE_GROUP_A_DIR)
}
# Loop over variant group B
cat("\n===== Processing variants from df_groupB =====\n")
for (i in seq_len(nrow(df_groupB))) {
  plot_snp_coverage(chrom = df_groupB$chrom[i],
                     start = df_groupB$start[i],
                     end = df_groupB$end[i],
                     snp_id = df_groupB$rsids[i],
                     out_dir = COVERAGE_GROUP_B_DIR)
}

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# TRANSCRIPTION FACTOR MOTIF ENRICHMENT WITH chromVAR
# =============================================================================
# Goal: compute, for every cell, an accessibility score for hundreds of
# transcription factor motifs, then integrate with gene expression (RNA) to
# identify master regulators of each population.
# NOTE: since we use the peaks called from this dataset, motifs are searched
# within those peaks only.
DefaultAssay(multiome_obj) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = SPECIES_TAXID, all_versions = FALSE))   # get motif PWMs from the JASPAR2020 database

# Build a binary (presence/absence) matrix indicating, for each peak (row) and
# each motif (column), whether the peak sequence contains that motif.
motif.matrix <- CreateMotifMatrix(features = granges(multiome_obj[["ATAC"]]), pwm = pwm_set, genome = GENOME_BUILD, use.counts = FALSE)   # binary, not counts
# motif.matrix is a temporary object before being stored in the Seurat object

motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)  # wrap the matrix into a Signac Motif object
multiome_obj <- SetAssayData(multiome_obj, assay = 'ATAC', layer = 'motifs', new.data = motif.object) # use `layer` instead of `slot` (depends on your Seurat version)

# Run chromVAR: computes an accessibility score for each motif in each cell
# (deviations from background).
multiome_obj <- RunChromVAR(object = multiome_obj, genome = BSgenome.Hsapiens.UCSC.hg38)  # after this step (30-60 min), a new "chromvar" assay is added containing a (motifs x cells) matrix of enrichment/depletion scores (normalized values).

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# Example: inspect a single motif/gene pair (replace "CEBPB" with a
# transcription factor relevant to your analysis)
example_tf <- "CEBPB"
motif.name <- ConvertMotifID(multiome_obj, name = example_tf)
gene_plot <- FeaturePlot(multiome_obj, features = paste0("sct_", example_tf), reduction = 'wnn.umap')
motif_plot <- FeaturePlot(multiome_obj, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

# =============================================================================
# CANDIDATE TF IDENTIFICATION PER CELL TYPE (RNA + motif evidence combined)
# =============================================================================
Idents(multiome_obj) <- multiome_obj$Cell_Type

# RNA-based markers
markers_rna <- presto::wilcoxauc(  # differential testing with presto (wilcoxauc); compares each cell type against all others, for both RNA and motif accessibility
  multiome_obj,
  group_by = 'Cell_Type',
  assay = 'data',           # layer to read data from (e.g. data, counts, scale.data)
  seurat_assay = 'SCT'      # assay to pull data from
)

# Motif-based markers (chromVAR assay)
markers_motifs <- presto::wilcoxauc(
  multiome_obj,
  group_by = 'Cell_Type',
  assay = 'data',
  seurat_assay = 'chromvar'
)

colnames(markers_rna)    <- paste0("RNA.", colnames(markers_rna))  # rename columns to distinguish the two result sets
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature  # add gene-name columns (for RNA) and TF names (converted from motif IDs, e.g. "MAXXXX.X" -> "TCF4")
markers_motifs$gene <- ConvertMotifID(multiome_obj, id = markers_motifs$motif.feature)  # if conversion fails (motif not linked to a gene), the original ID is kept

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# Function to aggregate results and find the top TFs for a given cell type.
# Input: group name (e.g. "GroupA"), padj threshold (default 0.01).
# Output: dataframe sorted by avg_auc (average of RNA.auc and motif.auc)
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  # Select positive (logFC > 0) and significant RNA markers
  ctmarkers_rna <- dplyr::filter(
    markers_rna,
    RNA.group == celltype,
    RNA.padj < padj.cutoff,
    RNA.logFC > 0
  ) %>% arrange(-RNA.auc)

  # Select positive and significant motif markers
  ctmarkers_motif <- dplyr::filter(
    markers_motifs,
    motif.group == celltype,
    motif.padj < padj.cutoff,
    motif.logFC > 0
  ) %>% arrange(-motif.auc)

  # Inner join: keep only TFs significant in both tests
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c("RNA.group", "gene", "RNA.auc", "RNA.pval")],
    y = ctmarkers_motif[, c("motif.group", "gene", "motif.auc", "motif.pval")],
    by = "gene"
  )

  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2  # average AUC (area under the curve) - a measure of the TF's discriminative power (both as RNA and as motif)
  top_tfs <- arrange(top_tfs, -avg_auc)  # sort by descending avg_auc

  return(top_tfs)
}

# Apply the function to each of your populations (PROJECT-SPECIFIC group names)
top_tfs_groupA <- topTFs("GroupA")
top_tfs_groupB <- topTFs("GroupB")
top_tfs_groupC <- topTFs("GroupC")

# Show the top 5 TFs for each cell type
cat("\n=== Top TFs for GroupA ===\n")
print(head(top_tfs_groupA, 5))
cat("\n=== Top TFs for GroupB ===\n")
print(head(top_tfs_groupB, 5))
cat("\n=== Top TFs for GroupC ===\n")
print(head(top_tfs_groupC, 5))

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# -----------------------------------------------------------------------
#               TO RELOAD THE OBJECT LATER, USE:
#               multiome_obj <- readRDS(RDS_FILE)
# -----------------------------------------------------------------------

# =============================================================================
# SNP / OPEN CHROMATIN OVERLAP AND TF IDENTIFICATION
# =============================================================================
# Prepare ATAC peaks and extract GRanges from the ATAC assay
peaks_gr <- granges(multiome_obj[["ATAC"]])  # all peaks remaining after filtering, converted to GRanges
peak_names <- rownames(multiome_obj[["ATAC"]])  # retrieve original peak names
if (is.null(peak_names)) {  # if peak names are missing, build them from genomic coordinates
  peak_names <- paste0(seqnames(peaks_gr), ":", start(peaks_gr), "-", end(peaks_gr))
}
names(peaks_gr) <- peak_names   # assign peak names to the GRanges object (so each peak has a text identifier)

# Build a GRanges for each of the two variant sets: chromosome, start, end, keep the rs ID as metadata
snp_gr_groupA <- GRanges(seqnames = df_groupA$chrom, ranges = IRanges(start = df_groupA$start, end = df_groupA$end), rsid = df_groupA$rsids)
snp_gr_groupB <- GRanges(seqnames = df_groupB$chrom, ranges = IRanges(start = df_groupB$start, end = df_groupB$end), rsid = df_groupB$rsids)

# Overlap variants and peaks. findOverlaps returns a pair of indices for each overlap.
overlaps_groupA <- findOverlaps(snp_gr_groupA, peaks_gr)
overlaps_groupB <- findOverlaps(snp_gr_groupB, peaks_gr)

# Build a dataframe listing, for each SNP, the overlapping peaks. Extracts the rs ID
# of the variant (from the original SNP) and the peak name (using the peak index).
snp_to_peak_groupA <- data.frame(snp = snp_gr_groupA$rsid[queryHits(overlaps_groupA)], peak = names(peaks_gr)[subjectHits(overlaps_groupA)], stringsAsFactors = FALSE)  # queryHits = variant index, subjectHits = peak index
snp_to_peak_groupB <- data.frame(snp = snp_gr_groupB$rsid[queryHits(overlaps_groupB)], peak = names(peaks_gr)[subjectHits(overlaps_groupB)], stringsAsFactors = FALSE)

# SNP-peak overlap (numeric indices). Extracts the numeric variant index and the peak index.
snp_peak_idx_groupA <- data.frame(snp = snp_gr_groupA$rsid[queryHits(overlaps_groupA)], peak_idx = subjectHits(overlaps_groupA), stringsAsFactors = FALSE)
snp_peak_idx_groupB <- data.frame(snp = snp_gr_groupB$rsid[queryHits(overlaps_groupB)], peak_idx = subjectHits(overlaps_groupB), stringsAsFactors = FALSE)

# Sanity check: number of variants falling within peaks
cat("Group A variants within peaks:", nrow(snp_peak_idx_groupA), "\n")
cat("Group B variants within peaks:", nrow(snp_peak_idx_groupB), "\n")

# Motif matrix (motifs x peaks) and motif -> TF mapping
DefaultAssay(multiome_obj) <- "ATAC"
motif_matrix <- GetAssayData(multiome_obj, assay = "ATAC", layer = "motifs") # motif_matrix is the same data retrieved from the Seurat object
dim(motif_matrix)     # peaks x motifs
motif_ids <- colnames(motif_matrix)   # column names are motif identifiers (e.g. MA0102.4)
motif_to_tf <- setNames(ConvertMotifID(multiome_obj, id = motif_ids), motif_ids)   # convert each motif ID to the corresponding transcription factor name (via Signac)
motif_matrix <- GetMotifData(multiome_obj, assay = "ATAC")
if (inherits(motif_matrix, "lgCMatrix")) {     # if logical (lgCMatrix), convert to numeric (dgCMatrix)
  motif_matrix <- as(motif_matrix, "dgCMatrix")
}
# Triplet: i = peak index (row), j = motif index (column)
triplet <- Matrix::summary(motif_matrix)   # extract triplets (row index, column index, value) of non-zero elements in the sparse matrix. summary() on a dgCMatrix returns a data.frame with columns i (peak index), j (motif index), x (value)

# Function to link each SNP (via its peak_idx) to the TFs present in that peak
merge_snp_tf <- function(snp_peak_idx, triplet, motif_to_tf, motif_ids) {
  merged <- merge(snp_peak_idx, triplet, by.x = "peak_idx", by.y = "i", all.x = TRUE)    # merge the SNP-peak-index dataframe with the triplet, using the peak index as key. A SNP can have multiple rows if its peak contains multiple motifs (one per row).
  merged$TF <- motif_to_tf[motif_ids[merged$j]]    # add the TF column: takes the motif ID (from the triplet's j column) and converts it to a TF name
  merged <- merged[, c("snp", "peak_idx", "TF")]     # keep only the columns of interest: SNP, peak index, TF (NA if no motif)
  return(merged)
}

# Apply the function
res_groupA <- merge_snp_tf(snp_peak_idx_groupA, triplet, motif_to_tf, motif_ids)
res_groupB <- merge_snp_tf(snp_peak_idx_groupB, triplet, motif_to_tf, motif_ids)

res_groupA$peak_name <- peak_names[res_groupA$peak_idx]  # add peak coordinates
res_groupB$peak_name <- peak_names[res_groupB$peak_idx]
head(res_groupA)
head(res_groupB)

# =============================================================================
# FINAL TABLE: SNP - TF - ACCESSIBILITY
# =============================================================================
# Mean per-peak accessibility across cell types
DefaultAssay(multiome_obj) <- "ATAC"
atac_mat <- GetAssayData(multiome_obj, layer = "data") # retrieve the normalized (TF-IDF) accessibility matrix: a sparse matrix, rows = peaks, columns = cells. Each cell contains a number indicating how "open" that peak is in that cell (higher = more accessible).
cell_types <- c("GroupA", "GroupB", "GroupC")   # PROJECT-SPECIFIC: your final cell type labels
acc_mean <- list()     # results stored in a list, e.g. acc_mean$GroupA is a vector with one value per peak: the mean accessibility of that peak in GroupA cells
for (ct in cell_types) {
  cells <- WhichCells(multiome_obj, expression = Cell_Type == ct)  # find the cells belonging to each cell type
  acc_mean[[ct]] <- Matrix::rowMeans(atac_mat[, cells, drop = FALSE])  # extract the sub-matrix restricted to those cells (columns) and compute, for each peak (row), the mean accessibility across cells of that type
}

# NOTE: this is an intermediate step that only keeps the first peak per SNP; the
# actual final aggregation is done in the clean_snp_set() function below.
# Otherwise you'd need: snp_unique_X <- res_X[!duplicated(res_X[, c("snp","peak_idx")]), c("snp","peak_idx","peak_name")]
snp_unique_groupA <- res_groupA[!duplicated(res_groupA$snp), c("snp", "peak_idx", "peak_name")]  # res_groupA/res_groupB (built above) contain the SNP -> peak mapping and, for each peak, the TFs binding it. Each row is a (snp, TF) pair, so the same SNP can appear on multiple rows (one per TF in its peak).
snp_unique_groupB <- res_groupB[!duplicated(res_groupB$snp), c("snp", "peak_idx", "peak_name")] # !duplicated(res_X$snp) keeps only the first occurrence of each SNP (removes duplicates).
# snp_unique_groupA has one row per SNP, with columns: snp (variant ID), peak_idx (numeric peak index, to access the matrix), peak_name (e.g. chr3-45890383-45891097)

# Add accessibility values
snp_unique_groupA$access_GroupA <- acc_mean$GroupA[snp_unique_groupA$peak_idx]
snp_unique_groupA$access_GroupB <- acc_mean$GroupB[snp_unique_groupA$peak_idx]
snp_unique_groupA$access_GroupC <- acc_mean$GroupC[snp_unique_groupA$peak_idx]

snp_unique_groupB$access_GroupA <- acc_mean$GroupA[snp_unique_groupB$peak_idx]
snp_unique_groupB$access_GroupB <- acc_mean$GroupB[snp_unique_groupB$peak_idx]
snp_unique_groupB$access_GroupC <- acc_mean$GroupC[snp_unique_groupB$peak_idx]

# Aggregate TFs per SNP using the original dataframes (res_groupA/res_groupB, which
# have multiple rows per SNP). For each SNP, if its peak contains multiple motifs,
# they are aggregated into a comma-separated string.
tf_agg_groupA <- aggregate(TF ~ snp, data = res_groupA, FUN = function(x) paste(unique(na.omit(x)), collapse = ", "))
tf_agg_groupB <- aggregate(TF ~ snp, data = res_groupB, FUN = function(x) paste(unique(na.omit(x)), collapse = ", "))

# Merge the TF aggregation with the unique SNP tables
snp_final_groupA <- merge(snp_unique_groupA, tf_agg_groupA, by = "snp", all.x = TRUE)
snp_final_groupB <- merge(snp_unique_groupB, tf_agg_groupB, by = "snp", all.x = TRUE)

write.csv(snp_final_groupA, file.path(OUTPUT_DIR, "SNP_groupA_summary.csv"), row.names = FALSE)   # intermediate results
write.csv(snp_final_groupB, file.path(OUTPUT_DIR, "SNP_groupB_summary.csv"), row.names = FALSE)   # intermediate results

saveRDS(multiome_obj, file = RDS_FILE)
gc()

# =============================================================================
# FINAL CLEANUP INTO A SINGLE SUMMARY TABLE
# =============================================================================
# Takes res_df (columns: snp, peak_idx, TF, peak_name) and acc_mean_list (with
# per-cell-type accessibility) and returns a clean, one-row-per-SNP dataframe.
# NOTE: this final table is NOT built from snp_final_groupA/snp_final_groupB;
# instead it re-aggregates all peaks for each SNP directly (group_by(snp) +
# summarise with max()/mean() across all peak_idx values).
clean_snp_set <- function(res_df, acc_mean_list) {

  # Add per-peak accessibility (mean value per cell type)
  res_df$acc_GroupA <- acc_mean_list$GroupA[res_df$peak_idx]
  res_df$acc_GroupB <- acc_mean_list$GroupB[res_df$peak_idx]
  res_df$acc_GroupC <- acc_mean_list$GroupC[res_df$peak_idx]

  # Group by SNP and compute aggregate statistics
  cleaned <- res_df %>% group_by(snp) %>% summarise(
      # Unique TFs (drop NA and duplicates)
      TF = paste(unique(na.omit(TF)), collapse = ", "),
      # All peak names involved
      peak_names = paste(unique(peak_name), collapse = "; "),
      # MAXIMUM accessibility across the SNP's peaks
      acc_GroupA_max = max(acc_GroupA, na.rm = TRUE),
      acc_GroupB_max = max(acc_GroupB, na.rm = TRUE),
      acc_GroupC_max = max(acc_GroupC, na.rm = TRUE),
      # MEAN accessibility across the SNP's peaks
      acc_GroupA_mean = mean(acc_GroupA, na.rm = TRUE),
      acc_GroupB_mean = mean(acc_GroupB, na.rm = TRUE),
      acc_GroupC_mean = mean(acc_GroupC, na.rm = TRUE),

      .groups = "drop"
    )

  return(cleaned)
}

# Apply the function to both variant sets (assumes res_groupA, res_groupB,
# acc_mean already exist)
cleaned_groupA <- clean_snp_set(res_groupA, acc_mean)   # considers all peaks, since res_groupA contains one row per SNP-peak-TF combination
cleaned_groupB <- clean_snp_set(res_groupB, acc_mean)

# Add a column to remember which set each SNP came from
cleaned_groupA$origin <- "GroupA"
cleaned_groupB$origin <- "GroupB"

# Merge the two tables
combined <- full_join(cleaned_groupA, cleaned_groupB, by = "snp", suffix = c(".GroupA", ".GroupB"))   # full join keeps all rows from both dataframes
final <- combined

final$origin_combined <- NA   # combined origin column: whether the SNP is found in GroupA, GroupB, or both
final$origin_combined[!is.na(final$origin.GroupA) & !is.na(final$origin.GroupB)] <- "both"
final$origin_combined[!is.na(final$origin.GroupA) & is.na(final$origin.GroupB)] <- "GroupA"
final$origin_combined[is.na(final$origin.GroupA) & !is.na(final$origin.GroupB)] <- "GroupB"

final$TF <- NA  # TF: union of GroupA and GroupB TFs (if both present)
for (i in 1:nrow(final)) {
  tfs <- c()
  if (!is.na(final$TF.GroupA[i])) tfs <- c(tfs, strsplit(final$TF.GroupA[i], ", ")[[1]])
  if (!is.na(final$TF.GroupB[i])) tfs <- c(tfs, strsplit(final$TF.GroupB[i], ", ")[[1]])
  tfs <- unique(tfs[!is.na(tfs)])
  final$TF[i] <- paste(tfs, collapse = ", ")
}

final$peak_names <- apply(final[, c("peak_names.GroupA", "peak_names.GroupB")], 1, function(x) {  # peak coordinates
  peaks <- unlist(strsplit(na.omit(x), "; "))
  peaks <- unique(peaks)
  if (length(peaks) == 0) return(NA)
  paste(peaks, collapse = "; ")
})

# Accessibility: take the max for *_max columns, the mean for *_mean columns.
# Maximum accessibility indicates the most open peak among those overlapping the
# SNP (i.e. the strongest regulatory signal), while mean accessibility reflects
# how open, on average, all peaks containing that SNP are.
for (ct in cell_types) {
  final[[paste0("acc_", ct, "_max")]] <- pmax(final[[paste0("acc_", ct, "_max.GroupA")]], final[[paste0("acc_", ct, "_max.GroupB")]], na.rm = TRUE)
  final[[paste0("acc_", ct, "_mean")]] <- rowMeans(final[, c(paste0("acc_", ct, "_mean.GroupA"), paste0("acc_", ct, "_mean.GroupB"))], na.rm = TRUE)
}

# Build the final dataframe by selecting the relevant columns.
# final_table is derived from cleaned_groupA/cleaned_groupB, not from the
# intermediate snp_final_groupA/snp_final_groupB tables.
final_table_cols <- c("snp", "origin_combined", "TF", "peak_names",
                       unlist(lapply(cell_types, function(ct) c(paste0("acc_", ct, "_max"), paste0("acc_", ct, "_mean")))))
final_table <- final[, final_table_cols]
colnames(final_table)[colnames(final_table) == "origin_combined"] <- "origin"   # rename origin_combined -> origin
final_table <- final_table[order(final_table$snp), ]  # sort by SNP
head(final_table)

write.csv(final_table, file.path(OUTPUT_DIR, "SNP_combined.csv"), row.names = FALSE)
write_xlsx(final_table, file.path(OUTPUT_DIR, "SNP_combined.xlsx"))

# NOTE: for each SNP, the final table lists all peaks containing it, all TFs
# binding those peaks, and all accessibility values across cell types, condensed
# into a single row. Biologically, the final table represents, for each SNP:
# where it falls, which TFs might bind it, how open that region is across cell
# types, whether the signal is strong (max), and whether the region is
# generally accessible (mean). This is a common strategy for prioritizing
# regulatory SNPs in GWAS + ATAC + motif-enrichment analyses.

# Merge the final table with a regulome annotation table to combine all
# available information
regulome_combined <- read.table(REGULOME_COMBINED_FILE, sep = "\t", header = TRUE, fill = TRUE)
colnames(regulome_combined)
colnames(final_table)
final_combined <- merge(regulome_combined, final_table, by.x = "rsids", by.y = "snp", all = TRUE)   # all = TRUE -> full outer join, keeping all rows from both tables

write_xlsx(final_combined, file.path(OUTPUT_DIR, "final_combined.xlsx"))

# =============================================================================
# TF EXPRESSION DOTPLOTS
# =============================================================================
# Clean up the TF names found in the final table
tf_col <- final_combined$TF[!is.na(final_combined$TF)]   # take the TF column, drop empty entries
all_tf <- unlist(strsplit(paste(tf_col, collapse = ", "), ", "))     # collapse everything into one string, then split into individual TFs
all_tf <- trimws(gsub("::.*", "", gsub("\\(.*", "", all_tf)))  # clean up names (remove "(var.2)", "::RXRA", extra spaces)
unique_tfs <- sort(unique(all_tf))   # unique, sorted list

writeLines(unique_tfs, file.path(OUTPUT_DIR, "unique_TFs_list.txt"))   # save the list

# Check which TFs are actually present in the dataset
DefaultAssay(multiome_obj) <- "SCT"
Idents(multiome_obj) <- multiome_obj$Cell_Type
tfs_present <- intersect(unique_tfs, rownames(multiome_obj))   # TFs actually present in the dataset
# tfs_present <- head(tfs_present, 50)   # uncomment to limit the number (e.g. top 50)

writeLines(tfs_present, file.path(OUTPUT_DIR, "TF_present_list.txt"))   # save the list

p_seurat <- DotPlot(multiome_obj, features = tfs_present, group.by = "Cell_Type") + RotatedAxis() + ggtitle("TF expression (Seurat)") + theme(axis.text.x = element_text(size = 7))
if (requireNamespace("SCpubr", quietly = TRUE)) {
  p_scpubr <- SCpubr::do_DotPlot(sample = multiome_obj, features = tfs_present, group.by = "Cell_Type", assay = "SCT", plot.title = "TF expression (SCpubr)", axis.text.x.angle = 45)
} else {
  message("SCpubr is not installed. Install it with: install.packages('SCpubr')")
}
p_seurat
p_scpubr

# The TF list can be long, so we intersect it with differentially expressed
# marker genes per cell type: this tells us which candidate regulatory TFs are
# also differentially expressed marker genes for that cell type.
DefaultAssay(multiome_obj) <- "SCT"    # extract significant marker genes per cluster (already computed in multiome_obj_markers_SCT)
tf_intersections <- lapply(cell_types, function(ct) {
  intersect(multiome_obj_markers_SCT$gene[multiome_obj_markers_SCT$cluster == ct], tfs_present)
})
names(tf_intersections) <- cell_types
print(tf_intersections)

# Build and save a summary dataframe
tf_intersection_df <- data.frame(
  Cell_Type = rep(cell_types, times = sapply(tf_intersections, length)),
  TF = unlist(tf_intersections)
)
write_xlsx(tf_intersection_df, file.path(OUTPUT_DIR, "TF_marker_intersection.xlsx"))

# Dotplot restricted to the shared/common TFs
common_tfs <- unique(tf_intersection_df$TF)
cat("Number of common TFs:", length(common_tfs), "\n")
p_seurat_common <- DotPlot(multiome_obj, features = common_tfs, group.by = "Cell_Type") + RotatedAxis() + ggtitle("TF expression") + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), axis.text.y = element_text(size = 8))
p_scpubr_common <- SCpubr::do_DotPlot(sample = multiome_obj, features = common_tfs, group.by = "Cell_Type", assay = "SCT", plot.title = "TF expression", axis.text.x.angle = 45, font.size = 9)
p_seurat_common
p_scpubr_common

# =============================================================================
# Integration with bulk ATAC-seq data
# =============================================================================
# (Placeholder — integrate with the bulk ATAC-seq pipeline results here, e.g.
# comparing single-cell peak accessibility against bulk ATAC-seq bigWig/peak
# calls from the companion pipeline.)
