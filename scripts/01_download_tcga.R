# ============================================
# scripts/01_download_tcga.R
# Single script to download TCGA-BRCA 450K methylation,
# RNA-seq counts, and clinical + PAM50 subtype.
# ============================================

# 0. Setup -------------------------------------------------------------------

project <- "TCGA-BRCA"
out_dir <- "data"
platform <- "450K"  # only 450K

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cleanup_gdc <- function() {
  tmp_dirs <- list.dirs(".", full.names = TRUE, recursive = FALSE)
  tmp_dirs <- tmp_dirs[grepl("^[0-9a-f\\-]{36}$", basename(tmp_dirs))]
  if (length(tmp_dirs)) unlink(tmp_dirs, recursive = TRUE)
  tmp_tars <- list.files(".", pattern = "\\.tar\\.gz$", full.names = TRUE)
  if (length(tmp_tars)) unlink(tmp_tars)
}

# 1. Install / load required packages ---------------------------------------

required_pkgs <- c("BiocManager", "Biobase", "TCGAbiolinks", "SummarizedExperiment",
                   "dplyr","readr","arrow","Matrix")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") install.packages("BiocManager")
    else if (pkg %in% c("arrow","Matrix")) install.packages(pkg)
    else BiocManager::install(pkg, ask = FALSE)
  }
}
suppressMessages({
  library(Biobase)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(arrow)
  library(Matrix)
})

# 2. Download DNA methylation (450K) ----------------------------------------

meth_csv     <- file.path(out_dir, "methylation_450k.csv")
meth_parquet <- file.path(out_dir, "methylation_450k.parquet")

if (!file.exists(meth_parquet)) {
  message("Downloading 450K methylation β-values…")
  query_met <- GDCquery(
    project = project,
    data.category = "DNA Methylation",
    data.type     = "Methylation Beta Value",
    platform      = "Illumina Human Methylation 450"
  )

  for (cs in c(5,2,1)) {
    cleanup_gdc()
    try({ GDCdownload(query_met, method="api", files.per.chunk=cs); break }, silent=TRUE)
  }
  message("Preparing methylation data…")
  meth_se <- GDCprepare(query_met)
  mat_den <- assay(meth_se)
  write.csv(mat_den, meth_csv, row.names=TRUE)         # optional backup
  df_meth <- data.frame(ProbeID=rownames(mat_den), mat_den, check.names=FALSE)
  write_parquet(df_meth, meth_parquet)
  cleanup_gdc()
} else {
  message("Methylation already downloaded: ", meth_parquet)
}

# 3. Download RNA-seq counts ------------------------------------------------

rna_csv     <- file.path(out_dir, "rna_counts.csv")
rna_parquet <- file.path(out_dir, "rna_counts.parquet")

if (!file.exists(rna_parquet)) {
  message("Downloading RNA-seq STAR counts…")
  query_rna <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query_rna)
  rna_se <- GDCprepare(query_rna)
  mat_rna <- assay(rna_se)
  write.csv(mat_rna, rna_csv, row.names=TRUE)           # optional backup
  df_rna <- data.frame(Gene=rownames(mat_rna), mat_rna, check.names=FALSE)
  write_parquet(df_rna, rna_parquet)
  cleanup_gdc()
} else {
  message("RNA counts already downloaded: ", rna_parquet)
}


message("01_download_tcga.R completed.")
