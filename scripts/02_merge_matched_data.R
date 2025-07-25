# ===============================
# scripts/tcga_merge_matched_data.R
# Fast, Verbose Merge of TCGA-BRCA RNA-seq & Methylation by Patient ID
# ===============================

# 1) Ensure data folder exists
if (!dir.exists("data")) {
  dir.create("data")
  cat("Created data/ folder\n")
} else {
  cat("data/ folder already exists\n")
}

# 2) Check inputs
if (!file.exists("data/rna_counts.csv") || !file.exists("data/methylation_450k.csv")) {
  stop("Missing RNA or methylation data. Run download scripts first.")
}

# 3) Skip if already merged
merged_files <- c("data/rna_matched.csv",
                  "data/meth_matched.csv",
                  "data/matched_patients.csv")
if (all(file.exists(merged_files))) {
  cat("Matched outputs already exist; skipping merge.\n")
  quit(status = 0)
}

library(data.table)

# 4) Read headers only to get sample barcodes
cat("4/9: Reading CSV headers to extract sample IDs…\n")
rna_hdr  <- names(fread("data/rna_counts.csv",   nrows = 0))
meth_hdr <- names(fread("data/methylation_450k.csv", nrows = 0))

rna_ids  <- substr(rna_hdr[-1],  1, 12)  # drop first column name
meth_ids <- substr(meth_hdr[-1], 1, 12)
cat("RNA headers found:", length(rna_ids), "samples\n")
cat("Meth headers found:", length(meth_ids), "samples\n")

# 5) Find common TFCA barcodes
cat("5/9: Finding common patient IDs…\n")
common_ids <- intersect(rna_ids, meth_ids)
cat("Matched patients:", length(common_ids), "\n")
if (length(common_ids) == 0) stop("No common patient IDs!")

# 6) Determine which columns to load
cat("6/9: Selecting only matched columns…\n")
# Include the first column (rownames) + those columns whose IDs match
rna_cols  <- c(rna_hdr[1],  rna_hdr[-1][rna_ids  %in% common_ids])
meth_cols <- c(meth_hdr[1], meth_hdr[-1][meth_ids %in% common_ids])

# 7) Load matched RNA-seq
cat("7/9: Loading RNA-seq subset (", length(rna_cols)-1, "samples)…\n")
dt_rna <- fread("data/rna_counts.csv", select = rna_cols)
# Build matrix
genes      <- dt_rna[[1]]
mat_rna    <- as.matrix(dt_rna[ , -1, with = FALSE])
rownames(mat_rna) <- genes

# 8) Load matched methylation
cat("8/9: Loading methylation subset (", length(meth_cols)-1, "samples)…\n")
dt_meth <- fread("data/methylation_450k.csv", select = meth_cols)
probes     <- dt_meth[[1]]
mat_meth   <- as.matrix(dt_meth[ , -1, with = FALSE])
rownames(mat_meth) <- probes

# 9) Write out matched objects
cat("9/9: Writing matched data to CSV\n")
# RNA
rna_df <- data.frame(Gene = rownames(mat_rna), mat_rna, check.names = FALSE)
fwrite(rna_df, "data/rna_matched.csv")
# Methylation
meth_df <- data.frame(ProbeID = rownames(mat_meth), mat_meth, check.names = FALSE)
fwrite(meth_df, "data/meth_matched.csv")
# Patient list
fwrite(data.table(patient_id = common_ids), "data/matched_patients.csv")

cat("Merge complete! Files saved in data/:\n")
cat(" -rna_matched.csv\n")
cat(" -meth_matched.csv\n")
cat(" -matched_patients.csv\n")
