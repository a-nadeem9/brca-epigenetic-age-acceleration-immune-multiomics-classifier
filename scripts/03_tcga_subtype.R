# Load required libraries
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)
library(dplyr)
library(readr)

# Step 1: Load matched patient IDs
matched <- read_csv("data/matched_patients.csv", show_col_types = FALSE)
matched$patient_id <- substr(matched$patient_id, 1, 12)

# Step 2: Get clinical data from GDC
clinical_data <- tryCatch({
  GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
}, error = function(e) {
  message("Error retrieving clinical data: ", e$message)
  stop("Exiting due to clinical data retrieval error.")
})

# Step 3: Get BRCA subtype annotations (PAM50)
subtypes <- TCGAquery_subtype("BRCA") 

# Step 4: Merge subtype into clinical data
clinical_data$submitter_id <- substr(clinical_data$submitter_id, 1, 12)
merged <- clinical_data %>%
  inner_join(subtypes[, c("patient", "BRCA_Subtype_PAM50")],
             by = c("submitter_id" = "patient"))

# Step 5: Filter for matched patients only
final_data <- merged %>%
  filter(submitter_id %in% matched$patient_id)

# Optional: rename for consistency
colnames(final_data)[colnames(final_data) == "BRCA_Subtype_PAM50"] <- "Subtype"

# Step 6: Flatten to remove list-columns before writing
flat_data <- as.data.frame(final_data)
flat_data <- flat_data[, sapply(flat_data, function(col) !is.list(col))]

# Step 7: Save final data
write.csv(flat_data, "data/tcga_brca_clinical_matched_with_subtypes.csv", row.names = FALSE)

cat("Saved clinical + subtype data for matched patients to: data/tcga_brca_clinical_matched_with_subtypes.csv\n")
