if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

# Install all the heavy Bioc/CRAN deps plus the pipeline itself
BiocManager::install(c(
  "ggplot2","glmmTMB","ggpubr","devtools","magick","reshape2",
  "minfi","sesame","methylclock"
))
remotes::install_github("CastellaniLab/EpigeneticAgePipeline")

#  Load the pipeline
library(EpigeneticAgePipeline)


# Run main() on β-matrix
main(
  directory      = "data",      # folder containing betaValues.csv
  normalize      = FALSE,       # already normalized β-values
  useBeta        = TRUE,        # read from betaValues.csv
  arrayType      = "450K",      # HM450 array
  useSampleSheet = FALSE,       # no Sample_Sheet.csv needed
  doParallel     = FALSE,
  writeBeta      = FALSE,
  useAdult       = FALSE
)

ages <- read.delim("epigeneticAge.txt", stringsAsFactors = FALSE, check.names = FALSE)
print(names(ages))


# Install/load dependencies
if (!requireNamespace("dplyr", quietly=TRUE)) BiocManager::install("dplyr", ask=FALSE)
library(dplyr)

# Read the pipeline’s output
ages <- read.table(
  "epigeneticAge.txt",
  header         = TRUE,
  sep            = "",
  stringsAsFactors = FALSE,
  check.names    = FALSE
)

# Rename & select the three clocks
three_clock_df <- ages %>%
  rename(
    SampleID = id,
    PhenoAge = Levine
  ) %>%
  select(SampleID, Horvath, Hannum, PhenoAge)

# Write out CSV
if (!dir.exists("results")) dir.create("results")
write.csv(
  three_clock_df,
  file      = "results/three_clocks_pheno.csv",
  row.names = FALSE
)

message("Wrote results/three_clocks_pheno.csv with ", nrow(three_clock_df), " samples")
