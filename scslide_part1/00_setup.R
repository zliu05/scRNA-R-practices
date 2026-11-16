# 00_setup.R
# Install (if missing) and load the packages used in scSLIDE vignette Part 1.
# Vignette: https://satijalab.org/seurat/articles/scslide_vignette

cran_pkgs <- c("Seurat", "princurve", "ggplot2", "ggridges", "remotes", "BiocManager")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("destiny", quietly = TRUE)) BiocManager::install("destiny")
if (!requireNamespace("BPCells", quietly = TRUE)) remotes::install_github("bnprks/BPCells/r")

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(princurve)
  library(ggplot2)
  library(ggridges)
})

# Shared settings used by all scripts
config <- list(
  data_file   = "Raw_count_Subset_104k_Ahern_2022_Cell.rds",  # download from the vignette's Dropbox link
  out_dir     = "output",
  donor_col   = "Donor",
  disease_col = "Disease",       # levels: HV, MILD, SEV, CRIT
  time_col    = "TimeSinceOnset",
  celltype_col = "azimuth_label"
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
