# 01_load_and_annotate.R
# Load the COVID-19 subset, normalize, annotate with Azimuth, pick HVGs.

source("00_setup.R")

object <- readRDS(config$data_file)
stopifnot(all(c(config$donor_col, config$disease_col) %in% colnames(object[[]])))

object <- NormalizeData(object)

# Pan-human Azimuth annotation through the cloud API. Needs network access and
# the AzimuthAPI package (remotes::install_github("satijalab/AzimuthAPI")).
# This creates the cell-type labels (azimuth_label) and the "azimuth_embed" reduction.
if (!requireNamespace("AzimuthAPI", quietly = TRUE)) {
  remotes::install_github("satijalab/AzimuthAPI")
}
object <- CloudAzimuth(object)
stopifnot(config$celltype_col %in% colnames(object[[]]),
          "azimuth_embed" %in% Reductions(object))

object <- FindVariableFeatures(object, nfeatures = 2000)

# Quick look at the design
print(table(object[[config$disease_col]][, 1], useNA = "ifany"))
print(length(unique(object[[config$donor_col]][, 1])))

saveRDS(object, file.path(config$out_dir, "01_annotated.rds"))
