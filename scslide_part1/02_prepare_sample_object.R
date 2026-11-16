# 02_prepare_sample_object.R
# Build the landmark assay and the cell-to-landmark WNN graph.
# Set use_wrapper <- TRUE for the one-call PrepareSampleObject() version,
# or FALSE to run the same steps manually.

source("00_setup.R")
use_wrapper <- FALSE

object <- readRDS(file.path(config$out_dir, "01_annotated.rds"))

if (use_wrapper) {
  object <- PrepareSampleObject(
    object = object, assay = "RNA", add.hvg = TRUE,
    group.by.CorTest = config$celltype_col, Y = config$disease_col,
    sketch.training = TRUE, group.by.Sketch = config$donor_col,
    ncells.per.group = 1000, ncells.landmark = 2000, ncomp = 20,
    pls.function = "cppls", k.nn = 5,
    name.reduction.1 = "azimuth_embed", dims.reduction.1 = 1:128,
    fix.wnn.weights = c(0.5, 0.5), verbose = TRUE
  )
} else {
  # 1. Uniform sketch: 1,000 cells per donor, balances donors for PLS training
  object <- SketchDataByGroup(
    object, assay = "RNA", ncells = 1000, sketched.assay = "TRAINING",
    method = "Uniform", group.by = config$donor_col, verbose = FALSE
  )

  # 2. Leverage-score sketch: 2,000 diverse landmark cells
  object <- SketchData(
    object, assay = "TRAINING", ncells = 2000, sketched.assay = "LANDMARK",
    method = "LeverageScore", verbose = FALSE
  )

  # 3. Disease-relevant genes (correlation with disease), merged with the HVGs
  deg_list <- QuickCorTest(
    object, assay = "TRAINING", layer = "data",
    group.by = config$celltype_col, Y = config$disease_col,
    cor.cutoff = 0.2, verbose = FALSE
  )
  all_HVG <- union(VariableFeatures(object), deg_list)
  message(length(all_HVG), " genes used for PLS")

  # 4. Supervised embedding: canonical PLS on the training cells
  DefaultAssay(object) <- "TRAINING"
  object <- ScaleData(object, features = all_HVG)
  object <- RunPLS(
    object = object, assay = "TRAINING", features = all_HVG,
    ncomp = 20, Y = c(config$disease_col), pls.function = "cppls", verbose = FALSE
  )

  # 5. Project the PLS model onto every cell
  proj.pls <- ProjectCellEmbeddings(
    query = object, query.assay = "RNA",
    reference = object, reference.assay = "TRAINING",
    reduction = "pls", dims = 1:20, scale = TRUE,
    normalization.method = "LogNormalize", verbose = TRUE
  )
  object[["proj.pls"]] <- CreateDimReducObject(
    embeddings = proj.pls, key = "pPLS_", assay = "RNA"
  )
  rm(proj.pls); invisible(gc())

  # 6. WNN graph from landmarks to all cells; Azimuth and PLS weighted 0.5 each
  DefaultAssay(object) <- "RNA"
  object <- FindmmNN(
    object, sketch.assay = "LANDMARK",
    reduction.list = list("azimuth_embed", "proj.pls"), k.nn = 5,
    nn.name = "weighted.nn", dims.list = list(1:128, 1:20),
    fix.wnn.weights = c(0.5, 0.5)
  )

  # Optional single-modality graph using PLS only
  object <- FindsmNN(
    object, sketch.assay = "LANDMARK", reduction = "proj.pls",
    k.nn = 5, dims = 1:20, nn.name = "single.nn"
  )
}

saveRDS(object, file.path(config$out_dir, "02_prepared.rds"))
