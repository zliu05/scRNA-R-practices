# 03_sample_level_analysis.R
# Sample-level embedding, diffusion map, and principal-curve pseudo-trajectory.

source("00_setup.R")

object <- readRDS(file.path(config$out_dir, "02_prepared.rds"))

# Landmark-by-sample matrix (cells near each landmark, per donor), Chi-squared normalized
sample_obj <- GenerateSampleObject(
  object = object, sketch.assay = "LANDMARK", nn.name = "weighted.nn",
  group.by = config$donor_col, rename.group.by = config$celltype_col,
  k.nn = 5, normalization.method = "ChiSquared",
  add.meta.data = TRUE, return.seurat = TRUE
)
sample_obj <- ScaleData(
  sample_obj, features = rownames(sample_obj), do.scale = FALSE, do.center = TRUE
)

# Diffusion map across donors: 10 components
sample_obj <- RunDiffusionMap(
  sample_obj, assay = "LMC", layer = "scale.data",
  features = rownames(sample_obj), ncomp = 10
)

# Add DC1-DC3 to the sample metadata
dc <- Embeddings(sample_obj, "DiffMap")[, 1:3]
sample_obj$DC1 <- dc[, 1]
sample_obj$DC2 <- dc[, 2]
sample_obj$DC3 <- dc[, 3]

# Plots: DC1 vs DC2 and DC1 vs DC3 colored by disease group
p1 <- DimPlot(sample_obj, reduction = "DiffMap", dims = c(1, 2), group.by = config$disease_col)
p2 <- DimPlot(sample_obj, reduction = "DiffMap", dims = c(1, 3), group.by = config$disease_col)
ggsave(file.path(config$out_dir, "diffmap_dc1_dc2.png"), p1, width = 6, height = 5)
ggsave(file.path(config$out_dir, "diffmap_dc1_dc3.png"), p2, width = 6, height = 5)

# DCs against time since onset (if the column exists)
if (config$time_col %in% colnames(sample_obj[[]])) {
  df <- sample_obj[[]]
  df[["DC3"]] <- sample_obj$DC3
  p3 <- ggplot(df, aes(x = .data[[config$time_col]], y = DC3,
                       color = .data[[config$disease_col]])) +
    geom_point() + theme_classic()
  ggsave(file.path(config$out_dir, "dc3_vs_time.png"), p3, width = 6, height = 4)
}

# Principal curve through DC1-DC2: a joint pseudo-trajectory across donors
princurve_res <- princurve::principal_curve(
  Embeddings(sample_obj, "DiffMap")[, 1:2],
  plot = FALSE, maxit = 20, trace = TRUE
)
sample_obj$traj_overall <- as.numeric(princurve_res$lambda)

saveRDS(sample_obj, file.path(config$out_dir, "03_sample_obj.rds"))
saveRDS(object, file.path(config$out_dir, "03_cell_obj.rds"))
