# 04_differential_expression.R
# Pseudobulk DE along each diffusion component, CD14 monocytes.

source("00_setup.R")

object     <- readRDS(file.path(config$out_dir, "03_cell_obj.rds"))
sample_obj <- readRDS(file.path(config$out_dir, "03_sample_obj.rds"))

# Donor x cell-type pseudobulk
bulk_obj <- AggregateExpression(
  object, assays = "RNA",
  group.by = c(config$donor_col, config$celltype_col), return.seurat = TRUE
)

sub_bulk_obj <- subset(bulk_obj, subset = azimuth_label == "CD14 monocyte")

# Merge sample-level covariates and trajectory values into the pseudobulk metadata.
# Match on donor id; AggregateExpression can rewrite "_" as "-", so normalize both sides.
norm_id <- function(x) gsub("_", "-", as.character(x))
bulk_meta <- sub_bulk_obj[[]]
sample_meta <- sample_obj[[]]
idx <- match(norm_id(bulk_meta[[config$donor_col]]), norm_id(rownames(sample_meta)))
stopifnot(!anyNA(idx))

for (v in c("Age", "Sex", config$disease_col)) {
  if (v %in% colnames(sample_meta)) sub_bulk_obj[[v]] <- sample_meta[[v]][idx]
}
sub_bulk_obj$traj_DC1 <- sample_meta$DC1[idx]
sub_bulk_obj$traj_DC2 <- sample_meta$DC2[idx]
sub_bulk_obj$traj_DC3 <- sample_meta$DC3[idx]
sub_bulk_obj$log_count <- log(Matrix::colSums(
  LayerData(sub_bulk_obj, assay = "RNA", layer = "counts")
))

# DC2 and DC3 are fit on patients only (healthy donors removed)
non_health_idx <- which(sub_bulk_obj[[config$disease_col]][, 1] != "HV")

lv <- c("Age", "Sex", "log_count")
de_dc1 <- TrajDETest(sub_bulk_obj, assay = "RNA", layer = "counts",
                     traj.var = "traj_DC1", latent.vars = lv)
de_dc2 <- TrajDETest(sub_bulk_obj, assay = "RNA", layer = "counts",
                     traj.var = "traj_DC2", latent.vars = lv, samples = non_health_idx)
de_dc3 <- TrajDETest(sub_bulk_obj, assay = "RNA", layer = "counts",
                     traj.var = "traj_DC3", latent.vars = lv, samples = non_health_idx)

for (nm in c("de_dc1", "de_dc2", "de_dc3")) {
  write.csv(get(nm), file.path(config$out_dir, paste0(nm, ".csv")))
}

# Top 20 up and 20 down genes per DC. The effect-size column name depends on the
# TrajDETest output, so it is found by pattern; check names(de_dc1) if this fails.
top_genes <- function(res, n = 20) {
  res <- as.data.frame(res)
  eff <- grep("coef|logfc|log2fc|estimate|beta", colnames(res), ignore.case = TRUE, value = TRUE)[1]
  stopifnot(!is.na(eff))
  ord <- res[order(res[[eff]], decreasing = TRUE), , drop = FALSE]
  unique(c(head(rownames(ord), n), tail(rownames(ord), n)))
}
DEG_DC1 <- top_genes(de_dc1)
DEG_DC2 <- top_genes(de_dc2)
DEG_DC3 <- top_genes(de_dc3)

# Heatmaps ordered by each DC. DC2/DC3 use the patient-only object.
sub_bulk_obj2 <- subset(sub_bulk_obj, cells = colnames(sub_bulk_obj)[non_health_idx])
cols_all <- c("#f8766d", "cornsilk", "#7cae02", "blue1", "#01bfc4", "#c77cff")

draw <- function(obj, feats, dc, cols, file) {
  obj <- ScaleData(obj, features = feats)
  p <- PlotAnnotatedHeatmap(
    object = obj, features = feats, sort.by = dc,
    anno.var = config$disease_col, anno.colors = list(cols),
    hm.colors = c("blue", "white", "red")
  )
  ggsave(file.path(config$out_dir, file), p, width = 7, height = 7)
}
draw(sub_bulk_obj,  DEG_DC1, "traj_DC1", cols_all,      "heatmap_dc1.png")
draw(sub_bulk_obj2, DEG_DC2, "traj_DC2", cols_all[-1],  "heatmap_dc2.png")
draw(sub_bulk_obj2, DEG_DC3, "traj_DC3", cols_all[-1],  "heatmap_dc3.png")
