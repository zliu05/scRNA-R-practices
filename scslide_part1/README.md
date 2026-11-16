# scSLIDE tutorial, Part 1

Scripts that follow Part 1 of the Seurat scSLIDE vignette
(https://satijalab.org/seurat/articles/scslide_vignette), sample-level analysis of
a COVID-19 scRNA-seq subset (Ahern et al. 2022, ~104k cells, 40 donors).

## Run order

Run from this folder, in R, in order:

| Script | What it does |
|---|---|
| `00_setup.R` | installs/loads packages, shared settings |
| `01_load_and_annotate.R` | load data, normalize, Azimuth annotation, HVGs |
| `02_prepare_sample_object.R` | sketching, PLS, projection, WNN graph (`use_wrapper` toggles `PrepareSampleObject()`) |
| `03_sample_level_analysis.R` | sample-level embedding, diffusion map, principal curve |
| `04_differential_expression.R` | CD14 monocyte pseudobulk trajectory DE and heatmaps |

Intermediate objects and figures go to `output/`.

## Notebook

`scslide_part1.Rmd` contains the same analysis as one R Markdown notebook. Knit it with
`rmarkdown::render("scslide_part1.Rmd")` to get the HTML report. No HTML is included yet
because the notebook hasn't been run.

## Data

Download `Raw_count_Subset_104k_Ahern_2022_Cell.rds` from the Dropbox link in the
vignette and put it in this folder.

## Status

**Not yet run.** These were written without an R installation available, so none of
the code has been executed. Expect to fix small things (argument names, metadata
column names such as `Age` and `Sex`, the effect-size column in `TrajDETest` output).
Step 01 also needs network access for `CloudAzimuth()`.
