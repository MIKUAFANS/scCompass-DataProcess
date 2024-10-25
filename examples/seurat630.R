#!/usr/bin/env Rscript
# usage: Rscript --vanilla seurat630.R --outdir . --data_dir "./GSM6508470/outs/filtered_feature_bc_matrix" --min_cells 3 --min_features 200 --min_rna 200 --percent_mt_max 15

# adjust to 630 new pipeline

# ---install---
# pkgs <- c("Seurat", "getopt")
# if (!require(pkgs, quietly = TRUE))
#     install.packages(pkgs)

# ---getopt---
library(getopt)

spec <- matrix(c(
  "verbose", "v", 2, "integer", "seurat filter", 
  "help", "h", 0, "logical", "help",
  "data_dir", "d", 1, "character", "Input directory",
  "min_cells", "cmin", 1, "integer", "Include features detected in at least this many cells",
  "min_features", "fmin", 1, "integer", "Include cells where at least this many features are detected",
  "min_rna", "rmin", 1, "integer", "Include rnas where at least this many features are detected",
  "percent_mt_max", "mmax", 1, "integer", "Include cells where at most this percentage of mitochondrial genes are detected",
  "outdir", "o", 1, "character", "Set output directory"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if(!is.null(opt$help)){
  cat(paste(getopt(spec, usage = TRUE), "\n"))
  q()
}

if(is.null(opt$min_cells)){opt$min_cells = 3}
if(is.null(opt$min_features)){opt$min_features = 200}
if(is.null(opt$min_rna)){opt$min_rna = 200}
if(is.null(opt$percent_mt_max)){opt$percent_mt_max = 15}

# ---analysis---
library(Seurat)

seurat_filter <- function(data_dir = data_dir, min_cells = min_cells, min_features = min_features, min_rna = min_rna, percent_mt_max = percent_mt_max, outdir = outdir){
  # ---input---
  wd <- getwd()
  project <- strsplit(data_dir, "/")[[1]]
  project <- project[length(project)-2]

  dat <- Read10X(data.dir = data_dir)
  seurat <- CreateSeuratObject(counts = dat, project = project, min.cells = min_cells, min.features = min_features)

  setwd(outdir)

  # ---filter mt percentage---
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  seurat_filter <- subset(x = seurat, subset = nFeature_RNA > min_rna & percent.mt < percent_mt_max)

  # ---output---
  write.csv(seurat_filter@assays$RNA@counts, paste0(project, "_count.csv"), row.names = T, quote = F)
}

seurat_filter(data_dir = opt$data_dir, min_cells = opt$min_cells, min_features = opt$min_features, min_rna = opt$min_rna, percent_mt_max = opt$percent_mt_max, outdir = opt$outdir)


