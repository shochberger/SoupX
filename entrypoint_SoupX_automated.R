#!/usr/bin/env Rscript

# Omnibenchmark wrapper for SoupX: loads an SCE (.sce.rds), builds a SoupChannel
# (tod=all droplets, toc=only cells via is.cell), uses colData$nn.clusters, estimates
# rho via Auto or applies a fixed --set_cont, runs adjustCounts(), and writes
# {name}.soupX_corrected.rds + {name}.soupX_percell.rds to --output_dir.

library(argparse)
library(SoupX)
library(SingleCellExperiment)
library(Matrix)
library(tibble)

parser <- ArgumentParser(description = "Run SoupX automated on a SingleCellExperiment")
parser$add_argument("--output_dir", "-o", required = TRUE, help = "Output directory by OB")
parser$add_argument("--name", "-n", required = TRUE, help = "name of the dataset")
parser$add_argument("--data.sce", dest= "sce_path", required = TRUE, help = "Path to SCE from data stage")
parser$add_argument("--set_cont", default = "Auto", help = "Contamination setting: Auto")
args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
cormat <- file.path(args$output_dir, paste0(args$name, ".soupX_corrected.rds"))
perCell <- file.path(args$output_dir, paste0(args$name, ".soupX_percell.rds"))

run_SoupX <- function(sce_path, set_cont, cormat, perCell){
  # read input (SCE object)
  sce <- readRDS(sce_path)

# sanity checks (fails fast if required )
if (!"counts" %in% assayNames(sce)) stop("Assay 'counts' missing.")
if (!("nn.clusters" %in% colnames(colData(sce)))) {
  stop("Clustering colum 'nn.clusters' not found in colData(sce).")
}
if (!"is.cell" %in% colnames(colData(sce))) stop("'is.cell' missing in colData(sce).")

## Manually construct SoupChannel object by providing table of counts and table of droplets
tod <- as(counts(sce), "dgCMatrix")
colnames(tod) <- sce@colData@listData[["Barcode"]]
toc <- tod[, sce$is.cell == TRUE, drop = FALSE]

# Adding extra metadata to SoupChannel -> Cluster-vector only for cells, human and no doublets
clusters <- setNames(sce[, sce$is.cell == T]@colData@listData[["nn.clusters"]],
                      sce[, sce$is.cell == T]@colData@listData[["Barcode"]])

# Create SoupChannel object
sc <- SoupX::SoupChannel(tod = tod, toc = toc)
sc <- SoupX::setClusters(sc, clusters)

# Estimate rho (contamination fraction): "Auto" via autoEstCont, else manually specify rho via setContaminationFraction
if(set_cont == "Auto"){
  sc <- autoEstCont(sc)
} else {
  sc <- setContaminationFraction(sc, contFrac = as.numeric(set_cont))
}

# Correct counts
corrected_counts <- adjustCounts(sc)

#Cell specific contamination, calculate for every cell how much soup was removed
perCell_cont <- data.frame(cell = colnames(corrected_counts),
                            cont = 1-(Matrix::colSums(corrected_counts) / Matrix::colSums(sc$toc)))

# Write outputs (OB will collect via configured paths)
saveRDS(corrected_counts, cormat)
saveRDS(perCell_cont, perCell)
}
cat(sprintf(
  "[soupx_automated] Done for %s\n  - corrected: %s\n  - perCell:   %s\n",
  args$name, cormat, perCell
))

run_SoupX(args$sce_path, args$set_cont, cormat, perCell)


