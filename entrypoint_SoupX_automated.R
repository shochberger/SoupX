#!/usr/bin/env Rscript

library(argparse)
library(SoupX)
library(SingleCellExperiment)
library(Matrix)


parser <- ArgumentParser(description = "Run SoupX automated on a SingleCellExperiment")
parser$add_argument("--output_dir", "-o", required = TRUE, help = "Output directory provided by omnibenchmark")
parser$add_argument("--name", "-n", required = TRUE, help = "EB_all (Dataset id; used in output filenames)")
parser$add_argument("--data.sce", dest= "sce_path", default = "/data/home/hochberger/background_benchmark/data/datasets/EB_all/sc_objects/test_data_all_features.rds", help = "Path to SCE")
parser$add_argument("--cluster_col", default = "nn.clusters", help = "colData, listdata colum with cluster labels")
parser$add_argument("--set_cont", default = "Auto", help = "Contamination setting: Auto")
args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
cormat <- file.path(args$output_dir, paste0(args$name, ".SoupX_automated_corrected_counts.rds"))
perCell <- file.path(args$output_dir, paste0(args$name, ".SoupX_automated_per_cell_contamination.rds"))


run_SoupX <- function(sce_path, cluster_col, set_cont, cormat, perCell){
  # read input (SCE object)
  sce <- readRDS(sce_path)

  # is there clustering, counts, is.cell
  if (!"counts" %in% assayNames(sce)) stop("Assay 'counts' missing.")
  if (!(cluster_col %in% colnames(colData(sce)))) stop("SoupX requires pre-computed clustering!.")
  if (!"is.cell" %in% colnames(colData(sce))) stop("'is.cell' missing in colData(sce).")

  # Manually construct SoupX channel
  tod <- as(counts(sce), "dgCMatrix")
  colnames(tod) <- sce@colData@listData[["Barcode"]]
  toc <- tod[, sce$is.cell == T]

  # Adding extra metadata to SoupChannel -> Cluster-vector only for cells, human and no doublets
  clusters <- setNames(sce[, sce$is.cell == T]@colData@listData[["nn.clusters"]],
                       sce[, sce$is.cell == T]@colData@listData[["Barcode"]])

  #SoupChannel
  sc <- SoupX::SoupChannel(tod = tod, toc = toc)
  #names(clusters) <- colnames(sc$toc)
  sc <- SoupX::setClusters(sc, clusters)


   # estimate rho (contamination)
  if(set_cont == "Auto"){
    sc <- autoEstCont(sc)
  } else {
    sc <- setContaminationFraction(sc, contFrac = as.numeric(set_cont))
  }

  # correct counts, calculate corrected matrix
  corrected_counts <- adjustCounts(sc)

  #Cell specific contamination, calculate for every cell how much soup was removed
  perCell_cont <- data.frame(cell = colnames(corrected_counts),
                             cont = 1-(Matrix::colSums(corrected_counts) / Matrix::colSums(sc$toc)))

  saveRDS(corrected_counts, cormat)
  saveRDS(perCell_cont, perCell)
}
cat(sprintf(
  "[soupx_automated] Done for %s\n  - corrected: %s\n  - perCell:   %s\n",
  args$name, cormat, perCell
))

run_SoupX(args$sce_path, args$cluster_col, args$set_cont, cormat, perCell)



