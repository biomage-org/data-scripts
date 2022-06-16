# install and load required libraries
install.packages("Matrix", "remotes")
remotes::install_github("MarioniLab/DropletUtils")

library(Matrix, quietly=T)
library(DropletUtils, quietly=T)

species <- commandArgs(trailingOnly = TRUE)

# set directory
kbcount_out_path <- "./"
setwd(fold_path)

# convert kb count output to 10x
for (fold in list.files(pattern="*_kbcount_output")) {
  base_path <- file.path(fold_path, fold)
  matrix_path <- file.path(base_path, "counts_unfiltered", "cells_x_genes.mtx")
  genes_path <- file.path(base_path, "/counts_unfiltered/cells_x_genes.genes.txt")
  barcodes_path <- file.path(base_path, "/counts_unfiltered/cells_x_genes.barcodes.txt")

  raw_mtx <- as(t(readMM(matrix_path)), "CsparseMatrix") # load mtx and transpose it
  rownames(raw_mtx) <- read.delim(genes_path, header = F)[,1] # attach genes
  colnames(raw_mtx) <- read.delim(barcodes_path, header = F)[,1] # attach barcodes

  t2g_path <- paste0("./", species, "-ref/", species, "-reft2g.txt")
  t2g <-  unique(read.delim(t2g_path, header=F)[,2:3]) # load t2g file
  t2g <- data.frame(t2g[,2], row.names = t2g[,1])
  gene_sym <- t2g[as.character(rownames(raw_mtx)),1] # get symbols for gene ids

  name = gsub("_kbcount_output","",fold)
  write10xCounts(paste0(fold_path,"/",name), gene.symbol = gene_sym, raw_mtx, overwrite = T, type = "sparse", version="3") # write results
}
