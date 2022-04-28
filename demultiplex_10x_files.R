# demultiplex features barcodes matrix files
# requirements: sample key encoded in barcode names, ex ACTCGCATCGACATC-1
# reads to dgcMatrix, and subset directly

library(dplyr)
library(fs)
library(purrr)
library(tibble)
library(tidyr)
library(stringr)
source("utils.R")

# functions

create_sample_key <- function(scmatrix, sample_pattern, barcode_pattern = ".*") {
    # careful with non-unique barcodes, it might be better to keep
    # barcodes with original names including sample

    column_names <- colnames(scmatrix)
    samples <- as.character(stringr::str_match(column_names, sample_pattern))
    barcodes <- as.character(stringr::str_match(column_names, barcode_pattern))

    tibble(sample = samples, barcodes = barcodes)
}

subset_count_matrix <- function(scmatrix, barcodes) {
    # helper to call purrr::map
    scmatrix[, barcodes]
}

read_demultiplexed <- function(sample_name, data_dir) {
    # read the 10x files from demultiplexed dir

    path <- file.path(data_dir, "out", sample_name)
    Seurat::Read10X(path)
}

################################################################################
# do things

# set directories
data_dir <- path()
chromium_dir <- path(data_dir, "<10x_folder>")

# read features file, important for ensemblIDs
features <- read_tsv(path(chromium_dir, "features.tsv.gz"), col_names = F)

# create sample - barcode table
sample_key <- create_sample_key(sc, "\\d+", ".*")

# subset sparse matrix by barcodes
sample_key %>%
    group_by(sample) %>%
    nest() %>%
    rowwise() %>%
    mutate(counts = purrr::map(data, ~ subset_count_matrix(sc, .x))) -> matrices

# write subsetted matrices to 10X format
matrices %>%
    select(sample, counts) %>%
    pwalk(call_writecounts, features, data_dir)



################################################################################
# tests

matrices %>%
    mutate(read = map(sample, read_demultiplexed, data_dir)) -> matrices

# test number of expressed genes
count_genes <- function(sparse_matrix) {
    sum(Matrix::rowSums(sparse_matrix) > 0)
}

# compare with the original count matrices
matrices %>%
    ungroup() %>%
    select(sample, counts, read) %>%
    mutate(
        all_equal = all.equal(counts, read),
        identical = identical(counts, read)
    ) %>%
    rowwise() %>%
    mutate(
        ori_n_genes = count_genes(counts),
        my_n_genes = count_genes(read)
    )
