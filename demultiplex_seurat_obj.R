# demultiplex Seurat object
# requirements: samples in seurat object metadata, column named "sample"

library(dplyr)
library(DropletUtils)
library(fs)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

################################################################################
# functions

create_sample_key <- function(seurat_obj) {
    # given a seurat object with sample metadata,
    # return tibble with sample - barcode reference

    samples <- seurat_obj$sample
    tibble(sample = samples, barcode = names(samples))
}

subset_sample <- function(seurat_obj, sample_name) {
    # subset seurat object by sample name

    # in case there's a sample named "NA"
    if (is.na(sample_name)) {
        subset(seurat_obj,
            cells = names(
                seurat_obj$sample[which(is.na(seurat_obj$sample))]
            )
        )
    } else {
        subset(seurat_obj,
            cells = names(
                seurat_obj$sample[which(seurat_obj$sample == sample_name)]
            )
        )
    }
}

get_count_matrix <- function(seurat_obj) {
    # wrapper to get count matrix
    seurat_obj@assays$RNA@counts
}

call_writecounts <- function(sample, counts, features, data_dir) {
    # write 10x files with some checks

    # make out and sample paths
    out_path <- path(data_dir, "out")
    sample_path <- path(out_path, sample)

    # create out dir, write10xCounts does not create dirs recursively
    dir_create(out_path)

    # check if row.names (which are unique) are the same that in features.tsv
    dedup_genes <- make.unique(features$X2)
    if (!identical(row.names(counts), dedup_genes)) stop("not the same")

    DropletUtils::write10xCounts(sample_path,
        counts,
        gene.id = features$X1,
        gene.symbol = features$X2,
        version = "3"
    )
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

# load rds file, containing pre-processed seurat obj
sc <- read_rds("./data/Yutao (HMS)/sc.Rds")

# read features file, important for ensemblIDs
features <- read_tsv(path(chromium_dir, "features.tsv.gz"), col_names = F)

# create sample - barcode table
sample_key <- create_sample_key(sc, "\\d+", ".*")

# get original seurat object count matrices
samples %>%
    group_by(sample) %>%
    nest() %>%
    rowwise() %>%
    mutate(counts = map(
        sample,
        ~ get_count_matrix(subset_sample(sc, .x))
    )) -> matrices

# write count matrices to files
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
