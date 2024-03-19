library(Seurat)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(fs)
# function definitions
subset_count_matrix <- function(scmatrix, barcodes) {
    # helper to call purrr::map
    scmatrix[, barcodes]
}
call_writecounts <- function(sample, counts, features, data_dir) {
    # write 10x files with some checks
    # make out and sample paths
    out_path <- path(data_dir, "out")
    sample_path <- path(out_path, sample)
    # create out dir, write10xCounts does not create dirs recursively
    dir_create(sample_path)
    # check if row.names (which are unique) are the same that in features.tsv
    dedup_genes <- make.unique(features$gene_name)

    if (!identical(row.names(counts), dedup_genes)) stop("not the same")
    Matrix::writeMM(Matrix::t(counts), path(sample_path, "DGE.mtx"))
    vroom::vroom_write(features, path(sample_path, "all_genes.csv"), delim = ",")
    vroom::vroom_write(
        data.frame(bc_wells = colnames(counts), sample = sample),
        path(sample_path, "cell_metadata.csv"),
        delim = ','
    )
}
# do stuff
data_dir <- "./"
# read data in

mtx <- file.path(data_dir, "DGE.mtx")
cells <- file.path(data_dir, "cell_metadata.csv")
features <- file.path(data_dir, "all_genes.csv")
scdata <- ReadMtx(mtx = mtx, cells = cells, features = features,
               cell.column = 1, feature.column = 1, cell.sep = ",",
               feature.sep = ",", skip.cell = 1, skip.feature = 1, mtx.transpose = TRUE)


mdata <- vroom::vroom(path(data_dir,  "cell_metadata.csv"))
features <- vroom::vroom(path(data_dir,"all_genes.csv"))
# extract the counts per sample to a tibble
mdata |>
    select(bc_wells, sample) |>
    rename(barcodes = bc_wells) |>
    group_by(sample) |>
    nest() |>
    rowwise() |>
    mutate(counts = purrr::map(data, ~ subset_count_matrix(scdata, .x))) -> matrices
# write the counts to 10x files
matrices %>%
    select(sample, counts) %>%
    pwalk(call_writecounts, features, data_dir)
