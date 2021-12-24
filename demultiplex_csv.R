# demultiplex csv count matrices
# requirements: samples encoded in barcode name

library(data.table)
library(Matrix)
library(fs)

# maybe implement as list, to avoid moving 2 objects around
# d <- list()
# d[[1]] <- fread(path(data_dir, "GSE145980_zebrafish_nonCM_wildtype_scRNA_raw_count_matrix.csv.gz"))
# d[[2]] <- preprocess_data(d[[1]])

#' clean original data.table CSV column names
#' 
#' Removes sample information from column names. It modifies in place!
#'
#' @param dt 
#' @param sample_barcode_tab 
#'
#' @return
#' @export
#'
#' @examples
clean_dt_colnames <- function(dt, clean_barcodes) {
  setnames(dt, base::colnames(dt), clean_barcodes)
}

#' make sample <-> barcode table
#'
#' Extracts sample name from "sample_barcode" encoded column names in GEO csv. 
#' Creates table with barcode - sample association.
#' Users should manually check if the regex is correct for the particular dataset
#' being demultiplexed.
#'
#' @param dt data.table original GEO dataset
#' @param sample_regex chr regex to parse column names for sample and barcodes
#'
#' @return data.table
#' @export
#'
#' @examples
make_sample_barcode_tab <- function(dt, sample_regex = "(.*)_([[:upper:]]+)") {
  samp_bc <- colnames(dt)

  sample_names <- gsub(sample_regex, "\\1", samp_bc)
  barcodes <- gsub(sample_regex, "\\2", samp_bc)

  clean_dt_colnames(dt, barcodes)

  # first var in dt is the gene_names var (data.tables don't have rownames)
  data.table(
    sample = sample_names[-1],
    barcode = barcodes[-1]
  )
}


#' Create list of barcodes in samples
#'
#' @param sample_barcode_tab data.table sample/barcode table
#'
#' @return list one element per sample, with every barcode in sample
#' @export
#'
#' @examples
list_barcodes_in_sample <- function(sample_barcode_tab) {
  # nest each barcode group to separate data.table
  nested_sample_dt <- sample_barcode_tab[, .(bc_list = list(.SD)), by = sample]

  # convert nested data table to list
  lapply(nested_sample_dt[["bc_list"]], unlist)
}


#' subset data.table
#'
#' Subsets cleaned (clean_dt_colnames) data.table, provided character vector of 
#' barcodes in sample. 
#' Helper function to simplify lapply calls. 
#'
#' @param dt data.table cleaned count csv
#' @param columns character vector 
#'
#' @return data.table subsetted data.table
#' @export
#'
#' @examples
sub_dt <- function(columns, dt) {
  # subset a data table by character vector, to ease lapply
  columns <- c("V1", columns)
  dt[, ..columns]
}


#' create list containing sparse matrix for each sample 
#'
#' CAUTION! only run if enough RAM is available. Otherwise, use the global 
#' environment.
#' 
#' transforms list of data tables per sample into list of sparse matrices.
#'
#' @param sample_list list list of barcodes per sample
#' @param dt data.table original csv data
#'
#' @return list each element a sample sparse matrix.
#' @export
#'
#' @examples
create_sparse_matrices <- function(sample_list, dt) {
  gc(verbose = F)
  # subset the original count data.table
  sub_dts <- lapply(sample_list, sub_dt, dt)
  rm(list = as.character(substitute(dt)), envir = .GlobalEnv)
  gc(verbose = F)

  # convert each subsetted count data.table to count matrix
  counts <- lapply(sub_dts, as.matrix, rownames = "V1")
  rm(sub_dts)
  gc(verbose = F)

  # convert each count matrix to sparse matrices
  sparse_counts <- lapply(counts, Matrix, sparse = T)
  rm(counts)
  gc(verbose = F)

  sparse_counts
}


#' export demultiplexed data
#' 
#' Creates `data_dir/out` folder, and exports 10X files in a folder per sample.
#' 
#' @param sample_dt data.table sample <-> barcode table
#' @param sparse_matrix_list list each element a sample sparse matrix.
#' @param data_dir chr root dir to export
#'
#' @return NULL
#' @export
#'
#' @examples
export_demultiplexed_data <- function(sample_dt, sparse_matrix_list, data_dir) {
  dir_create(path(data_dir, "out"))

  nested_sample_dt <- sample_dt[, .(bc_list = list(.SD)), by = sample]
  
  for (row in 1:nrow(nested_sample_dt)) {
    fname <- path(data_dir, "out", nested_sample_dt[row][["sample"]])

    # unnest barcodes in sample
    expected_barcodes_in_sample <- nested_sample_dt[row, bc_list[[1]]][["barcode"]]

    if (!identical(expected_barcodes_in_sample, colnames(sparse_matrix_list[[row]]))) {
      stop("not the same barcodes")
    }

    DropletUtils::write10xCounts(fname,
      sparse_matrix_list[[row]],
      version = "3"
    )
  }
}


################################################################################
# do things

data_dir <- path_abs("./data")

d <- fread(path(data_dir, "GSE145980_zebrafish_nonCM_wildtype_scRNA_raw_count_matrix.csv.gz"))
sample_tab <- make_sample_barcode_tab(d)


# these steps remove the previous object for memory reasons. couldn't manage 
# to abstract them away without hitting memory limit.

gc()
# subset the original count data.table
dt_subset <- lapply(list_barcodes_in_sample(sample_tab), sub_dt, d)
rm(d)
gc()

# convert each subsetted count data.table to count matrix
counts <- lapply(dt_subset, as.matrix, rownames = "V1")
rm(dt_subset)
gc()

# convert each count matrix to sparse matrices
sparse_counts <- lapply(counts, Matrix, sparse = T)
rm(counts)
gc()


export_demultiplexed_data(sample_tab, sparse_counts, data_dir)
