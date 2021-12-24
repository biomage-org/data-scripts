source("csv_helpers.R")

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
