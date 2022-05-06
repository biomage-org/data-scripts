library(fs)

call_writecounts <- function(sample, counts, features, data_dir) {
    # write 10x files with some checks

    # make out and sample paths
    out_path <- path(data_dir, "out")
    sample_path <- path(out_path, sample)

    # create out dir, write10xCounts does not create dirs recursively
    dir_create(out_path)

    # check if row.names (which are unique) are the same that in features.tsv
    dedup_genes <- make.unique(features$X1)
    if (!identical(row.names(counts), dedup_genes)) stop("not the same")

    DropletUtils::write10xCounts(sample_path,
                                 counts,
                                 gene.id = features$X1,
                                 gene.symbol = features$X2,
                                 version = "3"
    )
}

#' Download cell sets file
#'
#' Downloads cellset file from S3, given an experiment ID.
#' Uses aws-cli to do so, which must be configured and working.
#'
#' @param experiment_id character
#'
#' @return NULL
#' @export
download_cellset_file <- function(experiment_id) {
  remote_path <-
    paste("s3:/",
          "cell-sets-production",
          experiment_id,
          sep = "/")

  local_path <- file.path(experiment_id, "cellsets.json")
  args <- c("s3", "cp", remote_path, local_path)
  if (!file.exists(local_path)) {
    system2("aws", args)
  }
}


#' Download processed matrix
#'
#' Downloads processed matrix RDS file from S3, given an experiment ID.
#' Uses aws-cli to do so, which must be configured and working.
#'
#' @param experiment_id character
#'
#' @return NULL
#' @export
download_processed_matrix <- function(experiment_id) {
  remote_path <-
    paste0(paste("s3:/", "processed-matrix-production",
                 experiment_id,
                 sep = "/"),
           "/")
  local_path <- experiment_id

  args <- c("s3", "cp", remote_path, local_path, "--recursive")

  if (!file.exists(file.path(local_path, "r.rds"))) {
    system2("aws", args)
  }
}
