library(jsonlite)
library(Seurat)


download_cellset_file <- function(experiment_id) {
  remote_path <-
    paste("s3:/",
          "cell-sets-production",
          experiment_id,
          sep = "/")
  
  local_path <- paste0(experiment_id, "_cellset.json")
  args <- c("s3", "cp", remote_path, local_path)
  system2("aws", args)
}


download_processed_matrix <- function(experiment_id) {
  remote_path <-
    paste0(paste("s3:/", "processed-matrix-production",
                 experiment_id,
                 sep = "/"),
           "/")
  local_path <- paste0(experiment_id, "_r.rds")
  
  args <- c("s3", "cp", remote_path, local_path, "--recursive")
  system2("aws", args)
}


