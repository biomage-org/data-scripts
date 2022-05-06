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


#' Extract cell ids of cells in a cellset
#'
#' Extracts the cell ids of cells that belong to a cellset from the cellset object as
#' a numeric vector.
#' 
#' The `cellset_type` argument is an int that specifies the cellset class:
#'
#' \itemize{
#'   \item `1 = louvain`
#'   \item `2 = scratchpad`
#'   \item `3 = samples`
#'   \item `>= 4 = metadata tracks`
#'}
#' 
#' The `cellset_number` refers to a particular cellset within a class. 
#' 
#' So, to extract the third scratchpad cellset, we would set \code{cellset_type = 2},
#' \code{cellset_number = 3}.
#'
#' @param cellsets list parsed json cellset object
#' @param cellset_type int number corresponding to cellset type, see description
#' @param cellset_number int cellset position (cluster number, sample number, etc.)
#'
#' @return numeric vector of cell indices
#' @export
#'
get_cellset_cell_ids <-
  function(cellsets, cellset_type, cellset_number) {
    if (type == 2 & length(cellsets[["cellSets"]]) < 3) {
      # this check might not be necessary.
      # Added it bc I think I have seen cellset objects with two elements here
      stop("There are no scratchpad cell sets in this object")
    }
    
    as.numeric(cellsets[["cellSets"]][[cellset_type]][["children"]][[cellset_number]][["cellIds"]])
  }
