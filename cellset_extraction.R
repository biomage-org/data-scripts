library(Seurat)
library(dplyr)
source("utils.R")

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
#' @param cellset_type int number corresponding to cellset type, see details
#' @param cellset_number int cellset position (cluster number, sample number, etc.)
#'
#' @return numeric vector of cell indices
#' @export
#'
get_cellset_cell_ids <-
  function(cellsets, cellset_type, cellset_number) {
    if (cellset_type == 2 & length(cellsets[["cellSets"]]) < 3) {
      # this check might not be necessary.
      # Added it bc I think I have seen cellset objects with two elements here
      stop("There are no scratchpad cell sets in this object")
    }

    as.numeric(cellsets[["cellSets"]][[cellset_type]][["children"]][[cellset_number]][["cellIds"]])
}


#' Subset seurat object given cell_ids
#'
#' Cell Ids are an internal representation used by Cellenics. They are stored in
#' the `meta.data` slot of the Seurat object. This function subsets a Seurat
#' object given a vector of cell ids.
#'
#' @param processed_matrix Seurat rds object downloaded from S3
#' @param cell_ids numeric vector of cell ids to extract
#'
#' @return subsetted Seurat object
#' @export
#'
subset_seurat_cell_ids <- function(processed_matrix, cell_ids) {
  barcodes_to_keep <- processed_matrix@meta.data %>%
    dplyr::filter(cells_id %in% cell_ids) %>% # this is not a typo
    rownames()

  processed_matrix[, barcodes_to_keep]
}



#' extract cellset from a seurat object
#'
#' Given experiment ID and cellset identifiers (type and number), this function
#' downloads the required files, subsets and returns the seurat object.
#'
#' @param experiment_id character experiment ID
#' @param cellset_type int cellset type
#' @param cellset_number int cellset number
#'
#' @return subsetted seurat object
#' @export
#'
extract_cellset <- function(experiment_id, cellset_type, cellset_number) {

  # download stuff
  download_cellset_file(experiment_id)
  download_processed_matrix(experiment_id)

  # load stuff
  cellsets <- jsonlite::read_json(file.path(experiment_id, "cellsets.json"))
  processed_matrix <- readRDS(file.path(experiment_id, "r.rds"))

  # do stuff
  cell_ids <- get_cellset_cell_ids(cellsets, cellset_type, cellset_number)
  subset_seurat_cell_ids(processed_matrix, cell_ids)
}

