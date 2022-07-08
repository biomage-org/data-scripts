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


read_10x_with_ac <- function(config, input_dir) {
  counts_list <- list()
  annot_list <- list()
  
  samples <- config$samples
  
  
  for (sample in samples) {
    sample_dir <- file.path(input_dir, sample)
    sample_fpaths <- list.files(sample_dir)
    annot_fpath <- file.path(sample_dir, "features.tsv.gz")
    
    message("\nSample --> ", sample)
    message(
      "Reading files from ",
      sample_dir,
      " --> ",
      paste(sample_fpaths, collapse = " - ")
    )
    
    counts <- Seurat::Read10X(sample_dir, gene.column = 1)
    
    annot <- read.delim(annot_fpath, header = FALSE)
    
    # Equalizing number of columns in case theres no Gene Expression column
    annot <- annot[, c(1, 2)]
    
    message(
      sprintf(
        "Sample %s has %s genes and %s droplets.",
        sample, nrow(counts), ncol(counts)
      )
    )
    
    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annot
  }
  
  annot <- format_annot(annot_list)
  
  return(list(counts_list = counts_list, annot = annot))
}


format_annot <- function(annot_list) {
  annot <- unique(do.call("rbind", annot_list))
  colnames(annot) <- c("input", "name")
  
  message("Deduplicating gene annotations...")
  
  # add ENSEMBL ID for genes that are duplicated (geneNameDuplicated-ENSEMBL)
  # original name kept in 'original_name' column
  gname <- annot$name
  annot$original_name <- gname
  is.dup <- duplicated(gname) | duplicated(gname, fromLast = TRUE)
  
  # We need to convert the gene inputs from _ to - bc when we create the Seurat
  # object we do this, and the match would return NA values if any
  # of the inputs still has _.
  annot$input <- gsub("_", "-", annot$input)
  annot$name[is.dup] <- paste(gname[is.dup], annot$input[is.dup], sep = " - ")
  
  annot <- annot[!duplicated(annot$input), ]
  
  rownames(annot) <- annot$input
  return(annot)
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


#' Extract counts from Seurat object
#'
#' @param seurat_obj SeuratObject
#'
#' @return
#' @export
#'
extract_counts <- function(seurat_obj) {
  seurat_obj@assays$RNA@counts
  
}


#' Carefully extract names from gene_annotations
#'
#' Extracts gene names from the annotations data.frame taking care of matching
#' the input column.
#'
#' @param seurat_obj SeuratObject
#'
#' @return
#' @export
#'
extract_gene_symbols <- function(seurat_obj) {
  idx <-
    match(rownames(seurat_obj),
          seurat_obj@misc$gene_annotations$input)
  seurat_obj@misc$gene_annotations$name[idx]
}

#' add sample names to seurat object
#'
#' The seurat objects stored in S3 do not contain the sample names, which
#' complicates processing user support requests. This function adds them using the
#' sample_mapping.json file, which can be downloaded using `biomage experiment download`
#'
#' @param seurat_obj
#' @param sample_mapping_json_path
#'
#' @return
#' @export
#'
add_sample_names <- function(seurat_obj, sample_mapping_json) {
  sample_map <- jsonlite::read_json(sample_mapping_json) %>%
    tibble::enframe(x = ., name = "sample_name", value = "sample_id") %>%
    tidyr::unnest(sample_id)
  
  sample_names <- seurat_obj@meta.data %>%
    dplyr::left_join(sample_map, by = c("samples" = "sample_id")) %>%
    dplyr::pull(sample_name)
  
  Seurat::AddMetaData(seurat_obj, sample_names, col.name = "sample_name")
}
