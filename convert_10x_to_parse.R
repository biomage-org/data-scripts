library(Seurat)
library(dplyr)
library(fs)

# Set the input directory containing sample data
input_dir <- "./"
samples <- list.files(input_dir)

# Loop through each sample and convert data
for (sample in samples) {
  message("Converting sample ", sample)
  
  sample_dir <- file.path(input_dir, sample)
  sample_fpaths <- list.files(sample_dir)
  annot_fpath <- file.path(sample_dir, "features.tsv.gz")
  
  annotations <- read_10x_annotations(annot_fpath, sample)
  counts <- Seurat::Read10X(sample_dir, gene.column = 1, unique.features = TRUE)
  if (is(counts, "list")) {
    slot <- "Gene Expression"
    if (!(slot %in% names(counts))) slot <- names(counts)[1]
    counts <- counts[[slot]]
  }
  
  out_path <- path(input_dir, "converted", sample)
  if (!dir_exists(out_path)) dir_create(out_path)
  
  # write Parse files
  Matrix::writeMM(counts, path(out_path, "DGE.mtx"))
  vroom::vroom_write(annotations, path(out_path, "all_genes.csv"), delim = ",")
  vroom::vroom_write(
    data.frame(bc_wells = colnames(counts), sample = sample),
    path(out_path, "cell_metadata.csv"),
    delim = ","
  )
  
  message("Finished converting sample ", sample)
}

read_10x_annotations <- function(annot_fpath, sample) {
  gene_column <- 1
  
  annot <- read.delim(annot_fpath, header = FALSE)
  
  # Remove features that are not "Gene Expression"
  if (ncol(annot) > 2 && length(grep("Gene Expression", annot$V3)) > 0) {
    annot <- annot %>% dplyr::filter(V3 == "Gene Expression")
  }
}
