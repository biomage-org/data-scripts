#Convert a single parse sample to 10x.
#Asumes the samples have been demultiplexed
expected_files <- list("all_genes.csv","cell_metadata.csv","DGE.mtx")
convert_parse_sample_to_10x <- function(sample_path,sample_out_path){

  dir.create(sample_out_path)

  for (i in expected_files) {
    if (!(i %in% list.files(sample_path))) {
      message("Cant find ",i," in sample path.")
      return()
    }
  }

  all_genes <- read.csv(paste0(sample_path,"/all_genes.csv"))
  cell_metadata <- read.csv(paste0(sample_path,"/cell_metadata.csv"))
  mtx <- Matrix::readMM(paste0(sample_path,"/DGE.mtx"))

  all_genes[,3] <- "Gene Expression"
  write.table(all_genes, file=paste0(sample_out_path,"/features.tsv"), row.names=FALSE,col.names=FALSE, sep="\t")

  cell_metadata  <- as.data.frame(cell_metadata[,"bc_wells"])
  write.table(cell_metadata, file=paste0(sample_out_path,"/barcodes.tsv"), row.names=FALSE,col.names=FALSE, sep="\t")

  Matrix::writeMM(Matrix::t(mtx), paste0(sample_out_path,"/matrix.mtx"))
}

samples <- list("Gut HO", "Gut NO", "Lung HO", "Lung NO")
dir.create("./out")
for (sample in samples){
  sample_path <- paste0("./",sample,"/DGE_filtered")
  out_path <- paste0("./out/",sample)
  convert_parse_sample_to_10x(sample_path, out_path)
}

