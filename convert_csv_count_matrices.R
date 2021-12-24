library(data.table)
library(Matrix)
library(fs)

data_dir <- path_abs("./data/GSE150115_RAW/")


dir_create(path(data_dir, "out2"))

for (file in path_file(dir_ls(data_dir, glob = "*.gz"))) {
  dirname <- gsub("\\.fastq.*", "",file)
  
  dt <- fread(path(data_dir, file))
  dt <- as.matrix(dt, rownames = 1)
  dt <- Matrix::Matrix(dt, sparse = T)
  
  DropletUtils::write10xCounts(path(data_dir, "out2", dirname), x = dt)
  rm(dt)
  gc()
  
}
