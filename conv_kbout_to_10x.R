# load libraries
library(Matrix, quietly=T)
library(DropletUtils, quietly=T)

# set directory
fold_path <- "/Users/sara/Documents/Data/kaitilin"
setwd(fold_path)

# convert kb count output to 10x
for (fold in list.files(pattern="*_kbcount_output"))
{
  raw_mtx <- as(t(readMM(paste0(fold_path,"/",fold,'/counts_unfiltered/cells_x_genes.mtx'))), 'CsparseMatrix') # load mtx and transpose it
  rownames(raw_mtx) <- read.csv(paste0(fold_path,"/",fold,'/counts_unfiltered/cells_x_genes.genes.txt'), sep = '\t', header = F)[,1] # attach genes
  colnames(raw_mtx) <- read.csv(paste0(fold_path,"/",fold,'/counts_unfiltered/cells_x_genes.barcodes.txt'), header = F, sep = '\t')[,1] # attach barcodes

  t2g <-  unique(read.csv('/Users/sara/Documents/Tools/kallisto-bustools/mouse-ref/mouse-reft2g.txt', sep = '\t', header=F)[,2:3]) # load t2g file
  t2g <- data.frame(t2g[,2], row.names = t2g[,1])
  gene_sym <- t2g[as.character(rownames(raw_mtx)),1] # get symbols for gene ids

  name = gsub("_kbcount_output","",fold)
  write10xCounts(paste0(fold_path,"/",name), gene.symbol = gene_sym, raw_mtx, overwrite = T, type = "sparse", version="3") # write results
}

