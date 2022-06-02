source("filter_empty_drops.R")
source("load_user_files_multi.R")
source("utils.R")

# Demultiplex using HTOs 
hto_demux <- function() {
sampleIds <- list.files("./input/")
config <- list(input = list(type = "10x"), samples = sampleIds)

# Function load_user_files_multi doesn't take only the "Gene Expression" slot, but also the "Antibody Capture" slot (commented lines 76-81 of gem2s-2-load_user_files.R)
# Load unfiltered data
data_unf <- read_10x_with_ac(config, "./input")
# Load filtered data
data_fil <- read_10x_with_ac(config, "./out")

annot <- data_fil$output$annot
counts_list_unf <- data_unf$output$counts_list
counts_list_fil <- data_fil$output$counts_list
# To allow compatibility with call_writecounts
colnames(annot) <- c("X1", "X2", "original_name")

# Create Seurat object with Gene Expression and Antibody Capture for each sample
for (sample in sampleIds) {
  ge <- "Gene Expression"
  ac <- "Antibody Capture"
  counts_unf_ge <- counts_list_unf[[sample]][[ge]]
  counts_unf_ac <- counts_list_unf[[sample]][[ac]]
  counts_fil <- counts_list_fil[[sample]]

  scdata <- Seurat::CreateSeuratObject(counts = counts_unf_ge)
  scdata[["HTO"]] <- Seurat::CreateAssayObject(counts = counts_unf_ac)
  # Subset Seurat object by filtered barcodes
  scdata <- subset(scdata, cells = colnames(counts_fil))
  
  # Run HTODemux
  scdata <- Seurat::HTODemux(scdata, assay = "HTO", positive.quantile = 0.99)
  call_hto_demux(scdata,sample,annot)
  }

  call_hto_demux <- function(obj,sample,features){
  # Exclude cells classified as "Negative" and "Doublet"
    for(name in unique(obj$hash.ID))
    {
      if(name == "Negative" || name == "Doublet")
      {
        next
      }
    # Subset by HTO classification and write sample files
      obj.sub <- subset(x = obj, subset = hash.ID == name)
      counts <- obj.sub@assays$RNA@data
      sample_path <- paste0(sample,"-",name)
      message("Writing demultiplexed files for sample ",sample,"-",name)
      call_writecounts(sample_path, counts, features, "./hto_dem")
    }
  }
}

hto_demux()
