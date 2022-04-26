# set working directory and load rds object
  setwd("/Your/Directory")
  data <- readRDS("1867-counts_cells_cohort2.rds")

# load required package
  library("DropletUtils")
# comment
# look at data structure
  str(data)
# we can see that this is a dgCMatrix with gene symbols as rownames and cell barcodes as colnames
# let's look at the firs cell barcode as an example (BIOKEY_13_Pre_AAACCTGCAACAACCT-1)
  colnames(data)[1]
# in this particular case, cell barcodes consist of a prefix (BIOKEY_13_Pre_) and a sequence of bases (AAACCTGCAACAACCT-1)
# the prefixes are sample names, so we'll use them to demultiplex the data

# demultiplex data based on barcodes prefixes and export as 10X files
  # use a regular expression to extract prefixes
  data.pfx <- gsub("(.+)_[A-Z]+-1$", "\\1", colnames(data), perl=TRUE)
  # get unique sample names
  data.samples <- unique(data.pfx)
  # check  sample names (pay attention if are using this script to process a different dataset as the regular expression may need to be modified depending on the specific colnames)
  head(data.samples)
  tail(data.samples)

  # export as 10X files that can be directly uploaded to Cellenics
  # define the function
    # the function creates a subdirectory named "demultiplexed" inside the current working directory, and save 10X data for each sample in different subfolders
    # if a folder named "demultiplexed" already exists, it will stop and return an error to avoid overwriting files
  demultiplex_convert_to_10x <- function(obj, samples) {
          if(!dir.exists(file.path(getwd(), "demultiplexed"))) {
          dir.create(file.path(getwd(), "demultiplexed"))
        } else {
          print("WARNING! A demultiplexed directory already exists")
          return()
        }
        for (i in 1:length(samples)) {
        print(paste0("Converting sample ", samples[i]))
        DropletUtils::write10xCounts(path = paste0(getwd(),"/demultiplexed/",samples[i]), x = obj[,grep(paste0("^",samples[i],".*"),colnames(obj))], type = "sparse", version="3")
        }
  }

  # run the function
  demultiplex_convert_to_10x(obj = data, samples = data.samples)




