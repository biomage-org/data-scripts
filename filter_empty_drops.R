library(ggplot2)
source("utils.R")

filter_empty_drops <- function() {
  # Directories set up
  unlink("out", recursive = TRUE)
  unlink("qc_plots/emptyDrops")
  dir.create("qc_plots", showWarnings = FALSE)
  dir.create("qc_plots/emptyDrops", showWarnings = FALSE)
  dir.create("out", showWarnings = FALSE)


  sampleIds <- list.files("./input/")
  input_dir <- "./input"
  config <- list(input = list(type = "10x"), samples = sampleIds)
  filtered_samples <- list()

  # Run empty drops
  data <- pipeline::load_user_files(list(), list(), list(config = config), "./input")
  data <- pipeline::run_emptydrops(list(), list(), data$output)

  edrops <- data$output$edrops
  annot <- data$output$annot
  counts_list <- data$output$counts_list
  #To allow compatibility with call_writecounts
  colnames(annot) <- c("X1", "X2", "original_name")

  # NA FDR values should be set to 1 so all NA cells are filtered.
  edrops <- lapply(
    edrops,
    function(x) {
      x$FDR[is.na(x$FDR)] <- 1
      return(x)
    }
  )

  barcodes_to_keep <- lapply(
    edrops,
    function(x) {
      barcodes_to_keep <- rownames(x[x$FDR <= 0.001, ])
      return(barcodes_to_keep)
    }
  )

  #Populate filtered_samples list
  for(sample in sampleIds){
      message("Creating filtered sample matrices")
      counts <- counts_list[[sample]]
      counts_after_ed <- counts[, barcodes_to_keep[[sample]]]
      #save filtered samples in the filtered_samples list
      filtered_samples[[sample]] <- counts_after_ed
  }

  # Save plots for each sample
  for (sample in sampleIds) {
    counts <- counts_list[[sample]]
    ed_fdr <- edrops[[sample]]
    sample_barcodes_to_keep <- barcodes_to_keep[[sample]]

    bcrank <- DropletUtils::barcodeRanks(counts)
    filtered_bcrank <- bcrank[sample_barcodes_to_keep,]

    plot_data <- data.frame(rank = bcrank$rank, umis = bcrank$total, fdr = ed_fdr$FDR)

    #Make sure that we only plot desired fdr values
    filtered_ed_fdr <- ed_fdr[match(barcodes_to_keep[[sample]], rownames(ed_fdr)),]

    filtered_plot_data <- data.frame(rank = filtered_bcrank$rank, umis = filtered_bcrank$total, fdr = filtered_ed_fdr$FDR)

    g <- gridExtra::arrangeGrob(create_knee_plot(plot_data), create_knee_plot(filtered_plot_data), nrow = 2)
    ggsave(file = paste("qc_plots/emptyDrops/", sample, ".png"), g)
  }

  # Write sample files
  for (sample in sampleIds) {
    sample_counts <- filtered_samples[[sample]]

    if (all(rownames(sample_counts) %in% rownames(annot))) {
      sample_annot <- annot[rownames(sample_counts), ]
    } else {
      stop("Some features dont have annotations")
    }

    message("Writing files for sample ", sample)
    call_writecounts(sample, sample_counts, sample_annot, "./")
  }
}

#Requires df with columns rank, umis and fdr
create_knee_plot <- function(df){
    plot <- ggplot(unique_bcrank, aes(x = rank, y = umis, color = fdr)) +
        geom_point() +
        scale_y_continuous(trans = "log10") +
        scale_x_continuous(trans = "log10")
    return(plot)
}

filter_empty_drops()
