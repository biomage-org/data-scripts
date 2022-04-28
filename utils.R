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

