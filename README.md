# data-scripts
Scripts to demultiplex or pre-process user data to get into Cellenics.

# Utils.R
Contains utilities for the other functions.

# filter_empty_drops.R
Usage:
1. Download samples using aws s3 cp s3://biomage-originals-production/**PROJECTUUID** input --recursive
2. Copy samples table (Pure json as opposed to dynamodb JSON) into samples.json
3. run "python3 rename_samples.py"
4. Open data-scripts.rproj and load renv dependencies renv::restore()
5. Use the filter_empty_drops function to filter all samples in the input dir

# cellset extraction

To extract cellsets, you only need an experiment ID, and the index of the cellset
in the cellsets file. The cellset index is composed from the cellset class, as listed
below, and the cellset number inside each class (with 1-based indexing).

- 1 = louvain
- 2 = scratchpad
- 3 = samples
- 4 = metadata tracks

The easiest way to get this is to use Rstudio's list viewer: `View(parsed_json_object)`.
You can download the cellset file using the `download_cellset_file` function, import it
with `jsonlite::read_json` and explore it with the list viewer.

After getting the positions, The function `extract_cellset` will do everything
automagically, returning a subsetted seurat object.
