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

# hto_demultiplex.R
Usage:
1. Follow usage instructions for filter_empty_drops.R to filter all samples in the input dir
2. Use the hto_demux function to demultiplex all samples in the out dir, which were previously filtered by the filter_empty_drops function
