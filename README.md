# data-scripts
Scripts to demultiplex or pre-process user data to get into Cellenics.

# Utils.R
Contains utilities for the other functions.

# filter_empty_drops.R
Usage:
Download samples using aws s3 cp s3://biomage-originals-production/**PROJECTUUID** input --recursive
Copy samples table (Pure json as opposed to dynamodb JSON) into samples.json
run "python3 rename_samples.py"
Open data-scripts.rproj and load renv dependencies renv::restore()
Use the filter_empty_drops function to filter all samples in the input dir

