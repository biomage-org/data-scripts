#!/usr/bin/env python3
import pandas as pd
import loompy
import os
import sys
import scipy as sc
import gzip
import re
from anndata import AnnData

# This conversion script has been modified to specifically address column names and data types
# for Fly Cell Atlas (https://flycellatlas.org/)

def read_loom(file):

    # Validate because of this issue for loom files from Fly Cell Atlas:
    # https://github.com/scverse/anndata/issues/627
    with loompy.connect(file, validate=False) as ds:
        print(ds.shape)
        print(ds.ra.keys())
        print(ds.ca.keys())
        genes = ds.ra["Gene"]
        cells = ds.ca["CellID"]
        df = pd.DataFrame(ds[:, :], index=genes, columns=cells).T

    adata = AnnData(df)

    # FCA stores data as float32. Convert to integer so that it'll be stored as integers in the matrix file
    adata.X = adata.X.astype(int)

    return adata

def demultiplex(adata: AnnData, pattern: str):
    pattern = re.compile(pattern)
    row_indices = [i for i, name in enumerate(adata.obs_names) if pattern.match(name)]
    return adata[row_indices, :]


def write_counts(
    adata: AnnData,
    out: str,
    verbose: bool = False,
):
    os.makedirs(out, exist_ok=True)

    pd.DataFrame(adata.obs.index).to_csv(os.path.join(out, "barcodes.tsv"), sep="\t", index=False, header=None)
    pd.DataFrame(adata.var.index).to_csv(os.path.join(out, "features.tsv"), sep="\t", index=False, header=None)
    sc.io.mmwrite(os.path.join(out, "matrix.mtx"), sc.sparse.coo_matrix(adata.X.T), field="integer")

    if verbose:
        print("Wrote features.tsv, barcodes.tsv, and matrix.mtx files to {}".format(out))


def compress_and_delete(out):
    files = [os.path.join(out, f) for f in ['features.tsv', 'barcodes.tsv', 'matrix.mtx']]
    for file in files:
        with open(file, 'rb') as f_in, gzip.open(file + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)

        os.remove(file)

file = sys.argv[1]
project_name = os.path.splitext(os.path.basename(file))[0]

adata = read_loom(file)

os.makedirs(project_name, exist_ok=True)

# Insert the regex for sample_names here:
# demultiplex_meta = [{
#     "regex": "FCA_P31",
#     "sample_name": "FCA_P31"
# }, {...}]

demultiplex_meta = []

for sample_meta in demultiplex_meta:
    single_sample_data = demultiplex(adata, sample_meta["regex"])

    outpath = os.path.join(project_name, sample_meta["name"])
    write_counts(single_sample_data, outpath, True)
    compress_and_delete(outpath)
