#!/bin/python

import argparse
from os.path import join
import scanpy as sc
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--starsolodir', dest="starsolodir", default = None, required = True, help = 'Directory containing starsolo output to use to build anndata object.')
parser.add_argument('--barcodes', dest="barcodes", default = None, required = True, help = 'File containing a list of prefiltered barcodes to include.')
parser.add_argument('--qcstats', dest="qcstats", default = None, required = True, help = 'File containing cross modality qc metrics per barcode.')
parser.add_argument('--h5ad_file', dest="h5ad_file", default = None, required = True, help = 'Location to save anndata output.')
args = parser.parse_args()

### Testing
#ad = sc.read_10x_mtx("results/multiome/nf_gex_results/starsolo/Sample_5124-NM-1-hg38/Sample_5124-NM-1-hg38.Solo.out/GeneFull_ExonOverIntron/raw", make_unique=True)
#barcodes = pd.read_table("results/multiome/cross_modality_qc/Sample_5124-NM-1-hg38/barcodes_prefiltered.txt")
#qcstats = pd.read_table("results/multiome/cross_modality_qc/Sample_5124-NM-1-hg38/cross_modality_qc_prefiltered.txt")

## Read in AnnData object
print("Reading 10X GEX data.")
ad = sc.read_10x_mtx(args.starsolodir, make_unique=True)

## Read in barcodes and qcstats
print("Reading barcodes and qc stats.")
barcodes = pd.read_table(args.barcodes)
qcstats = pd.read_table(args.qcstats)
qcstats.set_index("barcode_gex", inplace=True)

## Filter anndata object for prefiltered barcodes
print("Filtering anndata object.")
ad_prefiltered = ad[barcodes["barcode_gex"],:]

## Add qc stats to anndata meta data
print("Adding qc stats to anndata object.")
ad_prefiltered.obs = pd.merge(ad_prefiltered.obs, qcstats, how="left", right_index=True, left_index=True)

## Write to file
print("Saving anndata object to h5ad file.")
ad_prefiltered.write(args.h5ad_file, compression="gzip")
#ad_prefiltered.write("results/multiome/cross_modality_qc/Sample_5124-NM-1-hg38/cross_modality_qc_prefiltered.h5ad", compression="gzip")
