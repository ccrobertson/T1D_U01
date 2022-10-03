#!/bin/python

import argparse
from os.path import join
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest


parser = argparse.ArgumentParser()
parser.add_argument('--input_h5ad', dest="input_h5ad", default = None, required = True, help = 'Input file.')
parser.add_argument('--metadata_columns', dest="metadata_columns", default = None, required = True, help = 'Per barcode QC metrics to use for detecting outliers.')
parser.add_argument('--contamination', dest="contamination", type = float, default = None, required = True, help = 'Percent contamination to remove (should be a fraction or the string "auto").')
parser.add_argument('--output_h5ad', dest="output_h5ad", default = None, required = True, help = 'Output file.')
parser.add_argument('--output_csv', dest="output_csv", default = None, required = True, help = 'Output csv file.')
args = parser.parse_args()

### Testing
#ad = sc.read("results/multiome/cross_modality_qc/Sample_5124-NM-1-hg38/cross_modality_qc_prefiltered.h5ad")
#metadata_columns = ["umis", "hqaa", "frac_mt_gex", "frac_mt_atac", "tss_enrichment", "max_fraction_reads_from_single_autosome", "dropkick_score"]
#metadata_columns = 'hqaa,tss_enrichment'

## Reading input
ad = sc.read(args.input_h5ad)
metadata_columns = args.metadata_columns.split(',')

## Define random forest parameters
clf = IsolationForest(
    #n_estimators=500,
    max_samples=0.1,
    warm_start=False,
    contamination=args.contamination,
    random_state=0,
    bootstrap=True
)

## Run IsolationForest on meta data
print("Running random forest.")
clf = clf.fit(ad.obs[metadata_columns].values)

## Add results back to anndata object
ad.obs["IsolationForest_label"] = clf.predict(ad.obs[metadata_columns].values) == 1

## Save results
print("Saving IsolationForest output to h5ad and csv files.")
ad.write(args.output_h5ad, compression="gzip")
ad.obs.to_csv(args.output_csv, index=True, header=True)

print("Done.")
