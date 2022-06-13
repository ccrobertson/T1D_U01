import argparse
from os.path import join
import scanpy as sc; sc.set_figure_params(color_map="viridis", frameon=False)
import dropkick as dk
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from skimage.filters import (
    threshold_li,
    threshold_otsu,
    threshold_mean,
    threshold_multiotsu,
)
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--indir', dest="indir", default = None, required = True, help = 'Directory containing starsolo output to use for dropkick.')
parser.add_argument('--outdir', dest="outdir", default = None, required = True, help = 'Location to save dropkick output.')
parser.add_argument('--cpus', dest="cpus", type=int, default = None, required = True, help = 'Number of threads for fitting glmnet model.')
parser.add_argument('--min_genes', dest="min_genes", type=int, default = None, required = True, help = 'Exclude barcodes with fewer than this many genes.')
parser.add_argument('--method', dest="method", default = None, required = True, help = 'Which thresholding method to use.')
args = parser.parse_args()

outdir = args.outdir

## Read in AnnData object
print("Reading 10X GEX data.")
ad = sc.read_10x_mtx(args.indir, make_unique=True)


## Preprocessing
print("Preprocessing.")
ad = dk.recipe_dropkick(ad, n_hvgs=None, X_final="raw_counts", min_genes=args.min_genes, mito_names='^mt-|^MT-|^Mt-')

## QC summary
print("Running dropkick qc.")
qc_plt = dk.qc_summary(ad)
qc_plt.savefig(join(args.outdir,"dropkick_qc.png"))
qc_plt.clf()

## Thresholding
print("Plotting thresholding.")
# generate threshold dicts
vals = np.array(ad.obs["arcsinh_n_genes_by_counts"])
thresh_multiotsu = {'arcsinh_n_genes_by_counts': {'direction': 'above', 'thresh': threshold_multiotsu(vals)}}
thresh_otsu = {'arcsinh_n_genes_by_counts': {'direction': 'above', 'thresh': threshold_otsu(vals)}}
thresh_li = {'arcsinh_n_genes_by_counts': {'direction': 'above', 'thresh': threshold_li(vals)}}
thresh_mean = {'arcsinh_n_genes_by_counts': {'direction': 'above', 'thresh': threshold_mean(vals)}}

# plot thresholds relative to arcsinh_n_genes_by_counts distribution
thresh_plot1 = dk.plot_thresh_obs(ad, thresholds = thresh_multiotsu, save_to=join(outdir,"dropkick_thresholding_multiotsu.png"))
thresh_plot2 = dk.plot_thresh_obs(ad, thresholds = thresh_otsu, save_to=join(outdir,"dropkick_thresholding_otsu.png"))
thresh_plot3 = dk.plot_thresh_obs(ad, thresholds = thresh_li, save_to=join(outdir,"dropkick_thresholding_li.png"))
thresh_plot4 = dk.plot_thresh_obs(ad, thresholds = thresh_mean, save_to=join(outdir,"dropkick_thresholding_mean.png"))


## Fit logistic model
print("Fitting dropkick model.")
ad_model = dk.dropkick(ad, n_jobs=args.cpus, mito_names='^mt-|^MT-|^Mt-', thresh_methods=[args.method], min_genes=args.min_genes)

thresh_auto = ad.uns["dropkick_thresholds"] #thresholds produced by running dropkick
thresh_plot0 = dk.plot_thresh_obs(ad, thresholds = thresh_auto, save_to=join(outdir,"dropkick_thresholding_auto.png"))


## Save results to file
print("Saving model output to h5ad and csv files.")
print(ad)
ad.write(join(args.outdir, "dropkick.h5ad"), compression="gzip")
ad.obs.to_csv(join(args.outdir, "dropkick.csv"), index=True, header=True)

print("Done.")
