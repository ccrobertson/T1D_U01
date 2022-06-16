#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
library(optparse)
library(ggplot2)
library(gplots)
library(ggExtra)
library(cowplot)
library(tidyr)
library(dplyr)

theme_set(theme_bw(base_size = 12))

# opts=list()
# opts$input_csv = "results/multiome/cross_modality_qc/Sample_5124-NM-1-hg38/cross_modality_qc_prefiltered_flagoutliers.csv"
# opts$outdir = "results/multiome/cross_modality_qc/Sample_5124-NM-1-hg38"
# opts$sample = "Sample_5124-NM-1-hg38"

option_list <- list(
  make_option(
    c("--input_csv"), type = "character", help = "Cross  modality QC metrics file"
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output directory for plots"
  ),
  make_option(
    c("--sample"), type = "character", help = "Sample name"
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


d = read.csv(opts$input_csv)
d$barcode_gex = d$X
d$IsolationForest_label = factor(d$IsolationForest_label, levels=c("True","False"), labels=c("Keep","Drop"))

dropkick_score_thresh = quantile(d$dropkick_score, probs=0.1)
d$dropkick_label = factor(d$dropkick_score > dropkick_score_thresh, levels=c(TRUE,FALSE), labels=c("Keep","Drop"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q5 Filtering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Testing what happens if I just remove any nucleus with a QC metric in the bottom
# 5th perentile of the distribution
percentile_cutoff = 0.05
thresholds = list()
thresholds[["umis"]] = quantile(d$umis, probs=percentile_cutoff)
thresholds[["hqaa"]] = quantile(d$hqaa, probs=percentile_cutoff)
thresholds[["frac_mt_gex"]] = quantile(d$frac_mt_gex, probs=1-percentile_cutoff)
thresholds[["frac_mt_atac"]] = quantile(d$frac_mt_atac, probs=1-percentile_cutoff)
thresholds[["max_fraction_reads_from_single_autosome"]] = quantile(d$max_fraction_reads_from_single_autosome, probs=1-percentile_cutoff)
thresholds[["tss_enrichment"]] = quantile(d$tss_enrichment, probs=percentile_cutoff)
thresholds[["dropkick_score"]] = quantile(d$dropkick_score, probs=percentile_cutoff)
d_filt = d %>%
  filter(umis > thresholds[["umis"]],
         hqaa > thresholds[["hqaa"]],
         frac_mt_gex < thresholds[["frac_mt_gex"]],
         frac_mt_atac < thresholds[["frac_mt_atac"]],
         max_fraction_reads_from_single_autosome < thresholds[["max_fraction_reads_from_single_autosome"]],
         tss_enrichment > thresholds[["tss_enrichment"]],
         dropkick_score > thresholds[["dropkick_score"]])
d$QuantileFilter_label = factor(d$droplet_id %in% d_filt$droplet_id, levels=c(TRUE, FALSE), labels=c("Keep","Drop"))


## Print thresholds for record:
print(thresholds)

# ### Explore some of the filtered out groups
# h = d[d$max_fraction_reads_from_single_autosome < 0.1 & d$hqaa>1e4 & d$QuantileFilter_label=="Drop",]
# mean(h$tss_enrichment < tss_enrichment_thresh)
# mean(h$dropkick_score < dropkick_score_thresh)
# mean(h$frac_mt_gex > frac_mt_gex_thresh)
# mean(h$frac_mt_atac > frac_mt_atac_thresh)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cross modality visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Visualizing qc metrics across modalities\n")
gex_metrics = list(
  reads = list(var="umis", label="UMIs in GEX"),
  mt = list(var="frac_mt_gex", label="Fraction MT in GEX")
)

atac_metrics = list(
  reads = list(var="hqaa", label="HQAA in ATAC"),
  mt = list(var="frac_mt_atac", label="Fraction MT in ATAC"),
  tss = list(var="tss_enrichment", label="TSS enrichment in ATAC"),
  autosome = list(var="max_fraction_reads_from_single_autosome", label="Max fraction HQAA from single autosome")
)

make2DDensityPlot = function(d, x, y) {
  p = ggplot(d, aes_string(x=x, y=y)) +
    geom_bin2d(bins=100) +
    geom_point(alpha=0) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10") +
    scale_fill_distiller(palette=4, direction=1)
  p2 = ggMarginal(p, type="density")
  return(p2)
}

addTitle = function(plots, titletext) {
  title <- ggdraw() + draw_label(titletext, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))
  p = plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins
  return(p)
}

makeDensityGrid = function(d, plotfile) {
  qcplots = list()
  qcplots[[1]] = make2DDensityPlot(d, x = gex_metrics[["reads"]]$var, y = atac_metrics[["reads"]]$var)
  qcplots[[2]] = make2DDensityPlot(d, x = gex_metrics[["mt"]]$var, y = atac_metrics[["mt"]]$var)
  qcplots[[3]] = make2DDensityPlot(d, x = gex_metrics[["reads"]]$var, y = gex_metrics[["mt"]]$var)
  qcplots[[4]] = make2DDensityPlot(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["mt"]]$var)
  qcplots[[5]] = make2DDensityPlot(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["tss"]]$var)
  qcplots[[6]] = make2DDensityPlot(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["autosome"]]$var)
  pngfile=plotfile
  cat("Writing plot to", pngfile, "\n")
  png(pngfile, width=900, height=700)
  plotsonly = plot_grid(plotlist = qcplots, ncol=2)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$sample))
  print(plotsPlusTitle)
  dev.off()
}

makeDensityGrid(d, file.path(opts$outdir, "qc_grid_prefiltered_density.png"))
makeDensityGrid(d[d$QuantileFilter_label=="Keep",], file.path(opts$outdir, "qc_grid_prefiltered_density_QuantileFilter_keep.png"))


make2DDensityPlot_by_factor = function(d, x, y, facet_by) {
  p = ggplot(d, aes_string(x=x, y=y)) +
    geom_bin2d(bins=100) +
    geom_point(alpha=0) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10") +
    scale_fill_distiller(palette=4, direction=1) +
    facet_grid(as.formula(paste(".~",facet_by)))
  return(p)
}

makeDensityGrid_by_factor = function(d, plotfile, facet_by) {
  qcplots = list()
  qcplots[[1]] = make2DDensityPlot_by_factor(d, x = gex_metrics[["reads"]]$var, y = atac_metrics[["reads"]]$var, facet_by=facet_by)
  qcplots[[2]] = make2DDensityPlot_by_factor(d, x = gex_metrics[["mt"]]$var, y = atac_metrics[["mt"]]$var, facet_by=facet_by)
  qcplots[[3]] = make2DDensityPlot_by_factor(d, x = gex_metrics[["reads"]]$var, y = gex_metrics[["mt"]]$var, facet_by=facet_by)
  qcplots[[4]] = make2DDensityPlot_by_factor(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["mt"]]$var, facet_by=facet_by)
  qcplots[[5]] = make2DDensityPlot_by_factor(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["tss"]]$var, facet_by=facet_by)
  qcplots[[6]] = make2DDensityPlot_by_factor(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["autosome"]]$var, facet_by=facet_by)
  pngfile=plotfile
  cat("Writing plot to", pngfile, "\n")
  png(pngfile, width=900, height=700)
  plotsonly = plot_grid(plotlist = qcplots, ncol=2)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$sample))
  print(plotsPlusTitle)
  dev.off()
}

makeDensityGrid_by_factor(d, file.path(opts$outdir, "qc_grid_prefiltered_density_by_dropkick_label.png"), facet_by="dropkick_label")
makeDensityGrid_by_factor(d, file.path(opts$outdir, "qc_grid_prefiltered_density_by_IsolationForest_label.png"), facet_by="IsolationForest_label")
makeDensityGrid_by_factor(d, file.path(opts$outdir, "qc_grid_prefiltered_density_by_QuantileFilter_label.png"), facet_by="QuantileFilter_label")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write final determination of nuclei to barcode file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.table(d[d$QuantileFilter_label=="Keep", c("barcode_gex", "barcode_atac")], file=file.path(opts$outdir, "barcodes_nuclei.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
