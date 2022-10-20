#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
library(optparse)
library(ggplot2)
library(gplots)
library(ggExtra)
library(cowplot)
library(tidyr)
library(dplyr)
library(yaml)

theme_set(theme_bw(base_size = 12))


option_list <- list(
  make_option(
    c("--gex_qc"), type = "character", help = "GEX qc metrics file"
  ),
  make_option(
    c("--droplet_utils_nuclei"), type = "character", help = "DropletUtils nuclei barcodes list"
  ),
  make_option(
    c("--droplet_utils_empty"), type = "character", help = "DropletUtils empty droplets barcodes list"
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

# ## Testing
# opts=list()
# opts$gex_qc = "results/fiveprime/qc/Sample_5125-NM-2.qc.txt"
# opts$droplet_utils_nuclei = "results/fiveprime/droplet_utils/Sample_5125-NM-2/barcodes_nuclei.txt"
# opts$droplet_utils_empty = "results/fiveprime/droplet_utils/Sample_5125-NM-2/barcodes_empty.txt"
# opts$outdir = "results/fiveprime/barcode_filtering/Sample_5125-NM-2"
# opts$sample = "Sample_5125-NM-2"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Reading GEX QC file", opts$gex_qc, "\n")
d_gex <- opts$gex_qc %>%
  read.table(header = TRUE) %>%
  select(barcode, umis, fraction_mitochondrial) %>%
  rename(frac_mt_gex = fraction_mitochondrial)

droplet_utils_nuclei = scan(opts$droplet_utils_nuclei, what="character")
droplet_utils_empty = scan(opts$droplet_utils_empty, what="character")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter barcodes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hard_thresholds = list()
hard_thresholds[["umis"]] = 3000
hard_thresholds[["frac_mt_gex"]] = 0.2
hard_thresholds[["umis_empty"]] = 50

### Filter for DropletUtils
d_droplet_utils = d_gex %>%
  filter(barcode %in% droplet_utils_nuclei)

### Filter for hard thresholds
d_filt = d_droplet_utils %>%
  filter(umis > hard_thresholds[["umis"]],
         frac_mt_gex < hard_thresholds[["frac_mt_gex"]])


### Emperical outlier thresholds
getOutlierThresholds = function(x) {
  M = median(x)
  SD = sd(x)
  lower = M - 2*SD
  upper = M + 2*SD
  return(c(lower=lower, upper=upper))
}
outlier_thresholds = list()
outlier_thresholds[["umis"]] = getOutlierThresholds(d_filt$umis)
outlier_thresholds[["frac_mt_gex"]] = getOutlierThresholds(d_filt$frac_mt_gex)


### Filter for outliers
d_final = d_filt %>%
  filter(umis > outlier_thresholds[["umis"]][["lower"]],
         umis < outlier_thresholds[["umis"]][["upper"]],
         frac_mt_gex > outlier_thresholds[["frac_mt_gex"]][["lower"]],
         frac_mt_gex < outlier_thresholds[["frac_mt_gex"]][["upper"]])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Barcode filtering for empty droplets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d_droplet_utils_empty = d_gex %>%
  filter(barcode %in% droplet_utils_empty)

d_empty = d_droplet_utils_empty %>%
  filter(umis < hard_thresholds[["umis_empty"]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save filtering results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Save filtering info for record:
filtering_summary = list(
  hard_thresholds = hard_thresholds,
  outlier_thresholds =outlier_thresholds,
  total_barcodes = nrow(d_gex),
  droplet_utils_nuclei = nrow(d_droplet_utils),
  droplet_utils_nuclei_passing_hard_threshold = nrow(d_filt),
  nuclei_after_outlier_filtering = nrow(d_final),
  droplet_utils_empty = nrow(d_droplet_utils_empty),
  droplet_utils_empty_passing_hard_thresholod = nrow(d_empty)
)
write_yaml(list(filtering_summary=filtering_summary), file=file.path(opts$outdir, "filtering_summary.yaml"))


### Save joint qc metrics
d_gex$nuclei = as.numeric(d_gex$barcode %in% d_final$barcode)
d_gex$empty =  as.numeric(d_gex$barcode %in% d_empty$barcode)
metrics = c("umis", "frac_mt_gex", "nuclei", "empty")
write.table(d_gex[,c("barcode", metrics)], file=file.path(opts$outdir, "barcode_filtering_qc.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

### Save barcode lists
write.table(d_final$barcode, file=file.path(opts$outdir, "barcodes_nuclei.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(d_empty$barcode, file=file.path(opts$outdir, "barcodes_empty.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cross modality visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Visualizing qc metrics across modalities\n")
gex_metrics = list(
  reads = list(var="umis", label="UMIs in GEX"),
  mt = list(var="frac_mt_gex", label="Fraction MT in GEX")
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
  qcplots = list(make2DDensityPlot(d, x = gex_metrics[["reads"]]$var, y = gex_metrics[["mt"]]$var))
  pngfile=plotfile
  cat("Writing plot to", pngfile, "\n")
  png(pngfile, width=900, height=700)
  plotsonly = plot_grid(plotlist = qcplots, ncol=1)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$sample))
  print(plotsPlusTitle)
  dev.off()
}

makeDensityGrid(d_gex, file.path(opts$outdir, "qc_density_all.png"))
makeDensityGrid(d_droplet_utils, file.path(opts$outdir, "qc_density_droplet_utils.png"))
makeDensityGrid(d_filt, file.path(opts$outdir, "qc_density_hard_threshold.png"))
makeDensityGrid(d_final, file.path(opts$outdir, "qc_density_outlier_filtering.png"))



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
  qcplots[[1]] = make2DDensityPlot_by_factor(d, x = gex_metrics[["reads"]]$var, y = gex_metrics[["mt"]]$var, facet_by=facet_by)
  pngfile=plotfile
  cat("Writing plot to", pngfile, "\n")
  png(pngfile, width=900, height=700)
  plotsonly = plot_grid(plotlist = qcplots, ncol=1)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$sample))
  print(plotsPlusTitle)
  dev.off()
}

## Create facet variables
d_droplet_utils$nuclei = as.numeric(d_droplet_utils$barcode %in% d_final$barcode)
d_filt$nuclei = as.numeric(d_filt$barcode %in% d_final$barcode)

makeDensityGrid_by_factor(d_gex, plotfile=file.path(opts$outdir, "qc_density_all_facet_nuclei.png"), facet_by="nuclei")
makeDensityGrid_by_factor(d_gex, plotfile=file.path(opts$outdir, "qc_density_all_facet_empty.png"), facet_by="empty")
makeDensityGrid_by_factor(d_droplet_utils, plotfile=file.path(opts$outdir, "qc_density_droplet_utils_facet_nuclei.png"), facet_by="nuclei")
makeDensityGrid_by_factor(d_filt, plotfile=file.path(opts$outdir, "qc_density_hard_threshold_facet_nuclei.png"), facet_by="nuclei")
