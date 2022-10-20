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
    c("--atac_qc"), type = "character", help = "ATAC qc metrics file"
  ),
  make_option(
    c("--barcode_map"), type = "character", help = "Multiome ATAC and GEX barcode map"
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

### Testing
# opts=list()
# opts$gex_qc = "results/multiome/nf_gex_results/qc/Sample_5124-NM-1-hg38.qc.txt"
# opts$atac_qc = "results/multiome/nf_atac_results/ataqv/single-nucleus/Sample_5124-NM-1-hg38.txt"
# opts$barcode_map = "resources/multiome_barcode_map.rds"
# opts$droplet_utils_nuclei = "results/multiome/droplet_utils/Sample_5124-NM-1/barcodes_nuclei.txt"
# opts$droplet_utils_empty = "results/multiome/droplet_utils/Sample_5124-NM-1/barcodes_empty.txt"
# opts$outdir = "results/multiome/cross_modality_qc/Sample_5124-NM-1"
# opts$sample = "Sample_5124-NM-1"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Reading GEX QC file", opts$gex_qc, "\n")
d_gex <- opts$gex_qc %>%
  read.table(header = TRUE) %>%
  select(barcode, umis, fraction_mitochondrial)

cat("Reading ATAC QC file:", opts$atac_qc, "\n")
d_atac <-  opts$atac_qc %>%
    read.table(col.names = c("barcode", "metric", "value")) %>%
    #mutate(value = case_when(value == "None" ~ "0", TRUE ~ value)) %>%  # deal gracefully with missing values
    pivot_wider(names_from = metric, values_from = value) %>%
    select(barcode, hqaa, tss_enrichment, percent_mitochondrial, max_fraction_reads_from_single_autosome)

droplet_utils_nuclei = scan(opts$droplet_utils_nuclei, what="character")
droplet_utils_empty = scan(opts$droplet_utils_empty, what="character")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get GEX --> ATAC barcode map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Getting GEX to ATAC barcode map:\n")
rna_to_atac_barcodes = readRDS(opts$barcode_map)
names(rna_to_atac_barcodes) <- c("rna_barcodes","atac_barcodes")
rna_to_atac_barcodes$droplet_id = paste0("droplet_", seq(1:nrow(rna_to_atac_barcodes)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine modalities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Merging QC data across modalities\n")
d_joint <- d_gex %>%
  inner_join(rna_to_atac_barcodes, by= c("barcode" = "rna_barcodes")) %>%
  rename(barcode_gex = barcode) %>%
  inner_join(d_atac, by = c("atac_barcodes" = "barcode")) %>%
  rename(barcode_atac = atac_barcodes) %>%
  mutate(frac_mt_atac = percent_mitochondrial/100) %>%
  rename(frac_mt_gex = fraction_mitochondrial) %>%
  select(droplet_id, barcode_gex, barcode_atac, umis, hqaa, frac_mt_gex, frac_mt_atac, tss_enrichment, max_fraction_reads_from_single_autosome)


cat("Total GEX barcode count:", nrow(d_gex),"\n")
cat("Total ATAC barcode count:", nrow(d_atac),"\n")
cat("Overlapping total barcode count :", nrow(d_joint),"\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Barcode filtering for nuclei
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Hard thresholds
hard_thresholds = list()
hard_thresholds[["umis"]] = 3000
hard_thresholds[["hqaa"]] = 5000
hard_thresholds[["frac_mt_gex"]] = 0.2
hard_thresholds[["frac_mt_atac"]] = 0.2
hard_thresholds[["max_fraction_reads_from_single_autosome"]] = 0.1
hard_thresholds[["tss_enrichment"]] = 2
hard_thresholds[["umis_empty"]] = 50
hard_thresholds[["hqaa_empty"]] = 100

### Filter for DropletUtils
d_droplet_utils = d_joint %>%
  filter(barcode_gex %in% droplet_utils_nuclei)

### Filter for hard thresholds
d_filt = d_droplet_utils %>%
  filter(umis > hard_thresholds[["umis"]],
         hqaa > hard_thresholds[["hqaa"]],
         frac_mt_gex < hard_thresholds[["frac_mt_gex"]],
         frac_mt_atac < hard_thresholds[["frac_mt_atac"]],
         max_fraction_reads_from_single_autosome < hard_thresholds[["max_fraction_reads_from_single_autosome"]],
         tss_enrichment > hard_thresholds[["tss_enrichment"]])


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
outlier_thresholds[["hqaa"]] = getOutlierThresholds(d_filt$hqaa)
outlier_thresholds[["frac_mt_gex"]] = getOutlierThresholds(d_filt$frac_mt_gex)
outlier_thresholds[["frac_mt_atac"]] = getOutlierThresholds(d_filt$frac_mt_atac)
outlier_thresholds[["max_fraction_reads_from_single_autosome"]] = getOutlierThresholds(d_filt$max_fraction_reads_from_single_autosome)
outlier_thresholds[["tss_enrichment"]] = getOutlierThresholds(d_filt$tss_enrichment)

### Filter for outliers
d_final = d_filt %>%
  filter(umis > outlier_thresholds[["umis"]][["lower"]],
         umis < outlier_thresholds[["umis"]][["upper"]],
         hqaa > outlier_thresholds[["hqaa"]][["lower"]],
         hqaa < outlier_thresholds[["hqaa"]][["upper"]],
         frac_mt_gex > outlier_thresholds[["frac_mt_gex"]][["lower"]],
         frac_mt_gex < outlier_thresholds[["frac_mt_gex"]][["upper"]],
         frac_mt_atac > outlier_thresholds[["frac_mt_atac"]][["lower"]],
         frac_mt_atac < outlier_thresholds[["frac_mt_atac"]][["upper"]],
         max_fraction_reads_from_single_autosome > outlier_thresholds[["max_fraction_reads_from_single_autosome"]][["lower"]],
         max_fraction_reads_from_single_autosome < outlier_thresholds[["max_fraction_reads_from_single_autosome"]][["upper"]],
         tss_enrichment > outlier_thresholds[["tss_enrichment"]][["lower"]],
         tss_enrichment < outlier_thresholds[["tss_enrichment"]][["upper"]])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Barcode filtering for empty droplets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d_droplet_utils_empty = d_joint %>%
  filter(barcode_gex %in% droplet_utils_empty)

d_empty = d_droplet_utils_empty %>%
  filter(umis < hard_thresholds[["umis_empty"]],
         hqaa < hard_thresholds[["hqaa_empty"]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save filtering results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Save filtering info for record:
filtering_summary = list(
  hard_thresholds = hard_thresholds,
  outlier_thresholds =outlier_thresholds,
  total_barcodes = nrow(d_joint),
  droplet_utils_nuclei = nrow(d_droplet_utils),
  droplet_utils_nuclei_passing_hard_threshold = nrow(d_filt),
  nuclei_after_outlier_filtering = nrow(d_final),
  droplet_utils_empty = nrow(d_droplet_utils_empty),
  droplet_utils_empty_passing_hard_thresholod = nrow(d_empty)
)
write_yaml(list(filtering_summary=filtering_summary), file=file.path(opts$outdir, "filtering_summary.yaml"))


### Save joint qc metrics
d_joint$nuclei = as.numeric(d_joint$droplet_id %in% d_final$droplet_id)
d_joint$empty =  as.numeric(d_joint$droplet_id %in% d_empty$droplet_id)
metrics = c("umis", "hqaa", "frac_mt_gex", "frac_mt_atac",
            "tss_enrichment", "max_fraction_reads_from_single_autosome", "nuclei", "empty")
write.table(d_joint[,c("droplet_id", "barcode_gex", "barcode_atac", metrics)], file=file.path(opts$outdir, "cross_modality_qc.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


### Save barcode lists
write.table(d_final$barcode_gex, file=file.path(opts$outdir, "barcodes_nuclei_gex.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(d_final$barcode_atac, file=file.path(opts$outdir, "barcodes_nuclei_atac.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(d_final[,c("droplet_id", "barcode_gex", "barcode_atac")], file=file.path(opts$outdir, "barcodes_nuclei.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

write.table(d_empty$barcode_gex, file=file.path(opts$outdir, "barcodes_empty_gex.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(d_empty$barcode_atac, file=file.path(opts$outdir, "barcodes_empty_atac.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(d_empty[,c("droplet_id", "barcode_gex", "barcode_atac")], file=file.path(opts$outdir, "barcodes_empty.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


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

makeDensityGrid(d_joint, file.path(opts$outdir, "qc_density_all.png"))
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

## Create facet variables
d_droplet_utils$nuclei = as.numeric(d_droplet_utils$droplet_id %in% d_final$droplet_id)
d_filt$nuclei = as.numeric(d_filt$droplet_id %in% d_final$droplet_id)

makeDensityGrid_by_factor(d_joint, plotfile=file.path(opts$outdir, "qc_density_all_facet_nuclei.png"), facet_by="nuclei")
makeDensityGrid_by_factor(d_joint, plotfile=file.path(opts$outdir, "qc_density_all_facet_empty.png"), facet_by="empty")
makeDensityGrid_by_factor(d_droplet_utils, plotfile=file.path(opts$outdir, "qc_density_droplet_utils_facet_nuclei.png"), facet_by="nuclei")
makeDensityGrid_by_factor(d_filt, plotfile=file.path(opts$outdir, "qc_density_hard_threshold_facet_nuclei.png"), facet_by="nuclei")
