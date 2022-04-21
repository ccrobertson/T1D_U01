#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
library(optparse)
library(config)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)

theme_set(theme_bw(base_size = 12))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script evaluates per-nucleus qc metrics from multiome GEX and ATAC
# and filters nuclei based on joint-modality thresholds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

option_list <- list(
  make_option(
    c("--input_gex"), type = "character", help = "GEX QC file"
  ),
  make_option(
    c("--input_atac"), type = "character", help = "ATAC QC file"
  ),
  make_option(
    c("--input_dropkick"), type = "character", help = "Dropkick csv file"
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output metrics file"
  ),
  make_option(
    c("--output_amulet"), type = "character", help = "Output csv file to use for running AMULET."
  ),
  make_option(
    c("--sample"), type = "character", help = "Sample name"
  ),
  make_option(
    c("--barcode_filtering_config"), type = "character", help = "Config file specifying pre-defined filtering thresholds for barcode qc metrics."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# # ### Test
# opts = list()
# opts$input_atac = "results/freeze1/nf_atac_results/ataqv/single-nucleus/Sample_test-hg38.txt"
# opts$input_gex = "results/freeze1/nf_gex_results/qc/Sample_test-hg38.qc.txt"
# opts$input_dropkick = "results/freeze1/dropkick/Sample_test-hg38/dropkick.csv"
# opts$outdir = "results/freeze1/cross_modality_qc/Sample_test-hg38"
# opts$output_amulet = "results/freeze1/amulet/Sample_test-hg38/amulet_input.csv"
# opts$sample <- "Sample_test-hg38"
# opts$barcode_filtering_config = "workflow/scripts/cross_modality_qc.yaml"

config <-  config::get(file = opts$barcode_filtering_config)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Reading GEX QC file", opts$input_gex, "\n")
d_gex <- opts$input_gex %>%
  read.table(header = TRUE) %>%
  select(barcode, umis, fraction_mitochondrial)

cat("Reading ATAC QC file:", opts$input_atac, "\n")
d_atac <-  opts$input_atac %>%
    read.table(col.names = c("barcode", "metric", "value")) %>%
    #mutate(value = case_when(value == "None" ~ "0", TRUE ~ value)) %>%  # deal gracefully with missing values
    pivot_wider(names_from = metric, values_from = value) %>%
    select(barcode, hqaa, tss_enrichment, percent_mitochondrial, max_fraction_reads_from_single_autosome)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get GEX --> ATAC barcode map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Peter does this here: /lab/work/porchard/islet-multiome/bin/rna-to-atac-nucleus-info.py
# def convert(nucleus, rna_to_atac):
#     library, genome, barcode = nucleus.split('-')
#     return f'{library}-{genome}-{rna_to_atac[barcode]}'
#
# # map RNA --> ATAC barcodes
# atac_barcodes = pd.read_csv(ATAC_BARCODE_LIST, header=None, names=['atac_barcode'])
# rna_barcodes = pd.read_csv(RNA_BARCODE_LIST, header=None, names=['rna_barcode'])
# barcodes = pd.concat([atac_barcodes, rna_barcodes], axis=1)
# barcode_rna_to_atac = dict(zip(barcodes.rna_barcode, barcodes.atac_barcode))

cat("Getting GEX to ATAC barcode map:\n")
RNA_BARCODE_LIST = 'resources/barcode_whitelist_multiome_GEX.txt'
ATAC_BARCODE_LIST = 'resources/barcode_whitelist_multiome_ATAC.txt'
rna_barcodes = scan(RNA_BARCODE_LIST, what="character")
atac_barcodes = scan(ATAC_BARCODE_LIST, what="character")
rna_to_atac_barcodes = data.frame(rna_barcodes, atac_barcodes)
rna_to_atac_barcodes$droplet_id = paste0("droplet_", seq(1:nrow(rna_to_atac_barcodes)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine modalities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Merging QC data across modalities\n")
#d_gex$barcode = sample(d_atac$barcode, size=nrow(d_gex), replace=FALSE) ### For testing plot code
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
cat("Overlapping barcode count:", nrow(d_joint),"\n")

#save joint metrics
write.table(d_joint, file=file.path(opts$outdir, "cross_modality_qc.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define barcode groups
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Sanity checks
sanityCheck = function(d, group_var) {
  cat("########### Summary by", group_var, "###########\n")
  cat("Barcodes total:", nrow(d),"\n")
  print(d %>%
   group_by(base::get(group_var)) %>%
   summarize(n = n(),
            umis_mean = mean(umis, na.rm=TRUE),
            umis_med = median(umis, na.rm=TRUE),
            umis_max = max(umis, na.rm=TRUE),
            hqaa_mean = mean(hqaa, na.rm=TRUE),
            hqaa_median = median(hqaa, na.rm=TRUE),
            hqaa_max = max(hqaa, na.rm=TRUE)))
  cat("\n\n")
}


### Deal with missing values
### Sanity checks -- missing TSS enrichment score
d_joint$tss_enrichment_missing = is.na(d_joint$tss_enrichment)
sanityCheck(d_joint, "tss_enrichment_missing")
d_joint <- d_joint[complete.cases(d_joint),]


################ Pre-defined filter thresholds

cat("Filtering barcodes\n")
#config <-  config::get(file = opts$barcode_filtering_config)

### Filtering for barcodes passing all filters
d_filter_pass = d_joint %>%
  #number of reads
  filter(umis > config$umi_min) %>%
  filter(hqaa > config$hqaa_min) %>%
  mutate(max_umi = unname(quantile(umis, config$pct_max_gex))) %>%
  mutate(max_hqaa = unname(quantile(hqaa, config$pct_max_atac))) %>%
  filter(umis < max_umi) %>%
  filter(hqaa < max_hqaa) %>%
  #%MT reads
  filter(frac_mt_gex < config$mt_max_gex) %>%
  filter(frac_mt_atac < config$mt_max_atac) %>%
  #ATAC filters
  filter(tss_enrichment > config$tss_min) %>%
  filter(tss_enrichment < config$tss_max) %>%
  filter(max_fraction_reads_from_single_autosome < config$autosome_max)

cat("Filter-pass barcode count:", nrow(d_filter_pass),"\n")

### Filtering for barcodes failing filters
d_filter_fail <- d_joint %>%
  filter(umis < config$umi_min) %>%
  filter(hqaa < config$hqaa_min) %>%
  #filter(frac_mt_gex > config$mt_max_gex | frac_mt_gex > config$mt_max_gex) %>%
  #filter(frac_mt_atac > config$mt_max_atac) %>%
  filter(tss_enrichment < config$tss_min)

cat("Filter-fail barcode count:", nrow(d_filter_fail),"\n")

d_joint$status_filters = "not_classified"
d_joint$status_filters[d_joint$barcode_gex %in% d_filter_pass$barcode_gex] <- "nuclei"
d_joint$status_filters[d_joint$barcode_gex %in% d_filter_fail$barcode_gex] <- "empty"

sanityCheck(d_joint, "status_filters")


################ K-means clustering on metrics

### Transform metrics to cluster on
transformMetric = function(x) {
  return(scale(log10(x + quantile(x[x>0],0.001, na.rm=TRUE) )))
}
metrics = d_joint[,c("umis","hqaa","tss_enrichment","max_fraction_reads_from_single_autosome", "frac_mt_gex", "frac_mt_atac")]
transformed_metrics = apply(metrics, MARGIN=2, FUN=transformMetric)

### Cluster by read numbers
d_joint$cluster2var = factor(kmeans(transformed_metrics[,c("umis","hqaa")], centers = 3)$cluster)

### Cluster by read numbers and TSS enrichment
d_joint$cluster3var = factor(kmeans(transformed_metrics[,c("umis","hqaa","tss_enrichment")], centers = 3)$cluster)

### Cluster by everything
d_joint$cluster6var = factor(kmeans(transformed_metrics[,c("umis","hqaa","tss_enrichment","max_fraction_reads_from_single_autosome", "frac_mt_gex", "frac_mt_atac")], centers = 3)$cluster)

sanityCheck(d_joint, "cluster2var")
sanityCheck(d_joint, "cluster3var")
sanityCheck(d_joint, "cluster6var")



############### Dropkick

### Get dropkick scores
cat("Getting dropkick scores\n")
dropkick = read.csv(opts$input_dropkick)
d_joint = d_joint %>%
  full_join(dropkick, by = c("barcode_gex" = "X"))

d_joint$status_dropkick = NA
d_joint$status_dropkick[is.na(d_joint$dropkick_label)] <- "not_scored"
d_joint$status_dropkick[d_joint$dropkick_label=="True"] <- "nuclei"
d_joint$status_dropkick[d_joint$dropkick_label=="False"] <- "empty"

### Some sanity checks -- dropkick score
sanityCheck(d_joint, "status_dropkick")
### How to see what dropkick score cutoff is?
# ad.uns["dropkick_thresholds"]["arcsinh_n_genes_by_counts"]["thresh"]
png(file.path(opts$outdir, "cross_modality_qc_plot_dropkick_score_dist.png"), width=500, height=800)
p = list()
p[[1]] = ggplot(dropkick) + geom_point(aes(x=arcsinh_n_genes_by_counts, y=pct_counts_ambient, colour=dropkick_score))
p[[2]] = ggplot(dropkick) + geom_point(aes(x=arcsinh_n_genes_by_counts, y=pct_counts_ambient, colour=dropkick_label))
plot_grid(plotlist = p, ncol=1, nrow=2)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cross modality visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Visualizing qc metrics across modalities\n")
gex_metrics = list(
  reads = list(var="umis", label="UMIs in GEX", thresh=config$umi_min),
  mt = list(var="frac_mt_gex", label="Fraction MT in GEX", thresh=config$mt_max_gex)
)

atac_metrics = list(
  reads = list(var="hqaa", label="HQAA in ATAC", thresh=config$hqaa_min),
  mt = list(var="frac_mt_atac", label="Fraction MT in ATAC", thresh=config$mt_max_atac),
  tss = list(var="tss_enrichment", label="TSS enrichment in ATAC", thresh=config$tss_min),
  autosome = list(var="max_fraction_reads_from_single_autosome", label="Max fraction HQAA from single autosome", thresh=config$autosome_max)
)

makeQCPlot = function(x, y, colour_by) {
  p = ggplot(d_joint) +
    geom_point(aes_string(x = x$var, y = y$var, colour = colour_by), shape=16, size=1, alpha=0.5) +
    xlab(x$label) + ylab(y$label) +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans="log10") +
    geom_vline(xintercept=x$thresh, linetype = "longdash") +
    geom_hline(yintercept=y$thresh, linetype = "longdash") +
    theme(legend.position="none")
  return(p)
}

addTitle = function(plots, titletext) {
  title <- ggdraw() + draw_label(titletext, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))
  p = plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins
  return(p)
}

makePlotGrid = function(colour_by) {
  qcplots = list()
  qcplots[[1]] = makeQCPlot(x = gex_metrics[["reads"]], y = atac_metrics[["reads"]], colour_by = colour_by)
  qcplots[[2]] = makeQCPlot(x = gex_metrics[["mt"]], y = atac_metrics[["mt"]], colour_by = colour_by)
  qcplots[[3]] = makeQCPlot(x = gex_metrics[["reads"]], y = gex_metrics[["mt"]], colour_by = colour_by)
  qcplots[[4]] = makeQCPlot(x = atac_metrics[["reads"]], y = atac_metrics[["mt"]], colour_by = colour_by)
  qcplots[[5]] = makeQCPlot(x = atac_metrics[["reads"]], y = atac_metrics[["tss"]], colour_by = colour_by)
  qcplots[[6]] = makeQCPlot(x = atac_metrics[["reads"]], y = atac_metrics[["autosome"]], colour_by = colour_by)
  pngfile=file.path(opts$outdir, paste0("cross_modality_qc_plot_",colour_by,".png"))
  cat("Writing plot to", pngfile, "\n")
  png(pngfile, width=900, height=500)
  plotsonly = plot_grid(plotlist = qcplots, ncol=3)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$sample, "-",colour_by))
  print(plotsPlusTitle)
  dev.off()
}

makePlotGrid(colour_by = "status_filters")
makePlotGrid(colour_by = "status_dropkick")
makePlotGrid(colour_by = "cluster2var")
makePlotGrid(colour_by = "cluster3var")
makePlotGrid(colour_by = "cluster6var")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save barcode lists
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Writing empty droplet barcodes to:", file.path(opts$outdir, "barcodes_empty.txt"),"\n")
d_joint %>%
  filter(status_dropkick == "empty") %>%
  select(barcode_gex, barcode_atac) %>%
  write.table(file = file.path(opts$outdir, "barcodes_empty.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)


cat("Writing quality nuclei barcodes to:", file.path(opts$outdir, "barcodes_nuclei.txt"),"\n")
d_joint %>%
  filter(status_dropkick == "nuclei") %>%
  select(barcode_gex, barcode_atac) %>%
  write.table(file = file.path(opts$outdir, "barcodes_nuclei.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create output for AMULET
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Which approach will we use?
dropkick_nuclei = d_joint$barcode_atac[d_joint$status_dropkick=="nuclei"]
amulet_csv = d_atac %>%
    left_join(rna_to_atac_barcodes, by= c("barcode" = "atac_barcodes")) %>%
    mutate(is_nucleus = as.numeric(barcode %in% dropkick_nuclei)) %>%
    select(barcode, droplet_id, is_nucleus)

write.csv(amulet_csv, file=opts$output_amulet, row.names=FALSE, quote=FALSE)
