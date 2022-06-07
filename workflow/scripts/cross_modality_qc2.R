#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
library(optparse)
#library(config)
library(ggplot2)
library(gplots)
library(ggExtra)
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
    c("--dropkick"), type = "character", help = "Dropkick results to use"
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output metrics file"
  ),
  make_option(
    c("--sample"), type = "character", help = "Sample name"
  ),
  make_option(
    c("--min_reads"), type = "numeric", help = "Min reads for prefiltering (before outlier detection)."
  ),
  make_option(
    c("--barcode_whitelist_gex"), type = "character", help = "10X Multiome GEX barcode whitelist to use."
  ),
  make_option(
    c("--barcode_whitelist_atac"), type = "character", help = "10X Multiome ATAC barcode whitelist to use."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# ### Test
# opts=list()
# opts$sample = "Sample_5124-NM-5-hg38"
# opts$input_atac = paste0("results/multiome/nf_atac_results/ataqv/single-nucleus/",opts$sample,".txt")
# opts$input_gex = paste0("results/multiome/nf_gex_results/qc/",opts$sample,".qc.txt")
# opts$dropkick = paste0("results/multiome/dropkick_multiotsu_100/", opts$sample,"/dropkick.csv")
# opts$outdir = file.path("results/multiome/cross_modality_qc", opts$sample)
# opts$min_reads = 100
# opts$barcode_whitelist_gex = 'resources/barcode_whitelist_multiome_GEX.txt'
# opts$barcode_whitelist_atac = 'resources/barcode_whitelist_multiome_ATAC.txt'


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
rna_barcodes = scan(opts$barcode_whitelist_gex, what="character")
atac_barcodes = scan(opts$barcode_whitelist_atac, what="character")
rna_to_atac_barcodes = data.frame(rna_barcodes, atac_barcodes)
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
# Add dropkick scores
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dropkick = read.table(opts$dropkick, sep=",", header=TRUE)
d_all = merge.data.frame(d_joint, dropkick, by.x="barcode_gex", by.y="X", all.x=TRUE)

#rename dropkick vars to keep track of which values come from dropkick
dropkick_heuristics = c("n_genes_by_counts", "total_counts",
                        "pct_counts_mito","pct_counts_ambient",
                        "arcsinh_total_counts", "arcsinh_n_genes_by_counts")
names(d_all)[names(d_all) %in% dropkick_heuristics] <- paste0("dropkick.", names(d_all)[names(d_all) %in% dropkick_heuristics])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save joint qc metrics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metrics_all = c("umis", "hqaa", "frac_mt_gex", "frac_mt_atac",
            "tss_enrichment", "max_fraction_reads_from_single_autosome",
            "dropkick.n_genes_by_counts", "dropkick.total_counts",
            "dropkick.pct_counts_mito", "dropkick.pct_counts_ambient",
            "dropkick.arcsinh_total_counts","dropkick.arcsinh_n_genes_by_counts",
            "dropkick_score")

write.table(d_all[,c("droplet_id", "barcode_gex", "barcode_atac", metrics_all, "dropkick_label")], file=file.path(opts$outdir, "cross_modality_qc.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prefiltering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove barcodes with zero reads for either modality
d_all = d_all[d_all[,"hqaa"]>0 & d_all[,"umis"]>0,]

# Remove barcodes with missing values
d <- d_all[complete.cases(d_all),]

# Remove barcodes with very low read counts (less than min_reads threshold)
d = d[d$hqaa>=opts$min_reads & d$umis>=opts$min_reads,]

# Remove barcodes with 100% MT counts
d = d[d$frac_mt_gex<1 & d$frac_mt_atac<1,]

# Remove barcodes with 100% atac reads from one autosome
d = d[d$max_fraction_reads_from_single_autosome<1,]

cat("Barcode count after prefiltering:", nrow(d),"\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualize QC metric correlation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Which metric is dropkick score most correlated with?
metrics_quant = c("umis", "hqaa", "frac_mt_gex", "frac_mt_atac",
            "tss_enrichment", "max_fraction_reads_from_single_autosome",
            "dropkick.n_genes_by_counts", "dropkick.total_counts",
            "dropkick.pct_counts_mito", "dropkick.pct_counts_ambient",
            "dropkick_score")

metrics_corr = cor(d[,metrics_quant], method="spearman")

png(file.path(opts$outdir, "corr_metrics_heatmap.png"), width=800, height=800)
heatmap.2(metrics_corr, trace="none", margins=c(20,20), density.info="none", breaks=seq(-1, 1, by=0.1), col=bluered(length(seq(-1, 1, by=0.1))-1))
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save prefiltered droplets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.table(d[,c("droplet_id", "barcode_gex", "barcode_atac", metrics_all, "dropkick_label")], file=file.path(opts$outdir, "cross_modality_qc_prefiltered.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(d[,c("barcode_gex", "barcode_atac")], file=file.path(opts$outdir, "barcodes_prefiltered.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

### Formatted for AMULET
d_all$is_nucleus = as.numeric(d_all$droplet_id %in% d$droplet_id)
write.csv(d_all[,c("barcode_atac", "droplet_id", "is_nucleus")], file=file.path(opts$outdir, "barcodes_prefiltered_for_amulet.csv"), row.names=FALSE, quote=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save empty droplets list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hqaa_empty_thresh = 50
umis_empty_thresh = 10
barcodes_empty = d_joint %>%
  filter(hqaa < hqaa_empty_thresh) %>%
  filter(umis < umis_empty_thresh) %>%
  select(barcode_gex, barcode_atac) %>%
  write.table(file = file.path(opts$outdir, "barcodes_empty.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)


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



make2DDensityPlot_with_lines = function(d, x, y) {
  #remove zeros
  d = d[d[,x]>0 & d[,y]>0,]
  p = ggplot(d, aes_string(x=x, y=y)) +
    geom_bin2d(bins=100) +
    geom_point(alpha=0) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10") +
    #highlight empty droplets
    geom_segment(aes(x = hqaa_empty_thresh, y = 0, xend = hqaa_empty_thresh, yend = umis_empty_thresh), linetype = "longdash", size=0.1, colour="grey") +
    geom_segment(aes(x = 0, y = umis_empty_thresh, xend = hqaa_empty_thresh, yend = umis_empty_thresh), linetype = "longdash", size=0.1, colour="grey") +
    #highlight prefiltering thresholds
    geom_segment(aes(x = opts$min_reads, y = opts$min_reads, xend = max(d[,x]), yend = opts$min_reads), linetype = "longdash", size=0.1, colour="grey") +
    geom_segment(aes(x = opts$min_reads, y = opts$min_reads, xend = opts$min_reads, yend = max(d[,y])), linetype = "longdash", size=0.1, colour="grey") +
    scale_fill_distiller(palette=4, direction=1)
  p2 = ggMarginal(p, type="density")
  return(p2)
}



png(file.path(opts$outdir,"qc_reads_unfiltered.png"))
make2DDensityPlot_with_lines(d_all, "hqaa","umis")
dev.off()

png(file.path(opts$outdir,"qc_reads_prefiltered.png"))
make2DDensityPlot_with_lines(d, "hqaa","umis")
dev.off()



make2DDensityPlot = function(d, x, y) {
  p = ggplot(d, aes_string(x=x, y=y)) +
    geom_bin2d(bins=100) +
    geom_point(alpha=0) +
    scale_x_continuous(trans="log10") +
    scale_y_continuous(trans="log10") +
    geom_hline(yintercept=umis_empty_thresh, linetype = "longdash") +
    scale_fill_distiller(palette=4, direction=1)
  p2 = ggMarginal(p, type="density")
  return(p2)
}

addTitle = function(plots, titletext) {
  title <- ggdraw() + draw_label(titletext, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))
  p = plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins
  return(p)
}

makePlotGrid = function(d, plotfile) {
  qcplots = list()
  qcplots[[1]] = make2DDensityPlot(d, x = gex_metrics[["reads"]]$var, y = atac_metrics[["reads"]]$var)
  qcplots[[2]] = make2DDensityPlot(d, x = gex_metrics[["mt"]]$var, y = atac_metrics[["mt"]]$var)
  qcplots[[3]] = make2DDensityPlot(d, x = gex_metrics[["reads"]]$var, y = gex_metrics[["mt"]]$var)
  qcplots[[4]] = make2DDensityPlot(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["mt"]]$var)
  qcplots[[5]] = make2DDensityPlot(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["tss"]]$var)
  qcplots[[6]] = make2DDensityPlot(d, x = atac_metrics[["reads"]]$var, y = atac_metrics[["autosome"]]$var)
  pngfile=file.path(opts$outdir, plotfile)
  cat("Writing plot to", pngfile, "\n")
  png(pngfile, width=900, height=500)
  plotsonly = plot_grid(plotlist = qcplots, ncol=3)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$sample))
  print(plotsPlusTitle)
  dev.off()
}

makePlotGrid(d_all, "qc_grid_unfiltered.png")
makePlotGrid(d, "qc_grid_prefiltered.png")
