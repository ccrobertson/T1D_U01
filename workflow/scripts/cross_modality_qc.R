#!/usr/bin/env Rscript
options(stringsAsFactors=F)
library(optparse)
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))


### This script evaluates per-nucleus qc metrics from multiome GEX and ATAC jointly

# sample_name = "Sample_5124-NM-1-hg38"
# atacdir = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_atac_results/ataqv/single-nucleus"
# gexdir = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_gex_results/qc"

sample_name = "Sample_test-hg38"
atacdir = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_atac_test2_results/ataqv/single-nucleus"
gexdir = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_gex_test_results/qc"

d_atac = read.table(file.path(atacdir, paste0(sample_name,".txt")), header=FALSE)
d_gex = read.table(file.path(gexdir, paste0(sample_name,".qc.txt")), header=TRUE)

# Convert atac to short format
table(d_atac[1:100,"V2"])

summary(d_gex$barcode %in% d_atac$V1)








option_list <- list(
  make_option(
    c("-i", "--input"), type = "character", help = "ATAC QC file"
  ),
  make_option(
    c("-o", "--output"), type = "character", help = "Output file"
  ),
  make_option(
    c("--sample"), type = "character", help = "Sample name"
  ),
  make_option(
    c("--tss_min"), type = "numeric", help = "Minimum TSS enrichment"
  ),
  make_option(
    c("--tss_max"), type = "numeric", help = "Maximum TSS enrichment"
  ),
  make_option(
    c("--hqaa_min"), type = "numeric", help = "Minimum number of HQAA reads"
  ),
  make_option(
    c("--autosome_max"), type = "numeric", help = "Maximum fraction from single autosome"
  ),
  make_option(
    c("--pct_max"), type = "numeric",
    help = "Maximum percentile of HQAA reads after filtering for hqaa_min"
  ),
  make_option(
    c("--mito_knees"), type = "character",
    help = "TSV with sample and percent_mitochondrial cutoffs"
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# # Test
# opts$input <- "/lab/work/porchard/2020-nih-islet-sn/work/downstream/results/nucleus-qc/metrics.txt"
# opts$output <- "work/liger_rerun/HPAP045/barcodes_atac.txt"
# opts$sample <- "HPAP045"
# opts$pct_max <- 0.9
# opts$hqaa_min <- 5000
# opts$tss_min <- 3
# opts$tss_max <- 10


# Main parameters
infile <- opts$input
outfile <- opts$output
sample_to_keep <- opts$sample
hqaa_min <- opts$hqaa_min
pct_max <- opts$pct_max
tss_max <- opts$tss_max
tss_min <- opts$tss_min
autosome_max <- opts$autosome_max
mito_knee_file <- opts$mito_knees

# Load % MT cutoffs
if (!is.null(mito_knee_file)) {
  message("Filtering sample based on % MT")
  mt_knee <- read.table(mito_knee_file, header = T) %>%
    filter(sample == sample_to_keep) %>%
    pull(knee)
} else {
  mt_knee <- 100
}

# Read QC file and filter according to parameters
infile %>%
  read.table(col.names = c("barcode", "metric", "value")) %>%
  mutate(value = case_when(value == "None" ~ "0", T ~ value)) %>%  # deal gracefully with missing values
  mutate(barcode = gsub("-([Mm]ock|Cytokine|Cyto|CVB4)", "_\\1", barcode)) %>%  # deal with perturbations
  mutate(value = as.numeric(value)) %>%
  spread(metric, value) %>%
  mutate(barcode = gsub("-atac", "", barcode)) %>%
  separate(barcode, c("sample", "barcode"), "-") %>%
  mutate(sample = gsub("_", "-", sample)) %>%  # deal with perturbations
  filter(sample == sample_to_keep) %>%
  filter(hqaa > hqaa_min) %>%
  mutate(max_reads = unname(quantile(hqaa, pct_max))) %>%
  filter(hqaa < max_reads) %>%
  filter(tss_enrichment > tss_min) %>%
  filter(tss_enrichment < tss_max) %>%
  filter(max_fraction_reads_from_single_autosome < autosome_max) %>%
  filter(percent_mitochondrial < mt_knee) %>%
  pull(barcode) %>%
  write.table(file = outfile, col.names = F, row.names = F, quote = F)
