#!/usr/bin/env Rscript
options(stringsAsFactors=F)
library(optparse)
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))


### This script evaluates per-nucleus qc metrics from multiome GEX and ATAC jointly

option_list <- list(
  make_option(
    c("-ia", "--input_atac"), type = "character", help = "ATAC QC file"
  ),
  make_option(
    c("-ig", "--input_gex"), type = "character", help = "GEX QC file"
  )
  make_option(
    c("-o", "--output"), type = "character", help = "Output file"
  ),
  make_option(
    c("--sample"), type = "character", help = "Sample name"
  ),
  make_option(
    c("--umi_min"), type = "numeric", help = "Minimum UMI counts"
  ),
  make_option(
    c("--pct_max"), type = "numeric",
    help = "Maximum percentile of UMI counts after filtering for umi_min"
  ),
  make_option(
    c("--mt_max"), type = "numeric", help = "Maximum MT fraction allowed"
  ),
  make_option(
    c("--pct_max_ribo"), type = "numeric",
    help = "Maximum percentile of % ribosomal after filtering for umi_min"
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


### Test
# opts$input_atac = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_atac_test2_results/ataqv/single-nucleus/Sample_test-hg38.txt"
# opts$input_gex = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_gex_test_results/qc/Sample_test-hg38.qc.txt"
# opts$output = "/lab/work/ccrober/T1D_U01/analysis/freeze1/joint_qc/Sample_test-hg38/barcodes_filtered.txt"
# opts$sample <- "Sample_test-hg38"
# opts$pct_max <- 0.9
# opts$hqaa_min <- 5000
# opts$tss_min <- 3
# opts$tss_max <- 10


# d_atac = read.table(opts$input_atac, header=FALSE)
# d_gex = read.table(opts$input_gex, header=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
infile_gex <- opts$input_gex
infile_atac <- opts$input_atac
outfile <- opts$output
sample <- opts$sample

#GEX filtering
umi_min <- opts$umi_min
pct_max <- opts$pct_max
mt_max <- opts$mt_max
pct_max_ribo <- opts$pct_max_ribo

#ATAC filtering
hqaa_min <- opts$hqaa_min
pct_max <- opts$pct_max
tss_max <- opts$tss_max
tss_min <- opts$tss_min
autosome_max <- opts$autosome_max
mito_knee_file <- opts$mito_knees

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Filtering GEX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- infile_gex %>%
  read.table(header = T, fill = T) %>%
  filter(number_umis > umi_min) %>%
  mutate(max_umi = unname(quantile(number_umis, pct_max)))

if(!is.null(pct_max_ribo)) {
  message("Also filtering by % ribosomal UMIs")
  ribo_qc <- sprintf("data/rna/ribosomal_qc/%s.txt", sample_to_keep) %>%
    read.table(header = T) %>%
    select(-number_umis)

  d <- d %>%
    left_join(ribo_qc, by = "barcode") %>%
    mutate(max_ribo = unname(quantile(fraction_ribosomal, pct_max_ribo))) %>%
    filter(fraction_ribosomal < max_ribo)
}

d %>%
  filter(number_umis < max_umi) %>%
  filter(fraction_mitochondrial < mt_max) %>%
  pull(barcode) %>%
  write.table(file = outfile, col.names = F, row.names = F, quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Filtering ATAC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
infile_atac %>%
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Cross modality filtering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
