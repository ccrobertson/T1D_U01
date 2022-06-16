options(stringsAsFactors=FALSE)
library(optparse)
library(celda)
library(rliger)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in 10x GEX data (starsolo output)
# and generates dgCMatrix objects for nuclei and empty droplets
# as input for decontx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ctrl_dge <- readRDS("../liger_tutorial_data/PBMC_control.RDS");
#stim_dge <- readRDS("../liger_tutorial_data/PBMC_interferon-stimulated.RDS");

option_list <- list(
  make_option(
    c("--input_10x_dir"), type = "character", help = "Directory of starsolo output containing matrix.mtx and genes.tsv files"
  ),
  make_option(
    c("--barcodes_nuclei"), type = "character", help = "File containing multiome barcodes from quality nuclei."
  ),
  make_option(
    c("--barcodes_empty"), type = "character", help = "File containing multiome barcodes from empty droplets."
  ),
  make_option(
    c("--counts_nuclei"), type = "character", help = "Output file name for storing quality nuclei counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--counts_empty"), type = "character", help = "Output file name for storing empty droplets counts matrix (dgCMatrix format)."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# ### Testing
# opts = list()
# opts$input_10x_dir = "results/freeze1/nf_gex_results/starsolo/Sample_test-hg38/Sample_test-hg38.Solo.out/GeneFull_ExonOverIntron/raw"
# opts$barcodes_nuclei = "results/freeze1/cross_modality_qc/Sample_test-hg38/barcodes_nuclei.txt"
# opts$barcodes_empty = "results/freeze1/cross_modality_qc/Sample_test-hg38/barcodes_empty.txt"
# opts$barcodes_doublets = "results/freeze1/amulet/Sample_5124-NM-1-hg38/MultipletBarcodes_01.txt"
# opts$counts_nuclei = "results/freeze1/decontx/Sample_test-hg38/counts_nuclei.rds"
# opts$counts_empty = "results/freeze1/decontx/Sample_test-hg38/counts_empty.rds"


### Read in data
matrix <- read10X(sample.dirs =c(opts$input_10x_dir), sample.names = c("dummy_name"), merge = FALSE)[[1]][[1]]

### Read in barcode lists
barcodes_nuclei = read.table(opts$barcodes_nuclei, header=TRUE)$barcode_gex
barcodes_empty = read.table(opts$barcodes_empty, header=TRUE)$barcode_gex
#barcodes_doublets = scan(opts$barcodes_doublets, what="character")

### Extract protein coding genes only
protein_coding_genes = read.table("resources/hg38/hg38.refGene.tss.bed.gz")$V4
protein_coding_genes = protein_coding_genes[protein_coding_genes %in% rownames(matrix)]
matrix_pc = matrix[protein_coding_genes,]

### Extract matrix for quality nuclei
#barcodes_to_keep = barcodes_nuclei[!barcodes_nuclei %in% barcodes_doublets]
barcodes_to_keep = barcodes_nuclei
matrix_nuclei = matrix_pc[,barcodes_to_keep]

### Extract matrix for empty droplets
barcodes_empty = barcodes_empty[barcodes_empty %in% colnames(matrix)]
matrix_empty = matrix_pc[,barcodes_empty]

### Save
saveRDS(matrix_nuclei, file=opts$counts_nuclei)
saveRDS(matrix_empty, file=opts$counts_empty)
