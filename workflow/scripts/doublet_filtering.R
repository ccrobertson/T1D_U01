options(stringsAsFactors=FALSE)
library(optparse)
library(Seurat)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script visualizes barcodes flagged as doublets by AMULET
# relative to provided Seurat clusters and then filters the
# provided counts matrix to remove the flagged doublet barcodes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


option_list <- list(
  make_option(
    c("--seurat_obj"), type = "character", help = "Preliminary seurat clustering object."
  ),
  make_option(
    c("--doublets"), type = "numeric", help = "List of barcodes flagged as doublets by AMULET."
  ),
  make_option(
    c("--map"), type = "numeric", help = "File to use to map atac barcodes to gex barcodes."
  ),
  make_option(
    c("--input_counts"), type = "character", help = "Input gex counts matrix."
  ),
  make_option(
    c("--output_counts"), type = "character", help = "Out gex counts matrix after removing doublets."
  ),
  make_option(
    c("--plotfile"), type = "character", help = "Seurat UMAP showing doublets."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# opts = list()
# opts$seurat_prelim = "results/multiome/seurat_prelim/Sample_5124-NM-2-hg38/seurat_obj.rds"
# #opts$seurat_round3 = "results/multiome/seurat_round3/Sample_5124-NM-2-hg38/seurat_obj.rds"
# opts$qc = "results/multiome/cross_modality_qc/Sample_5124-NM-2-hg38/cross_modality_qc.txt"
# opts$doublets = "results/multiome/amulet/Sample_5124-NM-2-hg38/MultipletBarcodes_01.txt"
# opts$outdir = "results/multiome/amulet/Sample_5124-NM-2-hg38"


#get doublet gex barcodes
doublets_barcode_atac = scan(opts$doublets, what="character")
map = read.table(opts$map, header=TRUE)
doublets_barcode_gex = map$barcode_gex[map$barcode_atac %in% doublets_barcode_atac]


#get Seurat objects pre and post decontamination
sobj = readRDS(opts$seurat_obj)


#UMAP highlighting doublets decontx
png(opts$plotfile)
DimPlot(sobj, reduction="umap", cells.highlight=doublets_barcode_gex)
dev.off()

#remove flagged doublets from count matrix
d = readRDS(opts$input_counts)
d2 = d[,!colnames(d) %in% doublets_barcode_gex]
saveRDS(d2, file=opts$output_counts)
