options(stringsAsFactors=FALSE)
library(optparse)
library(yaml)
library(dplyr)

option_list <- list(
  make_option(
    c("--barcode_to_cell_type"), type = "character", help = "Liger clusters collapsed into cell types"
  ),
  make_option(
    c("--barcode_map"), type = "character", help = "Multiome barcode GEX+ATAC map"
  ),
  make_option(
    c("--outdir"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# ## Testing
# opts = list()
# opts$barcode_to_cell_type = "results/liger/3GEX_noNM3/barcode_to_cluster_to_cell_type.csv"
# opts$barcode_map = "resources/multiome_barcode_map.rds"
# opts$outdir = "results/macs2/atac_barcodes_sample_cluster/noNM3"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format cluster to barcode map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
barcode_to_cell_type = read.csv(opts$barcode_to_cell_type)

barcode_map = readRDS(opts$barcode_map)
row.names(barcode_map) = barcode_map$GEX

barcode_to_cell_type$CB_atac = barcode_map[barcode_to_cell_type$CB_gex,"ATAC"]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract barcodes per sample per cluster
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sample_map = data.frame(nf_id = paste0("Sample_5124_NM_",seq(1:8),"_hg38"), samplename=paste0("Sample_5124-NM_NM-",seq(1:8)))
row.names(sample_map) = sample_map$nf_id


getBarcodesPerSamplePerCluster = function(d, sample, cell_type) {
  sub = d[d$SampleID==sample & d$cell_type==cell_type,]
  samplename = sample_map[sample,"samplename"]
  file = paste0("atac_barcodes_", samplename,"_", cell_type,".txt")
  write(sub$CB_atac, file=file.path(opts$outdir, file))
}

samples = unique(barcode_to_cell_type$SampleID)
cell_types = unique(barcode_to_cell_type$cell_type)

for (sample in samples) {
  for (type in cell_types) {
      cat("saving", sample, type,"\n")
      getBarcodesPerSamplePerCluster(barcode_to_cell_type, sample, type)
  }
}
