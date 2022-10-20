options(stringsAsFactors=FALSE)
library(dplyr)
library(optparse)


option_list <- list(
  make_option(
    c("--barcode_to_cell_type"), type = "character", help = "Based on Liger clustering."
  ),
  make_option(
    c("--run_design"), type = "character", help = " "
  ),
  make_option(
    c("--PSNG_threshold"), type = "numeric", help = "Min demuxlet PRB.SNG1 to include barcode in results."
  ),
  make_option(
    c("--prefix"), type = "character", help = "Output destination for liger results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser, positional_arguments=TRUE)
print(opts)

demuxlet_files = opts$args

opts <- opts$options

# ## Testing
# opts = list()
# opts$barcode_to_cell_type = "results/liger/multiome_5GEX/barcode_to_cluster_to_cell_type.csv"
# opts$run_design = "data/batch_design.tsv"
# opts$library_key = "data/library_key.tsv"
# opts$PSNG_threshold = 0.8
# opts$prefix = "results/master_barcode_maps/master_barcode_map_multiome_5GEX"
# demuxlet_files = c(
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-1/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-2/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-4/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-5/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-6/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-7/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-8/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-1/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-2/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-4/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-5/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-6/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-7/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/Sample_5125-NM-8/demuxlet.best")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get barcode to cell type map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Barcode to cell type map")
barcode_to_cell_type = read.csv(opts$barcode_to_cell_type)
barcode_to_cell_type$barcode = barcode_to_cell_type$CB_gex

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get barcode to donor map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Barcode to donor map")

processDemuxletFile = function(f) {
  #f = "results/demultiplex/demuxlet-unfiltered/Sample_5124-NM-1/demuxlet.best"
  cat("File:",f,"\n")
  donors = scan(gsub("demuxlet.best","chunked/demuxlet.aaaaaaaaaa.donors.txt", f, fixed=TRUE), what="character")
  vals = unlist(strsplit(f, split="/"))
  sample = vals[length(vals)-1]
  cat("Sample:",sample,"\n")
  cat("Donors:", donors, "\n")
  d = read.table(f, header=TRUE)
  d_cells = d %>%
    filter(BARCODE %in% barcode_to_cell_type$CB_gex,
      DROPLET.TYPE == "SNG",
      SNG.BEST.GUESS %in% donors)
  if (nrow(d_cells) > 0 ) {
    d_cells$Sample = sample
    d_cells$SampleID = gsub("-","_",sample)
  } else {
    cat("WARNING: No barcodes pass filter for", sample,"\n")
    d_cells = NULL
  }
  return(d_cells)
}

df = do.call("rbind", lapply(demuxlet_files, processDemuxletFile))

demuxlet_cols = names(df)
df$donor = df$SNG.BEST.GUESS
df$barcode = df$BARCODE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine into one map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Combine maps")
map = merge(df, barcode_to_cell_type, by=c("SampleID","barcode"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add condition info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Link to experimental condition")

library_to_batch = read.table(opts$library_key, header=TRUE)
row.names(library_to_batch) = library_to_batch$Master_SampleID

batch_design = read.table(opts$run_design, header=TRUE)
batch_design$unique_id = paste(batch_design$Batch_name, batch_design$Donor, sep="_")
row.names(batch_design) = batch_design$unique_id


map$batch = library_to_batch[map$Sample,"Batch"]
map$unique_id = paste(map$batch, map$donor, sep="_")

map$condition = batch_design[map$unique_id, "Condition"]
map$IEQ = batch_design[map$unique_id, "IEQ"]
map$unique_barcode = paste(map$Sample, map$barcode, sep="_")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Explore
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(map$condition)
table(map$condition, map$donor)
table(map$condition, map$batch)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Save")
demuxlet_cols_to_keep = demuxlet_cols[-1][1:(length(demuxlet_cols)-3)]
out = map %>%
  rename(sample = Sample) %>%
  select(unique_barcode, sample, barcode, donor, batch, condition, cell_type, unique_id, IEQ, all_of(demuxlet_cols_to_keep))

saveRDS(out, file=paste0(opts$prefix,".rds"))
write.table(out, file=paste0(opts$prefix,".tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
