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
# opts$barcode_to_cell_type = "results/liger/multimodal_noNM3/barcode_to_cluster_to_cell_type.csv"
# opts$run_design = "data/nandini_run_design.tsv"
# opts$PSNG_threshold = 0
# opts$prefix = "results/master_barcode_map/master_barcode_map_multimodal_noNM3"
# demuxlet_files = c(
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-1/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-2/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-4/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-5/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-6/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-7/demuxlet.best",
#   "results/demultiplex/demuxlet-unfiltered/5124-NM-8/demuxlet.best")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get barcode to cell type map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Barcode to cell type map")
barcode_to_cell_type = read.csv(opts$barcode_to_cell_type)

sample_map_nf = data.frame(nf_id = paste0("Sample_5124_NM_",seq(1:8),"_hg38"), samplename=paste0("Sample_5124-NM_NM-",seq(1:8)))
row.names(sample_map_nf) = sample_map_nf$nf_id

barcode_to_cell_type$SampleID2 = sample_map_nf[barcode_to_cell_type$SampleID, "samplename"]
barcode_to_cell_type$SampleID2[is.na(barcode_to_cell_type$SampleID2)] <- barcode_to_cell_type$SampleID[is.na(barcode_to_cell_type$SampleID2)]

barcode_to_cell_type$barcode = barcode_to_cell_type$CB_gex

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get barcode to donor map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Barcode to donor map")

processDemuxletFile = function(f) {
  cat("File:",f,"\n")
  donors = scan(gsub(".best",".donor_list.txt", f, fixed=TRUE), what="character")
  vals = unlist(strsplit(f, split="/"))
  sample = vals[length(vals)-1]
  cat("Sample:",sample,"\n")
  d = read.table(f, header=TRUE)
  d_cells = d[d$BARCODE %in% barcode_to_cell_type$CB_gex & d$BEST %in% paste0("SNG-",donors) & d$PRB.SNG1>opts$PSNG_threshold,]
  d_cells$SampleID = sample
  return(d_cells)
}

df = do.call("rbind", lapply(demuxlet_files, processDemuxletFile))

demuxlet_cols = names(df)

df$donor = sapply(df$BEST, function(x) gsub("SNG-","", x))

sample_map_demuxlet = data.frame(demux_id = paste0("5124-NM-",seq(1:8)), samplename=paste0("Sample_5124-NM_NM-",seq(1:8)))
row.names(sample_map_demuxlet) = sample_map_demuxlet$demux_id

df$SampleID2 = sample_map_demuxlet[df$SampleID,"samplename"]
df$SampleID2[is.na(df$SampleID2)] <- df$SampleID[is.na(df$SampleID2)]

df$barcode = df$BARCODE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine into one map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Combine maps")
map = merge(df, barcode_to_cell_type, by=c("SampleID2","barcode"))

map$run = sapply(map$SampleID2, function(x) unlist(strsplit(x, split="_"))[2])
map$batch = sapply(map$SampleID2, function(x) unlist(strsplit(x, split="_"))[3])
map$unique_id = paste(map$batch,map$donor, sep="_")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add condition info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Link to experimental condition")

design = read.table(opts$run_design, header=TRUE)
design$unique_id = paste(design$Batch_name, design$Subject_ID, sep="_")
row.names(design) = design$unique_id

design$condition = sapply(design$Condition, function(x) unlist(strsplit(x, split="_"))[1])
design$IEQ = as.numeric(sapply(design$Condition, function(x) unlist(strsplit(x, split="_"))[2]))

map$condition = design[map$unique_id, "condition"]
map$IEQ = design[map$unique_id, "IEQ"]

map$unique_barcode = paste(map$SampleID2, map$barcode, sep="_")

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
demuxlet_cols_to_keep = demuxlet_cols[-1][1:(length(demuxlet_cols)-2)]
out = map %>%
  rename(sample = SampleID2) %>%
  rename(demuxlet_sample_id = SampleID.x) %>%
  rename(nf_sample_id = SampleID.y) %>%
  select(unique_barcode, sample, barcode, run, donor, batch, condition, cell_type, unique_id, IEQ, all_of(demuxlet_cols_to_keep))

saveRDS(out, file=paste0(opts$prefix,".rds"))
write.table(out, file=paste0(opts$prefix,".tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
