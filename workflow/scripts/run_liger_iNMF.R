library(rliger)
library(optparse)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in nuclei counts matrices across samples and uses liger
# to integrate across batches and cluster nuclei
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


option_list <- list(
  make_option(
    c("--outdir"), type = "character", help = "Output destination for liger results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser, positional_arguments=TRUE)

## Show how input arguments were processed
print(opts)

## Testing
# opts = list()
# opts$outdir = "results/multiome/liger_clustering"
# opts$args = c("results/multiome/decontx_round2/Sample_5124-NM-1-hg38/counts_low_contamination_decontaminated.rds",
#   "results/multiome/decontx_round2/Sample_5124-NM-2-hg38/counts_low_contamination_decontaminated.rds")

## Read in count matrices and add sample names to barcodes
count_matrix_files = opts$args
getCounts = function(x) {
  sample = gsub("-","_",sub(".*/([\\w-]+)/[\\w-]+.rds$","\\1", x, perl=TRUE))
  cat("Loading count matrix for", sample, "\n")
  d = readRDS(x)
  colnames(d) <- paste0(sample, "_",colnames(d))
  return(list(sample, d))
}

counts = list()
for (i in 1:length(count_matrix_files)) {
  filename = count_matrix_files[i]
  out = getCounts(filename)
  sample = out[[1]]
  counts[[sample]] = out[[2]]
}


## Create Liger object
cat("Creating liger object\n")
liger_obj = createLiger(counts)

## Preprocessing
cat("Preprocessing\n")
liger_obj <- normalize(liger_obj)
liger_obj <- selectGenes(liger_obj)
liger_obj <- scaleNotCenter(liger_obj)

## iNMF
cat("Running matrix integration (iNMF)\n")
liger_obj <- optimizeALS(liger_obj, k = 20)

## Save
cat("Saving liger object\n")
saveRDS(liger_obj, file = file.path(opts$options$outdir, "liger_obj.rds"))
