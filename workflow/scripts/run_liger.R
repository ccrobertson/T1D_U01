library(liger)
library(optparse)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in nuclei counts matrices across samples and uses liger
# to integrate across batches and cluster nuclei 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


option_list <- list(
  make_option(
    c("--counts_raw"), type = "character", help = "Count matrix BEFORE correcting for ambient contamination."
  ),
  make_option(
    c("--counts_clean"), type = "character", help = "Count matrix AFTER correcting for ambient contamination."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for liger results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)



# ### Testing
# opts = list()
# opts$counts_raw = "results/freeze1/decontx/Sample_5124-NM-2-hg38/counts_nuclei.rds"
# opts$counts_empty = "results/freeze1/decontx/Sample_5124-NM-2-hg38/counts_nuclei_decontaminated.rds"
# opts$outdir = "results/freeze1/liger/Sample_5124-NM-2-hg38"
