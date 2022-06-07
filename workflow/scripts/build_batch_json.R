#!/usr/bin/env python
library(rjson)
library(optparse)

option_list <- list(
  make_option(
    c("--batchfile"), type = "character", help = "File defining sample pools and conditions in each batch."
  )
)


option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

## Testing
#opts = list()
#opts$batchfile = "Sample_islet_list_for_multiomics_Batches_long_format_with_libname.tsv"

batch_info = read.table(file=opts$batchfile, header=TRUE)

convertToList = function(x) {
  v = x$Condition
  names(v) = x$Subject_ID
  as.list(v)
}

convertBatchDFtoJSON = function(x) {
  x_by_batch = split(x, f=x$Libname)
  x_by_batch_by_subject = lapply(x_by_batch, FUN=convertToList)
  x_json = toJSON(list(batches=x_by_batch_by_subject), indent=4, method="C")
  return(x_json)
}

batch_json = convertBatchDFtoJSON(batch_info)
write(batch_json, file=gsub(".tsv$", ".json", opts$batchfile, perl=TRUE))
