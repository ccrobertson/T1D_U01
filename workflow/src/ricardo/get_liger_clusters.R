options(stringsAsFactors = F)

library(rliger)
library(dplyr)
library(tidyr)

args <- commandArgs(T)

infile <- args[1]
outfile <- args[2]

message("Reading first-pass clusters")
ligerex <- readRDS(infile)

message("Writing clusters to output file")
df <- rliger::plotByDatasetAndCluster(ligerex, return.plots = T)[[1]]$data %>%
  tibble::rownames_to_column("barcode") %>%
  select(barcode, Cluster)

write.table(df, outfile, col.names = T, row.names = F, sep = "\t", quote = F)
