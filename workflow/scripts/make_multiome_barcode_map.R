options(stringsAsFactors=FALSE)

map = data.frame(
  GEX=scan("resources/barcode_whitelist_multiome_GEX.txt", what="character"),
  ATAC=scan("resources/barcode_whitelist_multiome_ATAC.txt", what="character")
)
saveRDS(map, "resources/multiome_barcode_map.rds")
write.table(map, "resources/multiome_barcode_map.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

