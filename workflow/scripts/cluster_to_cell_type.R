options(stringsAsFactors=FALSE)
library(optparse)
library(yaml)
library(dplyr)

option_list <- list(
  make_option(
    c("--barcode_to_cluster"), type = "character", help = "Liger clusters tsv"
  ),
  make_option(
    c("--deg"), type = "character", help = "Liger DEG by cluster rds"
  ),
  make_option(
    c("--markers"), type = "character", help = "Liger DEG by cluster rds"
  ),
  make_option(
    c("--outdir"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# ### Testing
# opts = list()
# opts$barcode_to_cluster = "results/liger/3GEX_all/liger_clusters.tsv"
# opts$deg = "results/liger/3GEX_all/liger_deg_by_cluster.rds"
# opts$markers = "workflow/src/markers.yaml"
# opts$outdir = "results/liger/3GEX_all"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format cluster to barcode map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
barcode_to_cluster = read.table(opts$barcode_to_cluster, header=TRUE)

splitSampleIdBarcode = function(x) {
  vals = unlist(strsplit(x, split="_"))
  barcode = vals[length(vals)]
  sample = paste0(vals[1:(length(vals)-1)], collapse="_")
  return(c(barcode, sample))
}

barcode_to_cluster$CB_gex = t(sapply(barcode_to_cluster$barcode, splitSampleIdBarcode))[,1]
barcode_to_cluster$SampleID = t(sapply(barcode_to_cluster$barcode, splitSampleIdBarcode))[,2]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create liger cluster to cell type map
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clusters = sort(unique(barcode_to_cluster$Cluster))
deg = readRDS(opts$deg)
cell_type_markers = read_yaml(opts$markers)$cell_types


### Assing primary clusters - get top cluster per marker
getTopEnrichedClustersPerMarker = function(marker, deg, n=10) {
  d = deg[deg$feature==marker,]
  sorted = d[order(d$logFC, decreasing=TRUE),]
  return(sorted[1:n,])
}
#lapply(cell_type_markers, getTopEnrichedClustersPerMarker, deg=deg, n=10)
primary_clusters = do.call("rbind", lapply(cell_type_markers, getTopEnrichedClustersPerMarker, deg=deg, n=1))
primary_clusters$cell_type = row.names(primary_clusters)

#if two markers have highest enrichment for the same cluster, assign cluster
# to cell_type with highest logFC
primary_clusters_ordered = primary_clusters[order(primary_clusters$logFC, decreasing=TRUE),]
primary_clusters_dedup = primary_clusters_ordered[!duplicated(primary_clusters_ordered$group),]

### Assign secondary clusters - get sig markers per cluster
getSigMarkersPerCluster = function(cluster, deg) {
    d = deg[deg$group==cluster & deg$feature %in% cell_type_markers & deg$logFC>1 & deg$padj<0.01,]
    sorted = d[order(d$logFC, decreasing=TRUE),]
    return(sorted[1,])
}
secondary_clusters = do.call("rbind", lapply(clusters[!clusters %in% primary_clusters$group], getSigMarkersPerCluster, deg=deg))
secondary_clusters$cell_type = sapply(secondary_clusters$feature, function(x) {names(cell_type_markers)[cell_type_markers==x]})

### Combine primary and secondary
map = rbind(primary_clusters_dedup, secondary_clusters)
map = map[order(map$cell_type),]

row.names(map) = map$group
barcode_to_cluster$cell_type = map[as.character(barcode_to_cluster$Cluster), "cell_type"]

table(barcode_to_cluster$cell_type)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(map[,c("group", "cell_type", "feature","logFC","padj")], file=file.path(opts$outdir, "cluster_to_cell_type.csv"), row.names=FALSE, quote=FALSE)
write.csv(barcode_to_cluster[,c("SampleID", "CB_gex", "Cluster","cell_type")], file=file.path(opts$outdir, "barcode_to_cluster_to_cell_type.csv"), row.names=FALSE, quote=FALSE)


# ### Get top genes per cluster
# getTopFeaturesPerCluster = function(cluster, deg, n=10) {
#     d = deg[deg$group==cluster,]
#     sorted = d[order(d$logFC, decreasing=TRUE),]
#     return(sorted[1:n,])
# }
# lapply(clusters, getTopFeaturesPerCluster, deg=deg, n=10)
