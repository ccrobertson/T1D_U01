library(rliger)
library(optparse)


# ## Testing
# opts = list()
# opts$liger_obj = "results/multiome/liger_clustering/liger_obj.rds"
# opts$outdir = "results/multiome/liger_clustering"

#Exploting liger object class
#getSlots("liger")
#This seems to be the combined cell metadata: liger_obj_lv@cell.data
#head(liger_obj_lv@cell.data)
#I think if we add new columns here indicating subject and condition, this will
#allow us to look at these comparisons as well


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in liger iNMF results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Reading in liger object")
liger_obj = readRDS(opts$liger_obj)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nuclei clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quantile normalization
message("Quantile normalization")
liger_obj <- quantile_norm(liger_obj)

## Louvain clustering
liger_obj_lv <- louvainCluster(liger_obj, resolution=0.25)

## UMAP
liger_obj_lv <- runUMAP(liger_obj_lv, distance = "cosine", n_neighbors = 30, min_dist = 0.3)

## Visualization
png(file.path(opts$outdir, "umap_clusters_and_batch.png"), width=900, height=500)
all.plots <- plotByDatasetAndCluster(liger_obj_lv, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
all.plots[[1]] + all.plots[[2]]
dev.off()

gene_loadings <- plotGeneLoadings(liger_obj_lv, do.spec.plot = FALSE, return.plots = TRUE)
pdf(file.path(opts$outdir, "umap_gene_loadings.pdf"), width=9, height=7)
gene_loadings[[1]]
gene_loadings[[2]]
gene_loadings[[3]]
gene_loadings[[4]]
gene_loadings[[5]]
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gene plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotMarkerGenesByDataset = function(gene) {
  plotlist = plotGene(liger_obj_lv, gene, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
  plot_grid(plotlist, ncol=2)
}
#plotMarkerGenesByDataset("INS")

## Inferred marker genes
cluster.results <- runWilcoxon(liger_obj_lv, compare.method = "clusters")
#inferred_markers = ??
pdf(file.path(opts$outdir, "umap_inferred_markers_by_dataset.pdf"))
for ( i in 1:length(inferred_markers)) {
    plotMarkerGenesByDataset(inferred_markers[i])
}
dev.off()


## Known marker genes
known_markers = c("INS", "GCG", "SST", "PPY", "KRT19", "PRSS1", "PTPRC", "VWF", "SDS", "RSGS5")
pdf(file.path(opts$outdir, "umap_known_markers_by_dataset.pdf"))
for ( i in 1:length(known_markers)) {
    plotMarkerGenesByDataset(known_markers[i])
}
dev.off()

## Any batch associated genes?
#batch.results <- runWilcoxon(liger_obj_lv, compare.method = "dataset")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create pseudobulk matrices?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
