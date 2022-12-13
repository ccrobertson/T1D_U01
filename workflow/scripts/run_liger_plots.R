options(stringsAsFactors=FALSE)
library(optparse)
library(rliger)



# # ## Testing
# opts = list()
# opts$liger_obj = "results/multiome/liger_clustering/liger_obj_clusters.rds"
# opts$outdir = "results/multiome/liger_clustering"

#Exploring liger object class
#getSlots("liger")
#This seems to be the combined cell metadata: liger_obj_lv@cell.data
#head(liger_obj_lv@cell.data)
#I think if we add new columns here indicating subject and condition, this will
#allow us to look at these comparisons as well

option_list <- list(
  make_option(
    c("--liger_obj"), type = "character", help = "Liger iNMF results."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for liger results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

## Show how input arguments were processed
print(opts)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in liger iNMF+Louvain clustering results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Reading in liger object")
liger_obj = readRDS(opts$liger_obj)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Marker gene plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nsamples = length(liger_obj@raw.data)

## Known marker genes
known_markers = c("INS", "GCG", "SST", "PPY", "KRT19", "PRSS1", "PTPRC", "VWF", "SDS", "RGS5")
#pdf(file.path(opts$outdir, "umap_known_markers_by_dataset.pdf"))
for ( i in 1:length(known_markers)) {
    gene = known_markers[i]
    cat("Plotting", gene, "\n")
    png(file.path(opts$outdir, paste0("umap_",gene,".png")), width=1000, height=250*nsamples)
    plist = plotGene(liger_obj, gene, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
    p = plot_grid(plotlist=plist, ncol=2)
    print(p)
    dev.off()
}
#dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export UMAP coordinates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df = data.frame(liger_obj@cell.data, cluster=liger_obj@clusters, liger_obj@tsne.coords)
saveRDS(df, file=file.path(opts$outdir, "liger_umap_coordinates.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Differential expression analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Find cluster-specific genes
cluster_genes <- runWilcoxon(liger_obj, compare.method = "clusters")
saveRDS(cluster_genes, file=file.path(opts$outdir, "liger_deg_by_cluster.rds"))

## Find batch-specific genes
batch_genes <- runWilcoxon(liger_obj, compare.method = "datasets")
saveRDS(batch_genes, file=file.path(opts$outdir, "liger_deg_by_batch.rds"))
