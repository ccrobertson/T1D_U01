options(stringsAsFactors=FALSE)
library(optparse)
library(rliger)



# ## Testing
# opts = list()
# opts$liger_obj = "results/multiome/liger_clustering/liger_obj.rds"
# opts$outdir = "results/multiome/liger_clustering"

#Exploring liger object class
#getSlots("liger")
#This seems to be the combined cell metadata: liger_obj@cell.data
#head(liger_obj@cell.data)
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
liger_obj <- louvainCluster(liger_obj, resolution=0.25)

## UMAP
liger_obj <- runUMAP(liger_obj, distance = "cosine", n_neighbors = 30, min_dist = 0.3)
saveRDS(liger_obj, file=file.path(opts$outdir, "liger_obj_clusters.rds"))

## Visualization
png(file.path(opts$outdir, "umap_clusters_and_batch.png"), width=1100, height=500)
all.plots <- plotByDatasetAndCluster(liger_obj, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = TRUE)
all.plots[[1]] + all.plots[[2]]
dev.off()

# gene_loadings <- plotGeneLoadings(liger_obj, do.spec.plot = FALSE, return.plots = TRUE)
# pdf(file.path(opts$outdir, "umap_gene_loadings.pdf"), width=9, height=7)
# gene_loadings[[1]]
# gene_loadings[[2]]
# gene_loadings[[3]]
# gene_loadings[[4]]
# gene_loadings[[5]]
# dev.off()
