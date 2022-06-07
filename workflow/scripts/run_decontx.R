options(stringsAsFactors=FALSE)
library(optparse)
library(celda)
library(SingleCellExperiment)
library(scater)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in raw counts matrices (dgCMatrix) for nuclei and empty droplets
# and runs decontx to generate decontaminated counts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


option_list <- list(
  make_option(
    c("--counts_nuclei"), type = "character", help = "Input file name for storing quality nuclei counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--counts_empty"), type = "character", help = "Input file name for storing empty droplets counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for decontX results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)



# ### Testing
#opts = list()
#opts$counts_nuclei = "results/freeze1/decontx/Sample_5124-NM-1-hg38/counts_nuclei.rds"
#opts$counts_empty = "results/freeze1/decontx/Sample_5124-NM-1-hg38/counts_empty.rds"
#opts$outdir = "results/freeze1/decontx/Sample_5124-NM-1-hg38"

# opts = list()
# opts$counts_nuclei = "results/freeze1/decontx/Sample_test-hg38/counts_nuclei.rds"
# opts$counts_empty = "results/freeze1/decontx/Sample_test-hg38/counts_empty.rds"
# opts$outdir = "results/freeze1/decontx/Sample_test-hg38"


### Read data
counts_nuclei = readRDS(opts$counts_nuclei)
counts_empty = readRDS(opts$counts_empty)

### Message
cat("START: Training decontX model using", ncol(counts_nuclei), "nuclei and", ncol(counts_empty), "empty droplets.\n")

### Convert to SingleCellExperiment object
cat("Converting to SingleCellExperiment object.\n")
x_nuclei = SingleCellExperiment(assays = list(counts = counts_nuclei))
x_empty = SingleCellExperiment(assays = list(counts = counts_empty))

### Run decontx -- unclear if the background option is actually being used here
# Bug was supposedly fixed in dev branch: https://github.com/campbio/celda/issues/355
cat("Running decontX.\n")
#results <- decontX(x = x_nuclei, background = x_empty, seed = 8675)
results <- decontX(x = x_nuclei, delta=c(10,100), estimateDelta = FALSE, seed = 8675)

### UMAP plots
cat("Generating UMAP plots.\n")
umap <- reducedDim(results, "decontX_UMAP")

png(file.path(opts$outdir, "umap_clusters.png"))
plotDimReduceCluster(x = results$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])
dev.off()

png(file.path(opts$outdir, "umap_contamination.png"))
plotDecontXContamination(results)
dev.off()

#INS = beta
#GCG = alpha
#SST = delta
#KRT19 = ductal
#CPA1 = acinar
#PTPRC/CD45 = immune
results <- logNormCounts(results)
png(file.path(opts$outdir, "umap_marker_genes_pre.png"))
plotDimReduceFeature(as.matrix(logcounts(results)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = c("INS", "GCG","SST", "KRT19", "CPA1", "PTPRC", "HLA-A", "HLA-DPB1"),
    exactMatch = TRUE,
    useAssay = "counts"
  )
dev.off()

png(file.path(opts$outdir, "umap_marker_genes_post.png"))
plotDimReduceFeature(as.matrix(logcounts(results)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = c("INS", "GCG", "SST", "KRT19", "PRSS1", "PTPRC", "HLA-A", "HLA-DPB1"),
    exactMatch = TRUE,
    useAssay = "decontXcounts"
  )
dev.off()


markers <- list(
    beta_INS = c("INS"),
    alpha_GCG = c("GCG"),
    delta_SST = c("SST"),
    ductal_KRT19 = c("KRT19"),
    acinar_PRSS1 = c("PRSS1"),
    immune_CD45 = c("PTPRC")
  )
png(file.path(opts$outdir, "barplots_marker_genes_pre.png"))
plotDecontXMarkerPercentage(results,
    markers = markers,
    assayName = "counts")
dev.off()

png(file.path(opts$outdir, "barplots_marker_genes_post.png"))
plotDecontXMarkerPercentage(results,
    markers = markers,
    assayName = "decontXcounts")
dev.off()



### Saving results
cat("Saving results.\n")
saveRDS(results, file=file.path(opts$outdir, "decontx_results.rds"))
saveRDS(decontXcounts(results), file=file.path(opts$outdir, "counts_nuclei_decontaminated.rds"))



cat("DONE.\n")
