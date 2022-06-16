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
    c("--clusters"), type = "character", help = "File mapping barcodes to clusters."
  ),
  make_option(
    c("--max_contamination"), type = "numeric", help = "Maximum contamination fraction for barcodes to return."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for decontX results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)



# ### Testing
# opts = list()
# opts$counts_nuclei = "results/multiome/decontx/Sample_5124-NM-1-hg38/counts_nuclei.rds"
# opts$counts_empty = "results/multiome/decontx/Sample_5124-NM-1-hg38/counts_empty.rds"
# opts$clusters = "results/multiome/seurat_prelim/Sample_5124-NM-1-hg38/seurat_clusters.csv"
# opts$outdir = "results/multiome/decontx/Sample_5124-NM-1-hg38"
# opts$max_contamination = 0.3


### Read data
counts_nuclei = readRDS(opts$counts_nuclei)
counts_empty = readRDS(opts$counts_empty)
clusters = read.csv(opts$clusters)

### Message
cat("START: Training decontX model using", ncol(counts_nuclei), "nuclei and", ncol(counts_empty), "empty droplets.\n")

### Convert to SingleCellExperiment object
cat("Converting to SingleCellExperiment object.\n")
x_nuclei = SingleCellExperiment(assays = list(counts = counts_nuclei))
x_empty = SingleCellExperiment(assays = list(counts = counts_empty))

### Make sure cluster labels are in right order
row.names(clusters) = clusters$barcode_gex
z = clusters[colnames(counts_nuclei),"clusters"]

### Run decontx -- unclear if the background option is actually being used here
# Bug was supposedly fixed in dev branch: https://github.com/campbio/celda/issues/355
cat("Running decontX.\n")
#results <- decontX(x = x_nuclei, background = x_empty, seed = 8675)
results <- decontX(x = x_nuclei,
  z = z,
  background = x_empty,
  delta=c(10,30),
  estimateDelta = FALSE,
  seed = 8675)

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

known_markers = c("INS", "GCG", "SST", "PPY", "KRT19", "PRSS1", "PTPRC", "VWF", "SDS", "RSGS5")
results <- logNormCounts(results)
png(file.path(opts$outdir, "umap_marker_genes_pre.png"))
plotDimReduceFeature(as.matrix(logcounts(results)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = known_markers,
    exactMatch = TRUE,
    useAssay = "counts"
  )
dev.off()

png(file.path(opts$outdir, "umap_marker_genes_post.png"))
plotDimReduceFeature(as.matrix(logcounts(results)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = known_markers,
    exactMatch = TRUE,
    useAssay = "decontXcounts"
  )
dev.off()


markers <- list(
    beta_INS = c("INS"),
    alpha_GCG = c("GCG"),
    delta_SST = c("SST"),
    gamma_PPY = c("PPY"),
    ductal_KRT19 = c("KRT19"),
    acinar_PRSS1 = c("PRSS1"),
    immune_CD45 = c("PTPRC"),
    endothelial_VWF = c("VWF"),
    unclear_SDS = c("SDS"),
    stellate_RGS5 = c("RGS5")
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Saving results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Saving complete decontX results.\n")
saveRDS(results, file=file.path(opts$outdir, "decontx_results.rds"))

cat("Saving updated count matrices.\n")
raw_counts_filtered = counts(x_nuclei)[,metadata(results)[[1]][["contamination"]] < opts$max_contamination]
decontaminated_counts_filtered = decontXcounts(results)[,metadata(results)[[1]][["contamination"]]< opts$max_contamination]
saveRDS(raw_counts_filtered, file=file.path(opts$outdir, "counts_low_contamination_raw.rds"))
saveRDS(decontaminated_counts_filtered, file=file.path(opts$outdir, "counts_low_contamination_decontaminated.rds"))


cat("DONE.\n")
