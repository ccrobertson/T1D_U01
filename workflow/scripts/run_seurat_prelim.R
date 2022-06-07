options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(ggplot2)
theme_set(theme_bw(base_size = 12))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in raw counts matrices for separately for each batch
# and generates preliminary cell type clusters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


option_list <- list(
  make_option(
    c("--counts"), type = "character", help = "Input quality nuclei counts matrix."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for clustering results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# ### Testing
# opts = list()
# opts$counts = "results/freeze1/decontx/Sample_5124-NM-1-hg38/counts_nuclei.rds"
# opts$outdir = "results/freeze1/seurat_prelim/Sample_5124-NM-1-hg38"

# opts = list()
# opts$counts = "results/freeze1/decontx/Sample_test-hg38/counts_nuclei.rds"
# opts$outdir = "results/freeze1/seurat_prelim/Sample_test-hg38"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts = readRDS(opts$counts)

### Create Seurat object
#sce = SingleCellExperiment(assays = list(counts = counts))
#sce <- logNormCounts(sce)
#sobj = as.Seurat(sce, counts = "counts", data = "logcounts")
sobj <- CreateSeuratObject(counts = counts, project = "multiome", min.cells = 0, min.features = 0)
sobj

#vignette on converting from SingleCellExperiment to Seurat - https://satijalab.org/seurat/articles/conversion_vignette.html
#x_seurat <- as.Seurat(x_nuclei, counts = "counts", data = "logcounts")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Calculte QC metrics (note these counts are just protein coding genes, so MT percentage isn't meaningful)
#sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
sobj[["percent.INS"]] <- PercentageFeatureSet(sobj, pattern = "^INS$")
sobj[["percent.GCG"]] <- PercentageFeatureSet(sobj, pattern = "^GCG$")
head(sobj@meta.data, 5)
#VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.INS", "percent.GCG"), ncol = 4)

png(file.path(opts$outdir,"seurat_density_insulin.png"))
ggplot(sobj@meta.data) + geom_density(aes(x=percent.INS+0.01)) + scale_x_continuous(trans='log10') + xlab("Percent INS")
dev.off()

png(file.path(opts$outdir,"seurat_density_glucagon.png"))
ggplot(sobj@meta.data) + geom_density(aes(x=percent.GCG+0.01)) + scale_x_continuous(trans='log10') + xlab("Percent GCG")
dev.off()

# png(file.path(opts$outdir,"seurat_qc_scatter.png"))
# plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1
# dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normalization and scaling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Normalize
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

### Find highly variable genes
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sobj), 10)

# plot variable features with and without labels
png(file.path(opts$outdir,"seurat_variable_features.png"), width=800, height=400)
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

### Scale variable features
variable.genes = VariableFeatures(sobj)
all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = variable.genes)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
print(sobj[["pca"]], dims = 1:5, nfeatures = 5)

png(file.path(opts$outdir,"seurat_pca_dimloadings.png"), width=800, height=400)
VizDimLoadings(sobj, dims = 1:2, reduction = "pca")
dev.off()

png(file.path(opts$outdir,"seurat_pca_scatter.png"), width=600, height=400)
DimPlot(sobj, reduction = "pca")
dev.off()

pdf(file.path(opts$outdir,"seurat_pca_heatmap.pdf"))
DimHeatmap(sobj, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(sobj, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(sobj, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(sobj, dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(sobj, dims = 25:30, cells = 500, balanced = TRUE)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Determining number of PCs to use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#sobj <- JackStraw(sobj, num.replicate = 100)
# sobj <- ScoreJackStraw(sobj, dims = 1:20)
#
# png(file.path(opts$outdir,"seurat_jackstraw.png"), width=600, height=400)
# JackStrawPlot(sobj, dims = 1:20)
# dev.off()

png(file.path(opts$outdir,"seurat_elbow.png"), width=600, height=400)
ElbowPlot(sobj)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clustering nuclei
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sobj <- FindNeighbors(sobj, dims = 1:20)
sobj <- FindClusters(sobj, resolution = 0.5)

#view clusters for first 5 nuclei
head(Idents(sobj), 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UMAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sobj <- RunUMAP(sobj, dims = 1:20)

png(file.path(opts$outdir,"seurat_umap.png"), width=600, height=400)
DimPlot(sobj, reduction = "umap")
dev.off()

saveRDS(sobj, file = file.path(opts$outdir,"seurat_obj.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find/view marker genes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#sobj.markers <- FindAllMarkers(sobj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
# head(sobj.markers, n = 5)
#
# sobj.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)
#

known_markers = c("INS", "GCG", "SST", "KRT19", "PRSS1", "PTPRC")

png(file.path(opts$outdir,"seurat_markers_known_violin.png"), width=600, height=400)
VlnPlot(sobj, features = known_markers)
dev.off()

png(file.path(opts$outdir,"seurat_markers_known_umap.png"), width=600, height=600)
FeaturePlot(sobj, features = known_markers)
dev.off()

# png(file.path(opts$outdir,"seurat_markers_heatmap.png"), width=600, height=400)
# sobj.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(sobj, features = top10$gene) + NoLegend()
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign cell types to clusters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new.cluster.ids <- c("alpha", "beta", "delta", "gamma",
#   "acinar", "ductal", "immune",
#   "endothelial", "stellate", "other")
# names(new.cluster.ids) <- levels(sobj)
# sobj <- RenameIdents(sobj, new.cluster.ids)
#
# png(file.path(opts$outdir,"seurat_umap_annotated.png"), width=600, height=400)
# DimPlot(sobj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# dev.off()
