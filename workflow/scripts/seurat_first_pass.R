# Following the introductory tutorial here:
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(dplyr)
library(Seurat)
library(patchwork)

# # Load the PBMC dataset
# pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# # Initialize the Seurat object with the raw (non-normalized data).
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# pbmc

# Load test sample
#datadir = "/lab/work/ccrober/T1D_U01/analysis/freeze1/nf_gex_test_results/starsolo/Sample_test-hg38/Sample_test-hg38.Solo.out/Gene/raw"
run_name = "nf_gex"
sample_name = "Sample_5124-NM-1-hg38"
datadir = file.path("/lab/work/ccrober/T1D_U01/analysis/freeze1", paste0(run_name,"_results"), "starsolo", sample_name, paste0(sample_name, ".Solo.out/Gene/raw"))

# --soloFeatures Transcript3p Gene GeneFull SJ Velocyto
# #Gene = reads match the gene transcript
# #SJ = splice junctions
# #GeneFull = pre-mRNA (exons and introns)


seurobj.data = ReadSTARsolo(datadir)
seurobj = CreateSeuratObject(counts = seurobj.data, project = "multiome_gex", min.cells = 3, min.features = 200)


# Calculate mitochondrial read percentage
seurobj[["percent.mt"]] <- PercentageFeatureSet(seurobj, pattern = "^MT-")
VlnPlot(seurobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#nFeature_RNA = number of unique genes/features per cell
#nCount_RNA = number of total molecules per cell

#Following filtering criteria in Vivek's manuscript
#We filter cells that have total molecules over 1000 (should there be an upper bound?)
#We filter cells that have >5% mitochondrial counts
seurobj_filtered <- subset(seurobj, subset = nCount_RNA > 1000 & percent.mt < 5)
plot1 <- FeatureScatter(seurobj_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurobj_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
