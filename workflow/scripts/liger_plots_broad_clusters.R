options(stringsAsFactors=FALSE)
library(optparse)
library(rliger)
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 12))

option_list <- list(
  make_option(
    c("--liger_obj"), type = "character", help = "Liger clustering results."
  ),
  make_option(
    c("--barcode_to_cell_type"), type = "character", help = "Liger iNMF results."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for liger results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# ## Testing
# opts = list()
# opts$liger_obj = "results/liger/multimodal_noNM3/liger_obj_clusters.rds"
# opts$barcode_to_cell_type = "results/liger/multimodal_noNM3/barcode_to_cluster_to_cell_type.csv"
# opts$outdir = "results/plots/multimodal_noNM3"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overview of liger class
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# liger_obj@cell.data - data.frame containing meta data about cells
# liger_obj@clusters - vector (factor)
# liger_obj@tsne.coords - UMAP coordinates

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
liger_obj = readRDS(opts$liger_obj)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to data.frame
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df = data.frame(liger_obj@cell.data, cluster=liger_obj@clusters, liger_obj@tsne.coords)

map = read.csv(opts$barcode_to_cell_type)
row.names(map) = paste(map$SampleID, map$CB_gex, sep="_")


df$Cluster = map[row.names(df),"Cluster"]
df$cell_type = map[row.names(df),"cell_type"]

#sanity check
if (sum(!df$cluster==df$Cluster)>0) {
  message("ERROR files do not match")
  stop()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot UMAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df$modality = NA
df$modality[df$dataset == "Sample_5125_NM_1_NM_1"] <- "5GEX"
df$modality[df$dataset %in% c('Sample_5124_NM_1_hg38',
            'Sample_5124_NM_2_hg38','Sample_5124_NM_4_hg38',
            'Sample_5124_NM_5_hg38','Sample_5124_NM_6_hg38',
            'Sample_5124_NM_7_hg38','Sample_5124_NM_8_hg38')] <- "Multiome"


## Plot UMAP
#png(file.path(opts$outdir, "umap_cell_type.png"), width=1000, height=500)
pdf(file.path(opts$outdir, "umap_cell_type.pdf"), width=12, height=5)
ggplot(df, aes(x = X1, y = X2, color = factor(cell_type))) +
  geom_point(size = 0.3, alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 10, alpha=1))) +
  #theme(legend.key.size = 10, legend.key.alpa=1) +
  facet_wrap(~ factor(modality, levels=c("Multiome","5GEX"), labels=c("3' RNA (Multiome)", "5' RNA"))) +
  xlab("UMAP 1") + ylab("UMAP 2") + scale_color_discrete(name="")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot cell type proportions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

map_noNM3 = readRDS("results/master_barcode_map_3GEX_noNM3.rds")
map_all = readRDS("results/master_barcode_map_3GEX_all.rds")


#### Frequencies by donor
freqtable = map_noNM3 %>%
  group_by(donor, cell_type) %>%
  summarize(count=n())

total_cells = freqtable %>%
  group_by(donor) %>%
  summarize(total = sum(count))

freqtable = left_join(freqtable, total_cells)
freqtable$prop = freqtable$count/freqtable$total

freqtable$donor = factor(freqtable$donor,
  levels=c( "HPAP093" ,"HPAP105", "ICRH139" ,"ICRH142", "ICRH143", "HPAP107"),
  labels=c( "HPAP093" ,"HPAP105", "ICRH139" ,"ICRH142", "ICRH143", "HPAP107 (AAb+)"))

pdf(file.path(opts$outdir, 'barplots_by_donor.pdf'))
ggplot(freqtable, aes(x=donor, y=count, fill=cell_type)) +
  geom_bar(position="stack", stat="identity") +
  ylab("Number of nuclei") + xlab("Islet Donor") + scale_fill_discrete(name="") +
  theme(axis.text.x = element_text(size=16, angle=45, hjust=1), axis.title=element_text(size=20))

ggplot(freqtable, aes(x=donor, y=prop, fill=cell_type)) +
  geom_bar(position="stack", stat="identity") +
  ylab("Proportion of nuclei") + xlab("Islet Donor") + scale_fill_discrete(name="") +
  theme(axis.text.x = element_text(size=16, angle=45, hjust=1), axis.title=element_text(size=20))

dev.off()


##### Frequencies by condition
freqtable2 = map_noNM3 %>%
  group_by(condition, cell_type) %>%
  summarize(count=n())

total_cells2 = freqtable2 %>%
  group_by(condition) %>%
  summarize(total = sum(count))

freqtable2 = left_join(freqtable2, total_cells2)
freqtable2$prop = freqtable2$count/freqtable2$total

# freqtable2$condition = factor(freqtable2$condition,
#   levels=c("control", "cytokines", "mock", "infection", "AAb"),
#   labels = c("Control", "Cytokines", "Mock Infection", "CVB4 Infection", "AAb+")
# )

pdf(file.path(opts$outdir, 'barplots_by_condition.pdf'))
ggplot(freqtable2, aes(x=condition, y=count, fill=cell_type)) +
  geom_bar(position="stack", stat="identity") +
  ylab("Number of nuclei") + xlab("Treatment Condition") + scale_fill_discrete(name="") +
  theme(axis.text.x = element_text(size=16, angle=45, hjust=1), axis.title=element_text(size=20))

ggplot(freqtable2, aes(x=condition, y=prop, fill=cell_type)) +
  geom_bar(position="stack", stat="identity") +
  ylab("Proportion of nuclei") + xlab("Treatment Condition")+ scale_fill_discrete(name="") +
  theme(axis.text.x = element_text(size=16, angle=45, hjust=1), axis.title=element_text(size=20))

dev.off()
