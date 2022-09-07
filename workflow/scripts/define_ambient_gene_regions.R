options(stringsAsFactors=FALSE)
library(optparse)
library(rtracklayer)
library(ggplot2)
theme_set(theme_bw(base_size = 12))
library(reshape2)


option_list <- list(
  make_option(
    c("--counts_empty"), type = "character", help = "Input file name for storing empty droplets counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--counts_nuclei"), type = "character", help = "Input file name for storing nuclei droplets counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--gtf"), type = "character", help = "Genome annotation file used for generating count matrix."
  ),
  make_option(
    c("--p"), type = "character", help = "Filter out top p% of highest expressed genes in empty droplets. If set to auto, will determine p using empty droplet and nuclei count matrices."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for results."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

### Testing
# opts = list()
# opts$counts_empty =  "results/multiome/counts_by_sample_gex/Sample_5124-NM-2-hg38/counts_empty.rds"
# opts$counts_nuclei =  "results/multiome/counts_by_sample_gex/Sample_5124-NM-2-hg38/counts_nuclei.rds"
# opts$gtf = "resources/hg38/gencode.v39.annotation.gtf"
# opts$outdir = "results/multiome/bam_pass_qc_barcodes/Sample_5124-NM-2-hg38"
# opts$p = "auto"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Look at expression patterns in empty droplets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Read in empty droplet matrix
d = readRDS(opts$counts_empty)

## Prop
sorted_gene_counts = sort(rowSums(d), decreasing=TRUE)
summd = data.frame(gene = names(sorted_gene_counts), count = sorted_gene_counts)
summd$prop = summd$count/sum(summd$count)
summd$cumulative_prop = cumsum(summd$prop)
summd$rank_empty = seq(1:nrow(summd))

## Get nuclei counts for same genes
d_nuc = readRDS(opts$counts_nuclei)
d_nuc_ordered = d_nuc[row.names(summd),]
summd$count_nuc = rowSums(d_nuc_ordered)
summd$prop_nuc = summd$count_nuc/sum(summd$count_nuc)
summd$rank_nuc = rank(-1*summd$count_nuc)
summd$cumulative_prop_nuc = cumsum(summd$prop_nuc)
#cor(summd$rank_empty, summd$rank_nuc, method="spearman")

### Get gene annotation information
gtf = rtracklayer::import(opts$gtf)

## Determine ambient gene threshold
if (opts$p == "auto") {
  summd$cum_nuc_lt_cum_emp = as.numeric(summd$cumulative_prop_nuc > summd$cumulative_prop)
  summd$x = cumsum(summd$cum_nuc_lt_cum_emp)
  summd$ambient_indicator = summd$x == 0
  ambient_genes = row.names(summd)[summd$ambient_indicator==1]
  cat("Filtering out ", 100*(length(ambient_genes)/nrow(summd)), "%of genes\n")
} else {
  thresh = quantile(rowSums(d), p=(1-as.numeric(opts$p)))
  ambient_genes = row.names(d)[rowSums(d) > thresh]
  cat("Filtering out ", 100*(length(ambient_genes)/nrow(summd)), "%of genes\n")
}


### Output ambient gene regions to bed file
gtf_ambient = reduce(gtf[gtf$gene_name %in% ambient_genes & gtf$type=="gene",])
export.bed(gtf_ambient, con=file.path(opts$outdir,paste0("ambient_gene_regions_",opts$p,".bed")))


# Cumsum plot
summd_long = melt(summd[,c("rank_empty","cumulative_prop", "cumulative_prop_nuc")], id.vars="rank_empty")
summd_long$variable = factor(summd_long$variable, levels=c("cumulative_prop", "cumulative_prop_nuc"), labels=c("Empty", "Nuclei"))
png(file.path(opts$outdir,paste0("cumsum_ambient_genes_", opts$p, ".png")))
ggplot(summd_long[summd_long$rank_empty<=length(ambient_genes)*1.2,]) +
  geom_point(aes(x=rank_empty, y=value, color=variable), size=0.1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  geom_vline(xintercept=length(ambient_genes)+1,linetype = "longdash") +
  ggtitle(paste0("Excluding ", length(ambient_genes), " genes (", round(100*(length(ambient_genes)/nrow(summd)), digits=1),"%)")) +
  ylab("Cumulative proportion of GEX counts") + xlab("Gene rank in empty droplets")
dev.off()



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Look at expression patterns in nuclei
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ### Read in empty droplet matrix
# d_n = readRDS("results/multiome/counts_by_sample_gex/Sample_5124-NM-1-hg38/counts_nuclei.rds")
#
# # Prop
# sorted_gene_counts = sort(rowSums(d), decreasing=TRUE)
# summd = data.frame(gene = names(sorted_gene_counts), count = sorted_gene_counts)
# summd$prop = summd$count/sum(summd$count)
# summd$cummulative_prop = cumsum(summd$prop)
# summd$x = seq(1:nrow(summd))
#
# # Cumsum plot
# png('cumsum_nuclei_genes.png')
# ggplot(summd) +
#   geom_point(aes(x=x, y=cummulative_prop)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
# dev.off()
