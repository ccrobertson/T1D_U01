options(stringsAsFactors=FALSE)
library(optparse)
library(rtracklayer)
library(ggplot2)

option_list <- list(
  make_option(
    c("--counts_empty"), type = "character", help = "Input file name for storing empty droplets counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--gtf"), type = "character", help = "Genome annotation file used for generating count matrix."
  ),
  make_option(
    c("--p"), type = "numeric", help = "Filter out top p% of highest expressed genes in empty droplets."
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
# opts$gtf = "resources/hg38/gencode.v39.annotation.gtf"
# opts$outdir = "results/multiome/bam_pass_qc_barcodes/Sample_5124-NM-2-hg38"
# opts$p = 0.01
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Look at expression patterns in empty droplets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Read in empty droplet matrix
d = readRDS(opts$counts_empty)

### Define ambient genes
thresh = quantile(rowSums(d), p=(1-opts$p))
ambient_genes = row.names(d)[rowSums(d)>thresh]

### Get gene annotation information
gtf = rtracklayer::import(opts$gtf)

### Output ambient gene regions to bed file
gtf_ambient = reduce(gtf[gtf$gene_name %in% ambient_genes & gtf$type=="gene",])
export.bed(gtf_ambient, con=file.path(opts$outdir,paste0("ambient_gene_regions_",opts$p,".bed")))

# # Prop
sorted_gene_counts = sort(rowSums(d), decreasing=TRUE)
summd = data.frame(gene = names(sorted_gene_counts), count = sorted_gene_counts)
summd$prop = summd$count/sum(summd$count)
summd$cummulative_prop = cumsum(summd$prop)
summd$x = seq(1:nrow(summd))

# Cumsum plot
png(file.path(opts$outdir,paste0("cumsum_ambient_genes_", opts$p, ".png")))
ggplot(summd) +
  geom_point(aes(x=x, y=cummulative_prop)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  geom_vline(xintercept=length(ambient_genes)+1,linetype = "longdash")
dev.off()


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Look at expression patterns in nuclei
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ### Read in empty droplet matrix
# d = readRDS("results/multiome/counts_by_sample_gex/Sample_5124-NM-1-hg38/counts_nuclei.rds")
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
