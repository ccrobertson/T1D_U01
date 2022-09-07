options(stringsAsFactors=FALSE)
library(optparse)
library(rtracklayer)
library(plyranges)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script evaluates per-nucleus qc metrics from multiome GEX and ATAC
# and filters nuclei based on joint-modality thresholds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
option_list <- list(
  make_option(
    c("--summits"), type = "character", help = "ATAC-seq MACS2 narrowPeak file"
  ),
  make_option(
    c("--gtf"), type = "character", help = "Genome annotation file for general promoters."
  ),
  make_option(
    c("--outdir"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# ### Testing
# opts = list()
# opts$summits = "results/multiome/features/Sample_5124-NM-1-hg38_summits.bed"
# opts$summits = "results/multiome/macs2/Sample_5124-NM-1-hg38_summits.bed"
# opts$gtf = "resources/hg38/gencode.v39.annotation.gtf"
# opts$outdir = "results/multiome/features"


### Read in annotation information
gtf = rtracklayer::import(opts$gtf)
genes = gtf[gtf$type=="gene",]
transcripts = gtf[gtf$type=="transcript",]
promoters = promoters(genes) #upstream 2000, downstream 200
genes_promoters = union(genes, promoters)


### Define ATAC region subsets
summits = rtracklayer::import(opts$summits, format="BED")
summits = resize(summits,  fix="center", ignore.strand=TRUE, width=1)
features = resize(summits, fix="center", ignore.strand=TRUE, width=1000)
features_promoters = subsetByOverlaps(features, promoters)
features_genic = subsetByOverlaps(features, genes_promoters)
features_intergenic = subsetByOverlaps(features, genes_promoters, invert=TRUE)

### Save regions
export.bed(features, con=file.path(opts$outdir,"atac_summits_1kb_all.bed"))
export.bed(features_promoters, con=file.path(opts$outdir,"atac_summits_1kb_promoters.bed"))
export.bed(features_genic, con=file.path(opts$outdir,"atac_summits_1kb_genes.bed"))
export.bed(features_intergenic, con=file.path(opts$outdir,"atac_summits_1kb_intergenic.bed"))


### Define reference points for aggregate plots
tss = resize(genes, width=1, fix="start")
#exclude tss within 1kb of another tss
tss_1kb = resize(tss, width=2000, fix="center", ignore.strand=TRUE)
tss_overlaps = countOverlaps(tss, tss_1kb)
tss_keep = tss[tss_overlaps==1]
export.bed(reduce(tss), con=file.path(opts$outdir,"tss_all_genes.bed"))
export.bed(reduce(tss_keep), con=file.path(opts$outdir,"tss_nonoverlapping_genes.bed"))

#tss_open = subsetByOverlaps(tss, features)
#export.bed(reduce(tss_open), con=file.path(opts$outdir,"tss_open.bed"))

summits_intergenic = subsetByOverlaps(summits, features_intergenic)
export.bed(reduce(summits_intergenic), con=file.path(opts$outdir,"atac_summits_intergenic.bed"))

#tss_2 = rtracklayer::import("resources/hg38/hg38.refGene.tss.bed.gz")
#countOverlaps(tss, tss_2)

# cres = rtracklayer::import("results/fiveprime/5125-NM-1-rerun/scafe/annotate/Sample_5125-NM-1-5GEX/bed/Sample_5125-NM-1-5GEX.CRE.coord.bed.gz")
# cres_intergenic = subsetByOverlaps(cres, features_intergenic)
# export.bed(cres_intergenic, con=file.path(opts$outdir,"cres_overlapping_intergenic_atac_features.bed"))
#
#
# top_cres_intergenic = cres_intergenic[cres_intergenic$score>1000]
# top_cres_intergenic_region = resize(top_cres_intergenic, fix="center", ignore.strand=TRUE, width=5e5)
# export.bed(reduce(top_cres_intergenic), con=file.path(opts$outdir,"top_cres_overlapping_intergenic_atac_features.bed"))
# export.bed(reduce(top_cres_intergenic_region), con=file.path(opts$outdir,"top_cres_overlapping_intergenic_atac_features_loci.bed"))
