options(stringsAsFactors=FALSE)
library(rtracklayer)
library(data.table)
library(optparse)
library(dplyr)


option_list <- list(
  make_option(
    c("--gwas_file"), type = "character", help = " "
  ),
  make_option(
    c("--outfile"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser, positional_arguments=TRUE)

annot_files = opts$args
opts = opts$options

### Testing
# opts = list()
# opts$gwas_file = "results/fgwas/gwas_test_file.txt.gz"
# opts$outfile = "results/fgwas/fgwas_input_test.txt.gz"
# annot_files = c("results/fgwas/annotations/hg19/beta_ATAC_300b.bed", "results/fgwas/annotations/hg19/acinar_ATAC_300b.bed")


### Get annotations
annot_labels = gsub(".bed", "", sapply(annot_files, function(x) {vals=unlist(strsplit(x, split="/")); return(vals[length(vals)])}))
annot_sets = lapply(annot_files, FUN=rtracklayer::import.bed)
names(annot_sets) <- annot_labels

### Get gwas stats
gwas_stats = data.table::fread(opts$gwas_file, header=TRUE)

### Fix missing SNPID
gwas_stats$SNPID[is.na(gwas_stats$SNPID)] <- paste0(gwas_stats$CHR,"_", gwas_stats$POS)[is.na(gwas_stats$SNPID)]

### Convert gwas stats to GRanges object
gwas_gr <- with(gwas_stats, GRanges(seqnames = Rle(CHR),
                ranges = IRanges(start = POS - 1, end = POS, names = SNPID)))

### Get annot overlap
getAnnotationCol = function(set) {
    hits <- findOverlaps(gwas_gr, set)
    hits <- names(gwas_gr)[queryHits(hits)] %>%
       unique()
    col = as.numeric(gwas_stats$SNPID %in% hits)
    return(col)
}
cols = lapply(annot_sets, getAnnotationCol)
colsDF = as.data.frame(do.call(cbind, cols))

### Add to gwas_data
out = data.frame(gwas_stats, colsDF)

### Remove duplicate positions;
# NOTE must do this because fgwas will throw an error
# about unsorted SNPs if there are two variants at same position
out = out[!duplicated(paste0(out$CHR,"_",out$POS)),]
fwrite(out, file=opts$outfile, sep=" ", compress="gzip")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ricardo's code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Read input file
# message("Reading feature BED file...")
# feature <- read.table(feature_file) %>%
#   select(V1:V3)
#
# feature <- with(
#   feature, GRanges(seqnames = Rle(V1),
#                    ranges = IRanges(start = V2, end = V3))
# )
#
# # Read all SNPs
# message("Reading summary stats...")
# snps_trim <- fread(summstats_file, header = T, data.table = F)
#
# # Get regions that match feature
# snps_trim_gr <- with(
#   snps_trim, GRanges(seqnames = Rle(CHR),
#                  ranges = IRanges(start = POS, end = POS + 1, names = SNPID))
# )
# hits <- findOverlaps(snps_trim_gr, feature)
# hits <- names(snps_trim_gr)[queryHits(hits)] %>%
#   unique()
#
#
# ## Format annotation for fGWAS
# message("Writing output...")
# snps_out <- snps_trim %>%
#   mutate(is_feature = as.numeric(SNPID %in% hits))
#
# to_change <- which(colnames(snps_out) == "is_feature")
# colnames(snps_out)[to_change] <- feature_name
#
# # Write output
# fwrite(snps_out, file = outfile, col.names = T, row.names = F, quote = F, sep = "\t")
# message("Done!")
