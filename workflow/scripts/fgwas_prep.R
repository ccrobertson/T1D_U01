options(stringsAsFactors=FALSE)
library(rtracklayer)
library(data.table)
library(optparse)


option_list <- list(
  make_option(
    c("--gwas_file"), type = "character", help = " "
  ),
  make_option(
    c("--annot_dir"), type = "character", help = " "
  ),
  make_option(
    c("--outfile"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)


# ### Testing
# opts = list()
# opts$gwas_file = "results/fgwas/gwas_test_file.txt.gz"
# opts$annot_dir =  "results/fgwas/annotations/whole_islet-hg19"
# opts$outfile = "results/fgwas/fgwas_input.txt.gz"


# ### Get annotations
annot_files = scan(file.path(opts$annot_dir, "annotations.txt"), what="character")
annot_sets = lapply(file.path(opts$annot_dir,annot_files), FUN=rtracklayer::import)
names(annot_sets) <- sub(pattern="([\\w]+).bed", "\\1", x=annot_files, perl=TRUE)


### Get gwas stats
gwas_stats = data.table::fread(opts$gwas_file, header=TRUE)

### Convert gwas stats to GRanges object
gwas_gr <- with(gwas_stats, GRanges(seqnames = Rle(CHR),
                ranges = IRanges(start = POS - 1, end = POS, names = SNPID)))

### Get annot overlap
getAnnotationCol = function(set) {
    hits <- as.numeric(!is.na(findOverlaps(gwas_gr, set, select="first")))
}
cols = lapply(annot_sets, getAnnotationCol)
colsDF = as.data.frame(do.call(cbind, cols))

### Add to gwas_data
out = data.frame(gwas_stats, colsDF)
out$SNPID[is.na(out$SNPID)] <- "."

### Remove duplicate positions
out_uniq = out[!duplicated(paste0(out$CHR,":",out$P)),]
fwrite(out_uniq, file=opts$outfile, sep=" ", compress="gzip")
