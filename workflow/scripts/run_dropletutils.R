options(stringsAsFactors=FALSE)
library(optparse)
library(DropletUtils)

#https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html

option_list <- list(
  make_option(
    c("--input_10x_dir"), type = "character", help = "Directory of starsolo output containing matrix.mtx and genes.tsv files"
  ),
  make_option(
    c("--outdir"), type = "character", help = " "
  )
)


option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# # ### Testing
# opts = list()
# opts$input_10x_dir = "results/fiveprime/starsolo/Sample_5125-NM-1_NM-1_keep/starsolo.Solo.out/GeneFull_ExonOverIntron/raw"
# opts$outdir = "results/fiveprime/droplet_utils/Sample_5125-NM-1_NM-1"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sce <- read10xCounts(opts$input_10x_dir, col.names=TRUE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Barcode rank plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my.counts = counts(sce)
br.out <- barcodeRanks(my.counts)

knee = metadata(br.out)$knee
inflection = metadata(br.out)$inflection

png(file.path(opts$outdir, "barcode_rank.png"))
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=knee, col="dodgerblue", lty=2)
abline(h=inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
    legend=c(
      paste0("knee (", knee,")"),
      paste0("inflection (", inflection,")")
    ))
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preliminary filtering for empty droplets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
barcodes_empty = colnames(my.counts)[br.out$total < 0.5*inflection]
write.table(barcodes_empty, file=file.path(opts$outdir, "barcodes_empty.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Empty droplets vs nuclei
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.table(colnames(my.counts)[br.out$total > knee], file=file.path(opts$outdir, "barcodes_above_knee.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(colnames(my.counts)[br.out$total > inflection], file=file.path(opts$outdir, "barcodes_above_inflection.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

set.seed(100)
e.out <- emptyDrops(my.counts)

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)

table(Limited=e.out$Limited, Significant=is.cell)

png(file.path(opts$outdir, "emptydrops_diagnostic.png"))
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()

saveRDS(e.out, file=file.path(opts$outdir, "emptydrops.rds"))

barcodes_nuclei = colnames(my.counts)[!is.na(is.cell) & is.cell==TRUE]
write.table(barcodes_nuclei, file=file.path(opts$outdir, "barcodes_nuclei.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
