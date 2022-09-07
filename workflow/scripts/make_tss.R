# This script reads in the refGene annotation file downloaded from
# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
# and creates a TSS annotation file of 5' ends of all protein coding genes
#
# Code is adapted from https://github.com/porchard/ataqv/blob/master/data/tss/make_tss.R
#
library(dplyr)

### Read in refGene file
d = read.table("hg38.refGene.txt.gz")
cols = c("bin", "name","chrom","strand","txStart", "txEnd", "cdsStart", "cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
names(d) <- cols
d$score = 0
d2 <- d %>%
  dplyr::select(name, chrom, txStart, txEnd, strand, name2, score) %>%
  dplyr::rename(transcript = name, gene = name2)


### Define TSS
# if on - strand, start and end must be flipped
tss <- d2 %>%
  dplyr::mutate(tss_start = ifelse(strand == '+', txStart, txEnd - 1),
                tss_end = tss_start + 1) %>%
  dplyr::select(-txStart, -txEnd) %>%
  unique()


### Only keep protein-coding genes
tss <- tss[grep('^NM_', tss$transcript),]


### Eliminate redundant TSS, so that at a given position and strand, there is no more than one TSS
tss <- unique(tss[,c('chrom', 'tss_start', 'tss_end', 'gene', 'score', 'strand')])
tss$id <- with(tss, paste(chrom, tss_start, strand, sep = ':'))
number_unique <- length(unique(tss$id))
number_tss_at_position <- table(tss$id)
redundant <- names(number_tss_at_position[number_tss_at_position>1])
tss$keep <- T
for(i in redundant) {
  tss$keep[tss$id==i] <- c(T, rep(F, sum(tss$id==i)-1))
}
stopifnot(sum(tss$keep)==number_unique) # sanity check
tss <- tss[tss$keep,c('chrom', 'tss_start', 'tss_end', 'gene', 'score', 'strand')]


### Sort and write to file
tss <- tss[order(tss$chrom, tss$tss_start, tss$strand),]
write.table(x = tss, file="hg38.refGene.tss.bed", quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
