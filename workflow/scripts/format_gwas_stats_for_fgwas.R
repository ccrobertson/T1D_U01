options(stringsAsFactors=FALSE)
library(data.table)
library(optparse)


option_list <- list(
  make_option(
    c("--gwas_file"), type = "character", help = "Harmonized gwas catalog file"
  ),
  make_option(
    c("--outfile"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser, positional_arguments=FALSE)


# ### Testing
# opts = list()
# opts$gwas_file = "results/fgwas/gwas_stats/T1D_Crouch-hg38/GCST90013791.h.tsv.gz"
# opts$outfile = "results/fgwas/gwas_stats/T1D_Crouch-hg38/GCST90013791_for_fgwas.tsv.gz"


### Read file
cat("Reading file\n")
d = fread(opts$gwas_file)

### Rename variables
cat("Renaming variables\n")
d$SNPID = d$hm_variant_id
d$CHR = paste0("chr", d$hm_chrom)
d$POS = d$hm_pos
d$F = round(d$hm_effect_allele_frequency, digits=5)
d$Z = round(d$beta/d$standard_error, digits=5)
d$N = 181214
d$SE = round(d$standard_error, digits=5)

### Filter for MAF
cat("Filtering for MAF\n")
d = d[d$F>0.005 & d$F<0.995,]

### Extract only necessary columns
cat("Extracting columns\n")
d = d[,c("SNPID","CHR","POS","F","Z","N","SE","hm_chrom")]

### Sort
cat("Sort by chr then by pos\n")
d = d[order(d$hm_chrom, d$POS),]

### Remove any missing values
cat("Remove rows with missing values\n")
d = d[complete.cases(d),]

### Save
cat("Save results to file\n")
fwrite(d[,c("SNPID","CHR","POS","F","Z","N","SE")], file=opts$outfile, sep="\t", compress = "gzip")
