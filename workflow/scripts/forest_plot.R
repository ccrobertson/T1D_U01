library(ggplot2)
library(optparse)


option_list <- list(
  make_option(
    c("--input_file"), type = "character", help = " "
  ),
  make_option(
    c("--plotfile"), type = "character", help = " "
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser, positional_arguments=FALSE)

### Testing
# opts = list()
# opts$input_file = "results/fgwas/fgwas_run/T1D-hg19-Ricardo_clusters.txt"
# opts$plotfile = "results/fgwas/fgwas_run/T1D-hg19-Ricardo_clusters.png"

d = read.table(opts$input_file)
names(d) = c("Tissue", "lower", "estimate", "upper")
d$index = 1:nrow(d)

# Deal with cases where CI is <X value
d[grep("<",d$lower), c("lower","estimate","upper")] <- c(0,0,0)

d$lower = as.numeric(d$lower)
d$estimate = as.numeric(d$estimate)
d$upper = as.numeric(d$upper)

png(opts$plotfile)
ggplot(data=d, aes(y=index, x=estimate, xmin=lower, xmax=upper)) +
  geom_point() +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(d), labels=d$Tissue) +
  labs(title=' ', x='Enrichment Size', y = 'Study') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()
dev.off()
