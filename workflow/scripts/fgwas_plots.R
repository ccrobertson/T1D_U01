library(ggplot2)

files = c(
"results/fgwas/fgwas_run/T1D_Chiou-hg38-ATAC_1kb_cond.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-ATAC_1kb.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-ATAC_300b_cond.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-ATAC_300b.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-Ricardo_clusters_cond.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-Ricardo_clusters.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-tCRE_cond.txt",
"results/fgwas/fgwas_run/T1D_Chiou-hg38-tCRE.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-ATAC_1kb_cond.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-ATAC_1kb.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-ATAC_300b_cond.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-ATAC_300b.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-Ricardo_clusters_cond.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-Ricardo_clusters.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-tCRE_cond.txt",
"results/fgwas/fgwas_run/T1D_Crouch-hg38-tCRE.txt"
)


processFile = function(x) {
  df = read.table(x, header=FALSE)
  names(df) <- c("annotation", "lower", "estimate", "upper")
  # Deal with cases where CI is <X value
  df[grep("<",df$lower), c("lower","estimate","upper")] <- c(NA,NA,NA)
  df[grep(">",df$upper), c("lower","estimate","upper")] <- c(NA,NA,NA)
  df[grep("fail",df$lower), c("lower","estimate","upper")] <- c(NA,NA,NA)
  df[grep("fail",df$upper), c("lower","estimate","upper")] <- c(NA,NA,NA)
  # Make sure vals are numeric
  df$lower = as.numeric(df$lower)
  df$estimate = as.numeric(df$estimate)
  df$upper = as.numeric(df$upper)
  df$source = x
  return(df)
}
d = do.call("rbind", lapply(files, processFile))


d$cluster = sapply(d$annotation, function(x) strsplit(x, split="_")[[1]][1])
d$modality = sapply(d$annotation, function(x) strsplit(x, split="_")[[1]][2])
d$atac_type = sapply(d$annotation, function(x) strsplit(x, split="_")[[1]][3])
d$atac_type[d$atac_type%in%c("clusters", "noMHC")] <- NA
d$gwas = NA
d$gwas[grep("Chiou", d$source)] <- "T1D_Chiou"
d$gwas[grep("Crouch", d$source)] <- "T1D_Crouch"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Based on my annotations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my_clusters = rev(c('beta','alpha','gamma','delta','acinar','ductal','stellate', 'immune', 'endothelial'))
my_cluster_labels = rev(c('Beta','Alpha','Gamma','Delta','Acinar','Ductal','Stellate', 'Immune', 'Endothelial'))
my_cluster_colors = rev(c(
  Beta = rgb(255,87,51, maxColorValue=255),
  Alpha = rgb(255,191,0, maxColorValue=255),
  Gamma = rgb(242,140,40, maxColorValue=255),
  Delta = rgb(255,172,28, maxColorValue=255),
  Acinar = rgb(128,222,234, maxColorValue=255),
  Ductal = rgb(121,85,72, maxColorValue=255),
  Stellate = rgb(76,175,80, maxColorValue=255),
  Immune = rgb(171,71,188, maxColorValue=255),
  Endothelial = rgb(144,164,174, maxColorValue=255)
))
# my_clusters = c('alpha','beta','gamma','delta','acinar','ductal','stellate', 'immune', 'endothelial')
# my_cluster_labels = c('Alpha','Beta','Gamma','Delta','Acinar','Ductal','Stellate', 'Immune', 'Endothelial')
# my_cluster_colors = c(
#   Alpha = rgb(255,191,0, maxColorValue=255),
#   Beta = rgb(255,87,51, maxColorValue=255),
#   Gamma = rgb(242,140,40, maxColorValue=255),
#   Delta = rgb(255,172,28, maxColorValue=255),
#   Acinar = rgb(128,222,234, maxColorValue=255),
#   Ductal = rgb(121,85,72, maxColorValue=255),
#   Stellate = rgb(76,175,80, maxColorValue=255),
#   Immune = rgb(171,71,188, maxColorValue=255),
#   Endothelial = rgb(144,164,174, maxColorValue=255)
# )


d1 = d[d$cluster %in% my_clusters,]
d1$clusterF = factor(d1$cluster, levels=my_clusters, labels=my_cluster_labels)
d1$clusterN = as.numeric(d1$clusterF)
d1$category = paste0(d1$modality, "_", d1$atac_type)
d1_marginal = d1[!grepl("_cond.", d1$source),]
d1_conditional = d1[grepl("_cond.", d1$source),]

## Restricting to Crouch GWAS
plotdat = d1_marginal[d1_marginal$gwas=="T1D_Crouch" & (d1_marginal$atac_type=="1kb" | is.na(d1_marginal$atac_type)),]
plotdat$Modality = factor(plotdat$modality, levels=c("ATAC", "tCRE"), labels=c("ATAC", "5' RNA"))

pdf("fgwas_plot_ASHG22_v2.pdf", width=12, height=6)
ggplot(data=plotdat, aes(x=estimate, y=clusterF, xmin=lower, xmax=upper, color=clusterF)) +
  geom_point(size=5) +
  scale_color_manual(values = my_cluster_colors, name="") +
  geom_errorbarh(height=.3, size=2) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_bw() +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=20)) +
  facet_grid(.~Modality) +
  xlab("fGWAS ln(Enrichment)") + ylab("")
dev.off()

pdf("fgwas_plot_ASHG22_v3.pdf", width=12, height=6)
ggplot(data=plotdat, aes(x=estimate, y=clusterF, xmin=lower, xmax=upper, color=Modality)) +
  geom_point(aes(shape=Modality), size=5) +
  #scale_color_manual(values = my_cluster_colors, name="") +
  geom_errorbarh(height=.3, size=2) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  #facet_grid(clusterF~.) +
  #ggtitle("Marginal enrichments") +
  guides(color = guide_legend(override.aes = list(size = 10, alpha=1))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines")) +
  xlab("fGWAS ln(Enrichment)")  + ylab("")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Based on Ricardo's annotations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ricardo_clusters = rev(c("INS", "GCG", "KRT19", "PDGFRA", "PPY", "PRSS1","RGS5","SDS","SST", "TPSAB1", "VWF"))
ricardo_cluster_labels = rev(c("Beta(INS)", "Alpha(GCG)", "Ductal(KRT19)", "PDGFRA", "Gamma(PPY)", "Acinar(PRSS1)","Stellate(RGS5)","Immune(SDS)","Delta(SST)", "TPSAB1", "Endothelial(VWF)"))
d2 = d[d$cluster %in% ricardo_clusters,]
d2$clusterF = factor(d2$cluster, levels=ricardo_clusters, labels=ricardo_cluster_labels)
d2$clusterN = as.numeric(d2$clusterF)
d2$category = paste0(d2$modality, "_", d2$atac_type)
d2_marginal = d2[!grepl("_cond.", d2$source),]
d2_conditional = d2[grepl("_cond.", d2$source),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conditional analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot(data=d1_conditional, aes(x=estimate, y=clusterF, xmin=lower, xmax=upper)) +
  geom_point() +
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_bw() +
  facet_grid(gwas~category) +
  ggtitle("Conditional enrichments")


ggplot(data=d2_marginal, aes(x=estimate, y=clusterF, xmin=lower, xmax=upper)) +
  geom_point() +
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_bw() +
  facet_grid(gwas~modality) +
  ggtitle("Marginal enrichments")

ggplot(data=d2_conditional, aes(x=estimate, y=clusterF, xmin=lower, xmax=upper)) +
  geom_point() +
  geom_errorbarh(height=.1) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_bw() +
  facet_grid(gwas~modality) +
  ggtitle("Conditional enrichments")
