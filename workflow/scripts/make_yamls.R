library(yaml)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     demultiplex.yaml
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
batches = c("NM-1","NM-2","NM-3","NM-4","NM-5","NM-6","NM-7","NM-8")

multiome_samples = c("Sample_5124-NM_NM-1",
                      "Sample_5124-NM_NM-2",
                      "Sample_5124-NM_NM-3",
                      "Sample_5124-NM_NM-4",
                      "Sample_5124-NM_NM-5",
                      "Sample_5124-NM_NM-6",
                      "Sample_5124-NM_NM-7",
                      "Sample_5124-NM_NM-8")

fiveprime_samples = c("Sample_5125-NM-1_NM-1",
                          "Sample_5125-NM-2-8_NM-2",
                          "Sample_5125-NM-2-8_NM-3",
                          "Sample_5125-NM-2-8_NM-4",
                          "Sample_5125-NM-2-8_NM-5",
                          "Sample_5125-NM-2-8_NM-6",
                          "Sample_5125-NM-2-8_NM-7",
                          "Sample_5125-NM-2-8_NM-8")

demultiplex_list = list()
for (i in 1:length(multiome_samples)) {
  s = multiome_samples[i]
  b = batches[i]
  cat(s,"\n")
  demultiplex_list[[s]] = list()
  demultiplex_list[[s]][["batch"]] = b
  demultiplex_list[[s]][["bam_unfiltered"]] = gsub("{BATCH}", b, "results/multiome/nf_gex_results/prune/Sample_5124-NM-{BATCH}-hg38.before-dedup.bam", fixed=TRUE)
  demultiplex_list[[s]][["bam_auto"]] = gsub("{BATCH}", b, "results/multiome/bam_pass_qc_barcodes/Sample_5124-NM-{BATCH}-hg38/prune_barcodes_gex_exclude_ambient_genes_auto.bam", fixed=TRUE)
  demultiplex_list[[s]][["bam_p01"]] = gsub("{BATCH}", b, "results/multiome/bam_pass_qc_barcodes/Sample_5124-NM-{BATCH}-hg38/prune_barcodes_gex_exclude_ambient_genes_0.01.bam", fixed=TRUE)
  demultiplex_list[[s]][["barcodes"]] = gsub("{BATCH}", b, "results/multiome/bam_pass_qc_barcodes/Sample_5124-NM-{BATCH}-hg38/pass_qc_barcodes_gex.txt", fixed=TRUE)
}
for (i in 1:length(fiveprime_samples)) {
  s = fiveprime_samples[i]
  cat(s,"\n")
  demultiplex_list[[s]] = list()
  demultiplex_list[[s]][["batch"]] = batches[i]
  demultiplex_list[[s]][["bam_unfiltered"]] = gsub("{sample}", s, "results/fiveprime/starsolo/{sample}/starsolo.Aligned.sortedByCoord.out.bam", fixed=TRUE)
  demultiplex_list[[s]][["barcodes"]] = gsub("{sample}", s, "results/fiveprime/droplet_utils/{sample}/barcodes_nuclei.txt", fixed=TRUE)
}


write_yaml(list(samples=demultiplex_list), file="workflow/src/demultiplex_samples.yaml")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     fgwas.yaml
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clusters = c('acinar','alpha','beta','delta','ductal','endothelial','gamma','immune','stellate')

annot_list = list()
for (c in clusters) {
  cat(c,"\n")
  annot_list[[paste0(c, "_tCRE")]] = gsub("{cluster}", c, "results/scafe/scafe_by_cluster/{cluster}/annotate/{cluster}/bed/{cluster}.CRE.annot.bed", fixed=TRUE)
  annot_list[[paste0(c, "_ATAC")]] = gsub("{cluster}", c, "results/macs2/peakcalls/noNM3/{cluster}_summits_ext500_noblacklist.bed", fixed=TRUE)
}

write_yaml(list(fgwas_annotations=annot_list), file="workflow/src/fgwas_annotations.yaml")
