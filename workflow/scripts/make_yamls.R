library(yaml)
library(rjson)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     batch_design.json
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
batch_info = read.table(file="data/batch_design.tsv", header=TRUE)

split_by = function(x, by) {
  return(split(x, f=x[,by]))
}

convertBatchDFtoList = function(x) {
  x_by_batch = split_by(x, "Batch_name")
  x_by_batch_by_donor = lapply(x_by_batch, FUN=function(x) {y = split_by(x, by="Donor"); l = lapply(y, FUN=function(d) {as.list(d[,!names(d) %in% c("Batch_name","Donor")])}); return(list(donors=l))})
  #x_by_batch_by_donor_nest_list = lappy(x_by_batch_by_donor, FUN=function(x) {as.list(x[,!names(x) %in% c("Batch_name","Donor")])})
  return(x_by_batch_by_donor)
}

batch_list = convertBatchDFtoList(batch_info)
batch_json = toJSON(list(batches=batch_list), indent=4, method="C")
write(batch_json, file="data/batch_design.json")
write_yaml(list(batches=batch_list), file="data/batch_design.yaml")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     library_key.json
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library_key = read.table(file="data/library_key.tsv", header=TRUE)

convertLibraryDFtoJSON = function(x) {
  x_by_library = split(x, f=x$Master_SampleID)
  x_by_library_nested_list = lapply(x_by_library, as.list)
  x_json = toJSON(list(libraries=x_by_library_nested_list), indent=4, method="C")
  return(x_json)
}

library_json = convertLibraryDFtoJSON(library_key)
write(library_json, file="data/library_key.json")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     demultiplex.yaml
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
multiome_samples = c("Sample_5124-NM-1",
                      "Sample_5124-NM-2",
                      "Sample_5124-NM-3",
                      "Sample_5124-NM-4",
                      "Sample_5124-NM-5",
                      "Sample_5124-NM-6",
                      "Sample_5124-NM-7",
                      "Sample_5124-NM-8")

fiveprime_samples = c("Sample_5125-NM-1",
                          "Sample_5125-NM-2",
                          "Sample_5125-NM-3",
                          "Sample_5125-NM-4",
                          "Sample_5125-NM-5",
                          "Sample_5125-NM-6",
                          "Sample_5125-NM-7",
                          "Sample_5125-NM-8")

demultiplex_list = list()
for (i in 1:length(multiome_samples)) {
  s = multiome_samples[i]
  cat(s,"\n")
  demultiplex_list[[s]] = list()
  demultiplex_list[[s]][["barcodes"]] = gsub("{SAMPLE}", s, "results/multiome/droplet_utils/{SAMPLE}/barcodes_nuclei.txt", fixed=TRUE)
  demultiplex_list[[s]][["bam_unfiltered"]] = gsub("{SAMPLE}", s, "results/multiome/nf_gex_results/prune/{SAMPLE}-hg38.before-dedup.bam", fixed=TRUE)
  demultiplex_list[[s]][["bam_auto"]] = gsub("{SAMPLE}", s, "results/multiome/bam_pass_qc_barcodes/{SAMPLE}/prune_barcodes_gex_exclude_ambient_genes_auto.bam", fixed=TRUE)
  demultiplex_list[[s]][["bam_p01"]] = gsub("{SAMPLE}", s, "results/multiome/bam_pass_qc_barcodes/{SAMPLE}/prune_barcodes_gex_exclude_ambient_genes_0.01.bam", fixed=TRUE)
}
for (i in 1:length(fiveprime_samples)) {
  s = fiveprime_samples[i]
  cat(s,"\n")
  demultiplex_list[[s]] = list()
  demultiplex_list[[s]][["barcodes"]] = gsub("{SAMPLE}", s, "results/fiveprime/droplet_utils/{SAMPLE}/barcodes_nuclei.txt", fixed=TRUE)
  demultiplex_list[[s]][["bam_unfiltered"]] = gsub("{SAMPLE}", s, "results/fiveprime/starsolo/{SAMPLE}/starsolo.Aligned.sortedByCoord.out.bam", fixed=TRUE)
}
write_yaml(list(samples=demultiplex_list), file="workflow/src/demultiplex_samples.yaml")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     liger.yaml
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
multiome_samples = c("Sample_5124-NM-1",
                      "Sample_5124-NM-2",
                      "Sample_5124-NM-3",
                      "Sample_5124-NM-4",
                      "Sample_5124-NM-5",
                      "Sample_5124-NM-6",
                      "Sample_5124-NM-7",
                      "Sample_5124-NM-8")

fiveprime_samples = c("Sample_5125-NM-1",
                          "Sample_5125-NM-2",
                          "Sample_5125-NM-3",
                          "Sample_5125-NM-4",
                          "Sample_5125-NM-5",
                          "Sample_5125-NM-6",
                          "Sample_5125-NM-7",
                          "Sample_5125-NM-8")

threeprime_samples = c("Sample_HPAP036_gex",
                          "Sample_HPAP038_gex",
                          "Sample_HPAP039_gex",
                          "Sample_HPAP040_gex",
                          "Sample_HPAP044_gex",
                          "Sample_HPAP045_gex",
                          "Sample_HPAP055_gex",
                          "Sample_HPAP059-CVB4_gex",
                          "Sample_HPAP059-Cyto_gex",
                          "Sample_HPAP059-Mock_gex",
                          "Sample_ICRH122_gex",
                          "Sample_ICRH134-CVB4_gex",
                          "Sample_ICRH134-Cyto_gex",
                          "Sample_ICRH134-Mock_gex",
                          "Sample_ICRH135-CVB4_gex",
                          "Sample_ICRH135-Cyto_gex",
                          "Sample_ICRH135-Mock_gex")

#atac_samples = c()

counts_list = list()
for (i in 1:length(multiome_samples)) {
  s = multiome_samples[i]
  cat(s,"\n")
  counts_list[[s]] = gsub("{SAMPLE}", s, "results/multiome/decontx_round2/{SAMPLE}/counts_low_contamination_decontaminated.rds", fixed=TRUE)
}
for (i in 1:length(fiveprime_samples)) {
  s = fiveprime_samples[i]
  cat(s,"\n")
  counts_list[[s]] = gsub("{SAMPLE}", s, "results/fiveprime/decontx_round2/{SAMPLE}/counts_low_contamination_decontaminated.rds", fixed=TRUE)
}
for (i in 1:length(threeprime_samples)) {
  s = threeprime_samples[i]
  cat(s,"\n")
  counts_list[[s]] = gsub("{SAMPLE}", s, "results/threeprime/decontx_round2/{SAMPLE}/counts_low_contamination_decontaminated.rds", fixed=TRUE)
}

write_yaml(list(counts=counts_list), file="workflow/src/liger_samples.yaml")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     fgwas.yaml
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Annotation sets
annot_sets = list()

#NOTE: excluding endothelial and immune, too few cells
my_clusters = c('acinar','alpha','beta','delta','ductal','gamma','stellate', 'immune', 'endothelial')
annot_sets[["tCRE"]] = c(paste0(my_clusters, "_tCRE"), "wholeIslet_tCRE")
annot_sets[["ATAC_1kb"]] = c(paste0(my_clusters, "_ATAC_1kb"), "wholeIslet_ATAC_1kb")
annot_sets[["ATAC_300b"]] = c(paste0(my_clusters, "_ATAC_300b"), "wholeIslet_ATAC_300b")



ricardo_clusters = c("INS", "GCG", "KRT19", "PDGFRA", "PPY", "PRSS1","RGS5","SDS","SST", "TPSAB1", "VWF")
annot_sets[["Ricardo_clusters"]] = c(paste0(ricardo_clusters,"_summits.ext150.min2"),
                                    paste0(ricardo_clusters,"_summits.ext150.reproducible"),
                                    paste0(ricardo_clusters,"_summits.ext150.cluster_specific"),
                                    paste0(ricardo_clusters,"_summits.ext150.noblacklist"),
                                    "wholeIslet_Ricardo_clusters"
                                  )
## Annot source files
annot_list = list()
annot_list[["hg38"]] = list()
for (c in my_clusters) {
  cat(c,"\n")
  annot_list[["hg38"]][[paste0(c, "_tCRE")]] = gsub("{cluster}", c, "results/scafe_old/scafe_by_cluster/{cluster}/annotate/{cluster}/bed/{cluster}.CRE.annot.bed", fixed=TRUE)
  annot_list[["hg38"]][[paste0(c, "_ATAC_1kb")]] = gsub("{cluster}", c, "results/macs2/peakcalls/noNM3/{cluster}_summits_ext500_noblacklist.bed", fixed=TRUE)
  annot_list[["hg38"]][[paste0(c, "_ATAC_300b")]] = gsub("{cluster}", c, "results/macs2/peakcalls/noNM3/{cluster}_summits_ext150_noblacklist.bed", fixed=TRUE)
}
annot_list[["hg38"]][["wholeIslet_tCRE"]] = "results/fgwas/annotations/wholeIslet_tCRE.bed"
annot_list[["hg38"]][["wholeIslet_ATAC_300b"]] = "results/fgwas/annotations/wholeIslet_ATAC_300b.bed"
annot_list[["hg38"]][["wholeIslet_ATAC_1kb"]] = "results/fgwas/annotations/wholeIslet_ATAC_1kb.bed"
annot_list[["hg38"]][["wholeIslet_Ricardo_clusters"]] = "results/fgwas/annotations/wholeIslet_Ricardo_clusters.bed"


annot_list[["hg19"]] = list()
for (c in ricardo_clusters) {
  cat(c,"\n")
  annot_list[["hg19"]][[paste0(c, "_summits.ext150.min2")]] = gsub("{cluster}", c, "results/fgwas/annotations/ricardo_annots_hg19_reformatted/{cluster}_summits.ext150.min2.bed.bed3", fixed=TRUE)
  annot_list[["hg19"]][[paste0(c, "_summits.ext150.reproducible")]] = gsub("{cluster}", c, "results/fgwas/annotations/ricardo_annots_hg19_reformatted/{cluster}_summits.ext150.reproducible.bed3", fixed=TRUE)
  annot_list[["hg19"]][[paste0(c, "_summits.ext150.cluster_specific")]] = gsub("{cluster}", c, "results/fgwas/annotations/ricardo_annots_hg19_reformatted/{cluster}_summits.ext150.cluster_specific.bed.bed3", fixed=TRUE)
  annot_list[["hg19"]][[paste0(c, "_summits.ext150.noblacklist")]] = gsub("{cluster}", c, "results/fgwas/annotations/ricardo_annots_hg19_reformatted/{cluster}_summits.noblacklist.ext150.bed.bed3", fixed=TRUE)
  annot_list[["hg19"]][[paste0(c, "_peaks.narrowPeak")]] = gsub("{cluster}", c, "results/fgwas/annotations/ricardo_annots_hg19_reformatted/{cluster}_peaks.narrowPeak.bed3", fixed=TRUE)
}

write_yaml(list(fgwas_sets=annot_sets, fgwas_annotation_sources=annot_list), file="workflow/src/fgwas_annotations.yaml")
