library(yaml)
library(rjson)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     batch_design.json
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
batch_info = read.table(file="data/batch_design3.tsv", header=TRUE)

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
  demultiplex_list[[s]][["barcodes"]] = gsub("{SAMPLE}", s, "results/multiome/bam_pass_qc_barcodes/{SAMPLE}/pass_qc_barcodes_gex.txt", fixed=TRUE)
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
