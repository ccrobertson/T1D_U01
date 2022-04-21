
Starting a repository for storing code and documentation related to the T1D U01
collaborative project: Parker, Chen, and Collins labs

### Data location
* Multiome: /lab/data/seqcore/5124-NM
* 5' RNA-seq: /lab/data/seqcore/5125-NM
* Genotyping: /lab/data/nih/2022_02_16_T1D_genotypes

* From Ricardo's paper
* * Genotyping: /lab/data/nih/2021_07_02_T1D_genotypes/chip_genotypes
* * snRNA-seq libraries:
* * snATAC-seq libraries:
* * code for cvb4 gene quantification: Ricardo_code/notebooks/cvb4_counts.Rmd

### Environmental variables
```bash
export WORK=/lab/work/ccrober/T1D_U01
export freeze=${WORK}/results/freeze1
# export data=${WORK}/data
# export resources=${WORK}/resources
```

### Check md5
```bash
checkLibMD5sum () {
  local dir=$1
  local lib=$2
  cd ${dir}/${lib}
  md5sum -c ${lib}.md5 > ${lib}.md5check
}

cd ${freeze}
seqcore=/lab/data/seqcore
checkLibMD5sum ${seqcore} "5124-NM"
checkLibMD5sum ${seqcore} "5125-NM"
```


### Create symlinks to fastq files
```bash
createSymLinksATAC () {
  local target_dir=$1
  local link_dir=$2
  local sample=$3
  ln -s --force ${target_dir}/${sample}/*R1_001.fastq.gz ${link_dir}/${sample}_R1_001.fastq.gz
  ln -s --force ${target_dir}/${sample}/*R2_001.fastq.gz ${link_dir}/${sample}_R2_001.fastq.gz
  ln -s --force ${target_dir}/${sample}/*R3_001.fastq.gz ${link_dir}/${sample}_R3_001.fastq.gz
}

createSymLinksGEX () {
  local target_dir=$1
  local link_dir=$2
  local sample=$3
  ln -s --force ${target_dir}/${sample}/*R1_001.fastq.gz ${link_dir}/${sample}_R1_001.fastq.gz
  ln -s --force ${target_dir}/${sample}/*R2_001.fastq.gz ${link_dir}/${sample}_R2_001.fastq.gz
}


#symlinks for ATAC
for i in {1..8}; do
  for lane in a b c d; do
    sample=Sample_5124-NM-${i}-ATAC_${lane}
    echo ${sample}
    createSymLinksATAC "/lab/data/seqcore/5124-NM/ATAC" "${data}/fastq_atac" ${sample}
  done
done

#symlinks for GEX
for i in {1..8}; do
  sample=Sample_5124-NM-${i}-3GEX
  echo ${sample}
  createSymLinksGEX "/lab/data/seqcore/5124-NM/GEX/fastqs_5124-NM" "${data}/fastq_gex" ${sample}
done
```


### Create small test samples
```bash
#One million reads
cd ${data}/fastq_atac
zcat Sample_5124-NM-1-ATAC_a_R1_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R1_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R2_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R2_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R3_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R3_001.fastq.gz

#One thousand reads
cd ${data}/fastq_atac
zcat Sample_5124-NM-1-ATAC_a_R1_001.fastq.gz | head -n 4000 | gzip > Sample_test_xs_R1_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R2_001.fastq.gz | head -n 4000 | gzip > Sample_test_xs_R2_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R3_001.fastq.gz | head -n 4000 | gzip > Sample_test_xs_R3_001.fastq.gz

#One million reads
cd ${data}/fastq_gex
zcat Sample_5124-NM-1-3GEX_R1_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R1_001.fastq.gz
zcat Sample_5124-NM-1-3GEX_R2_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R2_001.fastq.gz

#One thousand reads
cd ${data}/fastq_gex
zcat Sample_5124-NM-1-3GEX_R1_001.fastq.gz | head -n 4000 | gzip > Sample_test_xs_R1_001.fastq.gz
zcat Sample_5124-NM-1-3GEX_R2_001.fastq.gz | head -n 4000 | gzip > Sample_test_xs_R2_001.fastq.gz
```


### Prepare nextflow input jsons
```bash
#cut -d"," -f 2 /lab/data/seqcore/5124-NM/ATAC/DemuxStats_5124-NM.csv | awk 'NR>1' | sed 's/"//g'
cut -d"," -f 2 /lab/data/seqcore/5124-NM/GEX/DemuxStats_5124-NM.csv | awk 'NR>1' | sed 's/"//g;s/-3GEX//g' > ${freeze}/multiome_sample_list.tsv

python ${scripts}/build_multiome_json.py --sample-file ${freeze}/multiome_sample_list.tsv --fastq-dir ${data}/fastq_atac --modality ATAC > ${freeze}/nf_atac_config.json
python ${scripts}/build_multiome_json.py --sample-file ${freeze}/multiome_sample_list.tsv --fastq-dir ${data}/fastq_gex --modality GEX > ${freeze}/nf_gex_config.json
```



### Process multiome ATAC data
For each readgroup, three fastq files are required:
* the first and second insert reads ('1' and '2')
* and the read with the nuclear barcode ('index')
```bash
module load golang/1.13.5 singularity/3.5.1
cd ${freeze}/nf_atac_results
nohup nextflow run -resume -params-file ${freeze}/nf_atac_config_test.json --results ${freeze}/nf_atac_results pipelines/snATACseq-NextFlow/main.nf
nohup nextflow run -resume -params-file ${freeze}/nf_atac_config.json --results ${freeze}/nf_atac_results pipelines/snATACseq-NextFlow/main.nf
```


### Process multiome GEX data
```bash
module load golang/1.13.5 singularity/3.5.1
cd ${freeze}/nf_gex_results
nohup nextflow run -resume -params-file ${freeze}/nf_gex_config_test.json --chemistry multiome --results ${freeze}/nf_gex_results pipelines/snRNAseq-NextFlow/main.nf &
nohup nextflow run -resume -params-file ${freeze}/nf_gex_config.json --chemistry multiome --results ${freeze}/nf_gex_results pipelines/snRNAseq-NextFlow/main.nf &
```


### Checking to see if any reads align to CVB4 viral genome
```bash
cd ${freeze}/nf_gex_results/starsolo/Sample_5124-NM-1-hg38
samtools index Sample_5124-NM-1-hg38.Aligned.sortedByCoord.out.bam
samtools view Sample_5124-NM-1-hg38.Aligned.sortedByCoord.out.bam AF311939.1
```


### Downstream multiome pipeline
```bash
mamba activate snakemake
cd ${WORK}
bash commands.sh

snakemake all -n --jobname "islets.{jobid}" --jobs 100 \
		--keep-going \
		--rerun-incomplete \
		--snakefile workflow/src/multiome.smk \
		--configfile workflow/src/multiome.yaml \
		--use-conda \
		--printshellcmds \
		--cluster-config workflow/envs/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus}" \
		> logs/snakemake.log 2>&1 &
```

## Jobs to keep
date            rule                                                     wildcards jobid   slurmid
31  Thu Apr 21 11:38:22 2022        dropkick       method=otsu, min_genes=30, sample=Sample_5124-NM-2-hg38    45  13710107
1   Thu Apr 21 11:38:20 2022        dropkick       method=otsu, min_genes=50, sample=Sample_5124-NM-2-hg38    46  13710077 Finished
6   Thu Apr 21 11:38:21 2022        dropkick  method=multiotsu, min_genes=30, sample=Sample_5124-NM-2-hg38    47  13710082
11  Thu Apr 21 11:38:21 2022        dropkick  method=multiotsu, min_genes=50, sample=Sample_5124-NM-2-hg38    48  13710087
2   Thu Apr 21 11:38:21 2022        dropkick  method=multiotsu, min_genes=50, sample=Sample_5124-NM-5-hg38    60  13710078 RUNNING


### WHAT DO I WANT TO DO?
1. create "nuclei" and "empty" GEX count matrices
  input: count matrix, nuclei barcode list, empty barcode list, doublet barcode list
  output: counts_nuclei.rds, counts_empty.rds
2. correct for ambient contamination in GEX
  input: counts_nuclei.rds, counts_empty.rds
  output: counts_nuclei_decontaminated.rds (result$decontXcounts)
3. barcode clustering on GEX (LIGER - http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html)
  input: counts_nuclei_decontaminated.rds  (and maybe also run on counts_nuclei.rds  for comparison), will have one count matrix per sample, which we will read into a list and integrate
  output: integrated umap across samples, clusters (barcode lists for each cluster)
  NOTE: once we have 5' GEX bam files, will want to do cross modality integration using http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html

For each cell type (pseudobulk analyses):
  4. generate atac-seq bam files using cluster barcode lists
  5. generate atac-seq peak calls and signal tracks
  6. generate pseudobulk count matrices
    - ATAC-seq: sample x peak
    - 3' GEX: sample x gene
    - 5' GEX: sample x tCRE (can further break into TSS proximal and distal)
  7. run paired differental analysis by condition/treatment
    -  Feature ~ subject + CVB4 + cytokine
  8. run caQTL, eQTL, tCRE analysis


### NEXT STEPS - Integrating across cell types and/or modalities
- Linking ATAC peaks to genes (single nucleus analysis) - Do nuclei with higher expression of gene X have more accessibility at peak Y? Existing tools?
- Linking peak/gene pairs to tCREs (pseudobulk analysis across individuals and/or cell types?)



### Process 5' scRNA-seq
```bash
scafe.workflow.sc.solo --help
```


### Genotype imputation
```bash
snakemake all --jobname "impute.{jobid}" --jobs 100 \
		--keep-going \
		--rerun-incomplete \
		--snakefile workflow/src/imputation.smk \
		--configfile workflow/src/imputation.yaml \
		--printshellcmds \
		--cluster-config workflow/envs/cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus}" \
		> logs/snakemake_imputation.log 2>&1 &
```
