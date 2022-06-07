
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

* Code from Arushi: /lab/work/arushiv/muscle-sn/analyses_hg3/cluster_sample_count_matrices/scripts/
* * main.nf (or sample.nf) - generating cell type specific ATAC-seq peak calls
* * quantify.nf - generate barcode x feature count matrix for ATAC-seq
* *

### Environmental variables
```bash
export WORK=/lab/work/ccrober/T1D_U01
export freeze=${WORK}/results/freeze1
# export data=${WORK}/data
# export resources=${WORK}/resources
```


### Batch data
```bash
awk -v OFS="\t" 'NR==1 {print "Libname", $0} NR>1 {print "Sample_5124-"$1"-hg38",$0}' Sample_islet_list_for_multiomics_Batches_long_format.tsv > Sample_islet_list_for_multiomics_Batches_long_format_with_libname.tsv

Rscript workflow/scripts/build_batch_json.R --batchfile Sample_islet_list_for_multiomics_Batches_long_format_with_libname.tsv
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


#One thousand droplets (i.e., barcodes)
barcodes_gex=${freeze}/barcodes_gex_1000.txt
bam_gex=${freeze}/nf_gex_results/prune/Sample_5124-NM-2-hg38.before-dedup.bam
bam_gex_1000=${freeze}/nf_gex_results/prune/Sample_test2-hg38.before-dedup.bam
awk 'NR>1' ${freeze}/cross_modality_qc/Sample_5124-NM-2-hg38/barcodes_nuclei.txt | head -n 100  | cut -d' ' -f 1 > ${barcodes_gex}
singularity exec workflow/envs/general.simg samtools view -h -b -@ 10 -D CB:${barcodes_gex} ${bam_gex} > ${bam_gex_1000}

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


### Split barcodes by origin individual
# NOTE:
# NM-1, NM-2, NM-3, and NM-4 pooled 5 individuals
# NM-5, NM-6, NM-7, and NM-8 pooled 6 individuals
#
# Number of individual in the pool may influence input parameter k?
```bash
## get common variants files
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_" -O common_variants_grch38.vcf && rm -rf /tmp/cookies.txt
## fix contig name
awk '$1 ~ /^#/ {print $0} $1 !~ /^#/ {print "chr"$0}' resources/common_variants_grch38.vcf > resources/common_variants_grch38_fixed.vcf

## run on test bam with only 1000 barcodes
sbatch --output=souporcell-1.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-1-hg38 5
#sbatch workflow/scripts/run_souporcell.slurm Sample_5124-NM-2-hg38 5
sbatch --output=souporcell-3.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-3-hg38 5
sbatch --output=souporcell-4.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-4-hg38 5
sbatch --output=souporcell-5.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-5-hg38 6
sbatch --output=souporcell-6.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-6-hg38 6
sbatch --output=souporcell-7.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-7-hg38 6
sbatch --output=souporcell-8.log workflow/scripts/run_souporcell.slurm Sample_5124-NM-8-hg38 6

```

### Troubleshooting specifying study design in yaml
```python

import yaml

with open("workflow/src/multiome.yaml") as file:
    try:
        config = yaml.safe_load(file)
        print(config)
    except yaml.YAMLError as exc:
        print(exc)

config = yaml.load("workflow/src/multiome.yaml")
```



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
