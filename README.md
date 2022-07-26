
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
```

### Build run design config
```bash
awk -v OFS="\t" 'NR==1 {print "Libname", $0} NR>1 {print "Sample_5124-"$1"-hg38",$0}' Sample_islet_list_for_multiomics_Batches_long_format.tsv > Sample_islet_list_for_multiomics_Batches_long_format_with_libname.tsv

Rscript workflow/scripts/build_batch_json.R --batchfile Sample_islet_list_for_multiomics_Batches_long_format_with_libname.tsv
```


x```bash
checkLibMD5sum () {
  local dir=$1
  local lib=$2
  cd ${dir}/${lib}
  md5sum -c ${lib}.md5 > ${lib}.md5check
}
checkLibMD5sum /lab/data/seqcore "5124-NM"
checkLibMD5sum /lab/data/seqcore "5125-NM"
```



### Process multiome
```bash
conda activate snakemake
bash commands.sh multiome
```


### Run nf on multiome ATAC
For each readgroup, three fastq files are required:
* the first and second insert reads ('1' and '2')
* and the read with the nuclear barcode ('index')
```bash
module load golang/1.13.5 singularity/3.5.1
cd results/multiome/nf_atac_results
nohup nextflow run -resume -params-file /lab/work/ccrober/T1D_U01/results/multiome/nf_atac_config.json --results /lab/work/ccrober/T1D_U01/results/multiome/nf_atac_results pipelines/snATACseq-NextFlow/main.nf
```


### Run nf on multiome GEX
```bash
module load golang/1.13.5 singularity/3.5.1
cd results/multiome/nf_gex_results
nohup nextflow run -resume -params-file /lab/work/ccrober/T1D_U01/results/multiome/nf_gex_config.json --chemistry multiome --results /lab/work/ccrober/T1D_U01/results/multiome/nf_gex_results pipelines/snRNAseq-NextFlow/main.nf &
```


### Process fiveprime
```bash
conda activate snakemake
bash dryrun.sh fiveprime
bash commands.sh fiveprime

## Run nf pipeline
module load golang/1.13.5 singularity/3.5.1
cd ${WORK}/results/fiveprime/5125-NM/nf_gex_results/
nohup nextflow run  -resume -params-file /lab/work/ccrober/T1D_U01/results/fiveprime/5125-NM/library_config.json -c ${WORK}/pipelines/snRNAseq-NextFlow/nextflow_5GEX.config --chemistry V2 --results /lab/work/ccrober/T1D_U01/results/fiveprime/5125-NM/nf_gex_results ${WORK}/pipelines/snRNAseq-NextFlow/main.nf
```


### Process genotype data
```bash
ln -s /lab/data/nih/2022_02_16_T1D_genotypes results/imputation/2022_02_16_T1D_genotypes
bcftools query --list-samples results/imputation/2022_02_16_T1D_genotypes/chip_genotypes/Erdos_T1DPancreas_021622_Omni2.5Exome-8_v1.3.vcf.gz > results/imputation/2022_02_16_T1D_genotypes/donors.txt

conda activate snakemake
bash commands.sh imputation
```


### Cynthia fiveprime
```bash
ln -s /lab/work/herrerab/6108-NM_resequenced /lab/work/ccrober/T1D_U01/results/fiveprime/

module load golang/1.13.5 singularity/3.5.1
singularity run workflow/envs/scafe_latest.sif pipelines/SCAFE/scripts/scafe.workflow.sc.solo --help

 scafe.download.demo.input
