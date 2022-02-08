
Starting a repository for storing code and documentation related to the T1D U01
collaborative project: Parker, Chen, and Collins lab

### Data location
* Multiome --> /lab/data/seqcore/5124-NM
* 5' RNA-seq --> /lab/data/seqcore/5125-NM


### Universal environmental variables
```bash
export pipelines=/lab/work/ccrober/pipelines
export WORK=/lab/work/ccrober/T1D_U01
export freeze=${WORK}/analysis/freeze1
export scripts=${WORK}/scripts
export data=${WORK}/data
export resources=${WORK}/resources
```


### Check md5
```bash
checkLibMD5sum () {
  local dir=$1
  local lib=$2
  cd ${dir}/${lib}
  md5sum -c ${lib}.md5 > ${WORK}/${lib}.md5check
}

seqcore=/lab/data/seqcore
checkLibMD5sum ${seqcore} "5124-NM"
checkLibMD5sum ${seqcore} "5125-NM"
```


### Create symlinks to fastq files
```bash
createSymLinks () {
  local target_dir=$1
  local link_dir=$2
  local sample=$3
  ln -s --force ${target_dir}/${sample}/*R1_001.fastq.gz ${link_dir}/${sample}_R1_001.fastq.gz
  ln -s --force ${target_dir}/${sample}/*R2_001.fastq.gz ${link_dir}/${sample}_R2_001.fastq.gz
  ln -s --force ${target_dir}/${sample}/*R3_001.fastq.gz ${link_dir}/${sample}_R3_001.fastq.gz
}

target_dir=/lab/data/seqcore/5124-NM/ATAC
link_dir=${data}/fastq_atac
createSymLinks ${target_dir} ${link_dir} Sample_5124-NM-1-ATAC_a
createSymLinks ${target_dir} ${link_dir} Sample_5124-NM-1-ATAC_b
```

### Create test sample
```bash
cd ${data}/fastq_atac
zcat Sample_5124-NM-1-ATAC_a_R1_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R1_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R2_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R2_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R3_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R3_001.fastq.gz
```


### Process multiome data
For each readgroup, three fastq files are required:
* the first and second insert reads ('1' and '2')
* and the read with the nuclear barcode ('index')

```bash
module load golang/1.13.5 singularity/3.5.1
export NXF_SINGULARITY_CACHEDIR=${WORK}/singularity

#multiome snATAC-seq
#cd ${freeze}/nf_results_atac
#sbatch --time="48:00:00" --mem=105G --ntasks=12 --output=${freeze}/nf_results_atac/run_nf_atac.log --wrap="nextflow run -params-file ${scripts}/nf_atac_config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf"

cd ${freeze}/nf_results_atac
nohup nextflow run -params-file ${scripts}/nf_atac_config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf

cd ${freeze}/nf_results_atac
nohup nextflow run -params-file ${scripts}/nf_atac_config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf -resume


#troubleshooting test run
samtools view -c /lab/work/ccrober/T1D_U01/analysis/freeze1/nf_results_atac/work/20/6e24a3880ba63d21b85f5c7c346ce8/Sample_test-L001-hg38.bam
2000024
samtools view -c /lab/work/ccrober/T1D_U01/analysis/freeze1/nf_results_atac/work/06/4fcf1add3565c0d2b7a3be2442b17a/Sample_test-hg38.pruned.bam
0



#multiome snRNA-seq
nextflow run -params-file library-config.json --results ${freeze}/nf_results_rna ${pipelines}/snRNAseq-NextFlow/main.nf
```


### Process 5' scRNA-seq
```bash
scafe.workflow.sc.solo --help
```
