
Starting a repository for storing code and documentation related to the T1D U01
collaborative project: Parker, Chen, and Collins lab

### Data located in
# Multiome --> /lab/data/seqcore/5124-NM
# 5' RNA-seq --> /lab/data/seqcore/5125-NM


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
seqcore=/lab/data/seqcore

createSymLinks () {
  target_dir=$1
  link_dir=$2
  while read filepath; do
     filepath_array=($(echo "${filepath}" | tr '/' '\n'))
     file=${filepath_array[-1]}
     #echo ${filepath} ${file}
     ln -s --force ${filepath} ${project_folder}/${run_id}/fastq_files/${file}
  done < <(ls -1 ${fastq_folder}/*/*fastq.gz)

}


lib_mult=5124-NM
lib_fprime=5125-NM
```



### To run Multiome nextflow pipelines

For each readgroup, three fastq files are required:
* the first and second insert reads ('1' and '2')
* and the read with the nuclear barcode ('index')


```bash
module load golang/1.13.5 singularity/3.5.1

#multiome snATAC-seq
nextflow run -params-file library-config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf

#multiome snRNA-seq
nextflow run -params-file library-config.json --results ${freeze}/nf_results_rna ${pipelines}/snRNAseq-NextFlow/main.nf

```


### To process 5' scRNA-seq
```bash
scafe.workflow.sc.solo --help
```
