
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

### Prepare nextflow input jsons
```bash
#cut -d"," -f 2 /lab/data/seqcore/5124-NM/ATAC/DemuxStats_5124-NM.csv | awk 'NR>1' | sed 's/"//g'
cut -d"," -f 2 /lab/data/seqcore/5124-NM/GEX/DemuxStats_5124-NM.csv | awk 'NR>1' | sed 's/"//g;s/-3GEX//g' > ${freeze}/multiome_sample_list.tsv

python ${scripts}/build_multiome_json.py --sample-file ${freeze}/multiome_sample_list.tsv --fastq-dir ${data}/fastq_atac --modality ATAC > ${freeze}/nf_atac_config.json
python ${scripts}/build_multiome_json.py --sample-file ${freeze}/multiome_sample_list.tsv --fastq-dir ${data}/fastq_gex --modality GEX > ${freeze}/nf_gex_config.json

```


### Create test sample
```bash
cd ${data}/fastq_atac
zcat Sample_5124-NM-1-ATAC_a_R1_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R1_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R2_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R2_001.fastq.gz
zcat Sample_5124-NM-1-ATAC_a_R3_001.fastq.gz | head -n 4000000 | gzip > Sample_test_R3_001.fastq.gz
```


### Process multiome ATAC data
For each readgroup, three fastq files are required:
* the first and second insert reads ('1' and '2')
* and the read with the nuclear barcode ('index')

```bash
module load golang/1.13.5 singularity/3.5.1

#multiome snATAC-seq
#cd ${freeze}/nf_results_atac
#sbatch --time="48:00:00" --mem=105G --ntasks=12 --output=${freeze}/nf_results_atac/run_nf_atac.log --wrap="nextflow run -params-file ${scripts}/nf_atac_config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf"

cd ${freeze}/nf_results_atac
nohup nextflow run -params-file ${scripts}/nf_atac_config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf

cd ${freeze}/nf_results_atac
nohup nextflow run -params-file ${scripts}/nf_atac_config.json --results ${freeze}/nf_results_atac ${pipelines}/snATACseq-NextFlow/main.nf -resume

```


### Process multiome GEX data
```bash

module load golang/1.13.5 singularity/3.5.1

cd ${freeze}/nf_results_gex
nextflow run -params-file library-config.json --results ${freeze}/nf_results_rna ${pipelines}/snRNAseq-NextFlow/main.nf
```


### Process 5' scRNA-seq
```bash
scafe.workflow.sc.solo --help
```
