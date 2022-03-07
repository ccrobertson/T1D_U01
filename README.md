
Starting a repository for storing code and documentation related to the T1D U01
collaborative project: Parker, Chen, and Collins labs

### Data location
* Multiome: /lab/data/seqcore/5124-NM
* 5' RNA-seq: /lab/data/seqcore/5125-NM

* From Ricardo's paper
* * Genotyping: /lab/data/nih/2021_07_02_T1D_genotypes/chip_genotypes
* * snRNA-seq libraries:
* * snATAC-seq libraries:


### Universal environmental variables
```bash
export pipelines=/lab/work/ccrober/pipelines
export WORK=/lab/work/ccrober/T1D_U01
export freeze=${WORK}/results/freeze1
export scripts=${WORK}/scripts
export data=${WORK}/data
export resources=${WORK}/resources
```

### Copy code
```bash
mkdir -p ${resources}/Ricardo_code
rsync -r -v /home/albanus/analyses/2020_nih_islets_sn_t1d/commands ${resources}/Ricardo_code
rsync -r -v /home/albanus/analyses/2020_nih_islets_sn_t1d/src ${resources}/Ricardo_code
rsync -r -v /home/albanus/analyses/2020_nih_islets_sn_t1d/scripts ${resources}/Ricardo_code
rsync -r -v /home/albanus/analyses/2020_nih_islets_sn_t1d/notebooks ${resources}/Ricardo_code

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



### Prepare nextflow input jsons
```bash
#cut -d"," -f 2 /lab/data/seqcore/5124-NM/ATAC/DemuxStats_5124-NM.csv | awk 'NR>1' | sed 's/"//g'
cut -d"," -f 2 /lab/data/seqcore/5124-NM/GEX/DemuxStats_5124-NM.csv | awk 'NR>1' | sed 's/"//g;s/-3GEX//g' > ${freeze}/multiome_sample_list.tsv

python ${scripts}/build_multiome_json.py --sample-file ${freeze}/multiome_sample_list.tsv --fastq-dir ${data}/fastq_atac --modality ATAC > ${freeze}/nf_atac_config.json
python ${scripts}/build_multiome_json.py --sample-file ${freeze}/multiome_sample_list.tsv --fastq-dir ${data}/fastq_gex --modality GEX > ${freeze}/nf_gex_config.json
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


### Process multiome ATAC data
For each readgroup, three fastq files are required:
* the first and second insert reads ('1' and '2')
* and the read with the nuclear barcode ('index')
```bash
module load golang/1.13.5 singularity/3.5.1

run=nf_atac_test2
cd ${freeze}/${run}_results
nohup nextflow run -resume -params-file ${freeze}/${run}_config.json --results ${freeze}/${run}_results ${pipelines}/snATACseq-NextFlow/main.nf

run=nf_atac
cd ${freeze}/${run}_results
nohup nextflow run -resume -params-file ${freeze}/${run}_config.json --results ${freeze}/${run}_results ${pipelines}/snATACseq-NextFlow/main.nf



## Exploring what caused job to fail
#the full run made it through alignment but then failed during pruning and/or peak calling
cd ${freeze}/nf_atac_results
grep -l "CANCELLED" work/*/*/.command.log > failed_logs
while read -r line;
do
  echo ${line};
  head ${line};
  echo " ";
  tail ${line};
  echo " ";
  echo " ";
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
done < failed_logs

#Seems like macs2 ran out of memory at 2:03 AM which triggered all other jobs to be cancelled at 2:06 AM
#cd ${freeze}/nf_atac_results
#tail work/02/f4b069320e6e8a4d528d3ce5eab84b/.command.log
# INFO  @ Sat, 12 Feb 2022 01:41:13: #3 Pre-compute pvalue-qvalue table...
# slurmstepd-wolverine: error: Job 13406085 exceeded memory limit (10860173 > 8192000), being killed
# slurmstepd-wolverine: error: Exceeded job memory limit
# slurmstepd-wolverine: error: *** JOB 13406085 ON wolverine CANCELLED AT 2022-02-12T02:03:03 ***

# #what version of macs2 does the pipeline use?
# singularity exec /lab/work/ccrober/singularity-cache/porchard-default-general-20220107.img macs2 --version
# macs2 2.2.7.1
# singularity exec /lab/work/ccrober/singularity-cache/porchard-default-general-20220107.img macs2 callpeak --help
```


### Process multiome GEX data
```bash
module load golang/1.13.5 singularity/3.5.1

run=nf_gex_test
cd ${freeze}/${run}_results
nohup nextflow run -resume -params-file ${freeze}/${run}_config.json --chemistry multiome --results ${freeze}/${run}_results ${pipelines}/snRNAseq-NextFlow/main.nf

#Can find starsolo QC info for gene-level counts here:
#starsolo/Sample_test-hg38/Sample_test-hg38.Solo.out/Gene/Summary.csv

run=nf_gex
cd ${freeze}/${run}_results
nohup nextflow run -resume -params-file ${freeze}/${run}_config.json --chemistry multiome --results ${freeze}/${run}_results ${pipelines}/snRNAseq-NextFlow/main.nf
```

Checking to see if any reads align to CVB4 viral genome
```bash
cd ${freeze}/nf_gex_results/starsolo/Sample_5124-NM-1-hg38
samtools index Sample_5124-NM-1-hg38.Aligned.sortedByCoord.out.bam
samtools view Sample_5124-NM-1-hg38.Aligned.sortedByCoord.out.bam AF311939.1
```


### Filtering nuclei
```bash
snakemake -kpr -j 500 -s src/multiome.smk --configfile src/multiome.yml --cluster-config src/cluster.config.yml --cluster "sbatch -t {cluster.time} -n {cluster.N} -c {cluster.threads} --mem={cluster.mem} -o ${logdir}/slurm-%j.out --parsable" --cluster-status slurm_status.py 2>&1




### Process 5' scRNA-seq
```bash
scafe.workflow.sc.solo --help
```
