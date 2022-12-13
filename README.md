
Starting a repository for storing code and documentation related to the T1D U01
collaborative project: Parker, Chen, and Collins labs

### Locations
* Multiome: /lab/data/seqcore/5124-NM
* 5' RNA-seq: /lab/data/seqcore/5125-NM
* Genotyping: /lab/data/nih/2022_02_16_T1D_genotypes

* From Ricardo's paper
* * Genotyping: /lab/data/nih/2021_07_02_T1D_genotypes/chip_genotypes
* * code for cvb4 gene quantification: Ricardo_code/notebooks/cvb4_counts.Rmd
* * Rmd notebook where he does a lot of the clustering and plotting: workflow/src/ricardo/liger_integration_all_samples_134-5.Rmd

* Code from Arushi: /lab/work/arushiv/muscle-sn/analyses_hg3/cluster_sample_count_matrices/scripts/
* * main.nf (or sample.nf) - generating cell type specific ATAC-seq peak calls
* * quantify.nf - generate barcode x feature count matrix for ATAC-seq
* *

### Load workspace
```bash
#function defined in bashrc
loadU01 () {
  export WORK=/lab/work/ccrober/T1D_U01
  module load golang/1.13.5 singularity/3.5.1
}
```

### Symlinks
``bash
ln -s /lab/data/seqcore/5124-NM
ln -s /lab/data/seqcore/5125-NM-1
ln -s /lab/data/seqcore/5125-NM_resequenced_2-8
```

### MD5 sums
```bash
cd data/5124-NM
md5sum -c 5124-NM.md5 > ${WORK}/5124-NM.md5.check

cd data/5125-NM-1
md5sum -c 5125-NM-1.md5 > ${WORK}/5125-NM-1.md5.check

cd data/5125-NM_resequenced_2-8
md5sum -c 5125-NM_resequenced_2-8.md5 > ${WORK}/5125-NM_resequenced_2-8.md5.check
```


### Master library key and run design
When new libraries come in, decide on "master sample ID" and add lines to library_key.tsv to link the seqcore sampleID to the master ID
that will be used elsewhere in the pipelines. Also update the batch_design.tsv file, which specifies which donors/conditions are
included in each batch.
```bash

#keys to update
data/batch_design.tsv
data/library_key.tsv

#convert keys to configs
workflow/scripts/make_yamls.R
```



### Make library jsons
```bash
python workflow/scripts/make_library_json.py \
  --bcl2fq_report data/5124-NM/ATAC/DemuxStats_5124-NM.csv \
  --fastq_dir data/5124-NM/ATAC \
  --library_key data/library_key.tsv \
  --modality multiome_ATAC > data/5124-NM_atac_config.json

python workflow/scripts/make_library_json.py \
  --bcl2fq_report data/5124-NM/GEX/DemuxStats_5124-NM.csv \
  --fastq_dir data/5124-NM/GEX/fastqs_5124-NM \
  --library_key data/library_key.tsv \
  --modality multiome_3GEX > data/5124-NM_gex_config.json

python workflow/scripts/make_library_json.py \
  --bcl2fq_report data/5125-NM-1/DemuxStats_5125-NM.csv \
  --fastq_dir data/5125-NM-1/fastqs_5125-NM \
  --library_key data/library_key.tsv \
  --modality 5GEX > data/5125-NM-1.json

python workflow/scripts/make_library_json.py \
  --bcl2fq_report data/5125-NM_resequenced_2-8/DemuxStats_5125-NM.csv \
  --fastq_dir data/5125-NM_resequenced_2-8/fastqs_5125-NM \
  --library_key data/library_key.tsv \
  --modality 5GEX > data/5125-NM_resequenced_2-8.json

```

### Run nf on 2020 ChenLab ATAC
#!!!!! MISSING HPAP055 !!!!!!!!!!!
```bash
cd $WORK/results/2020_Chenlab_atac/nf_atac_results
nohup nextflow -bg run -resume \
  -params-file $WORK/data/2020-Chenlab_atac_hg38.json \
  -c $WORK/results/2020_Chenlab_atac/nextflow.config \
  -work-dir $WORK/results/2020_Chenlab_atac/nf_atac_results/work \
  --results $WORK/results/2020_Chenlab_atac/nf_atac_results \
  /lab/work/ccrober/sw/snATACseq-NextFlow/main.nf

```


### Run nf on 2020 ChenLab 3GEX
# Note that the cellranger web summary report provides which chemistry was
# used for each library
```bash
cd $WORK/results/2020_Chenlab_gex_V3/nf_gex_results
nohup nextflow -bg run -resume \
  -params-file $WORK/data/2020-Chenlab_gex_V3_hg38.json \
  -c $WORK/results/2020_Chenlab_gex_V3/nextflow.config \
  -work-dir $WORK/results/2020_Chenlab_gex_V3/nf_gex_results/work \
  --results $WORK/results/2020_Chenlab_gex_V3/nf_gex_results \
  --chemistry V3 \
  /lab/work/ccrober/sw/snRNAseq-NextFlow/main.nf


cd $WORK/results/2020_Chenlab_gex_V2/nf_gex_results
nohup nextflow -bg run -resume \
  -params-file $WORK/data/2020-Chenlab_gex_V2_hg38.json \
  -c $WORK/results/2020_Chenlab_gex_V2/nextflow.config \
  -work-dir $WORK/results/2020_Chenlab_gex_V2/nf_gex_results/work \
  --results $WORK/results/2020_Chenlab_gex_V2/nf_gex_results \
  --chemistry V2 \
  /lab/work/ccrober/sw/snRNAseq-NextFlow/main.nf
```


### Run nf on multiome ATAC
```bash
nohup nextflow run -resume -params-file results/multiome/nf_atac_config.json --results results/multiome/nf_atac_results ./snATACseq-NextFlow/main.nf &
```

### Run nf on multiome GEX
```bash
nohup nextflow run -resume -params-file results/multiome/nf_gex_config.json --chemistry multiome --results results/multiome/nf_gex_results ./snRNAseq-NextFlow/main.nf &
```

### Process multiome
```bash
bash run.sh multiome
```

### Process fiveprime
```bash
bash run.sh fiveprime
```

### Demultiplex donors (barcode --> donor)
```bash
bash run.sh demultiplex.smk
```

### Liger clusters (barcode --> cell_type)
```bash
bash run.sh liger.smk
```

### Integrate to master barcode map (barcode --> donor, condition, cell_type)
```bash
bash run.sh integration.smk
```

### SCAFE
```bash
bash run.sh scafe.smk
```

### MACS2
```bash
bash run.sh macs2.smk
```

### fGWAS
```bash
#get whole islet annotations
awk 'BEGIN {OFS="\t"} {print $1, $2-150, $2+150}' results/multiome/macs2/Sample_5124-NM-1-hg38_summits.bed > results/fgwas/annotations/hg38/wholeIslet_ATAC_300b.bed

awk 'BEGIN {OFS="\t"} {print $1, $2-500, $2+500}' results/multiome/macs2/Sample_5124-NM-1-hg38_summits.bed > results/fgwas/annotations/hg38/wholeIslet_ATAC_1kb.bed


snakemake all -n \
  --snakefile workflow/src/reformat_ricardo_bed_for_fgwas.smk \
  --printshellcmds \
  --cores 1 \
  > logs/reformat_ricardo_bed_for_fgwas.log 2>&1 &

bash run.sh fgwas.smk
```

### Plots for HIRN/ASHG
```bash
sbatch workflow/scripts/run_liger_plots.slurm

Rscript workflow/scripts/liger_plots_broad_clusters.R \
  --liger_obj results/liger/multiome_5GEX/liger_obj_clusters.rds \
  --barcode_to_cell_type results/liger/multiome_5GEX/barcode_to_cluster_to_cell_type.csv \
  --outdir results/plots/multiome_5GEX
```

### Browser colors
#Endocrine == red/orange/yellow
alpha = 255,191,0
beta = 255,87,51
gamma = 242,140,40
delta = 255,172,28

#Exocrine
acinar = 128,222,234

#Other
ductal = 121,85,72
endothelial = 144,164,174
stellate = 76,175,80
immune = 171,71,188








<!-- ### Comparison of 5GEX and 3GEX for grant
```bash
#get atac features to use as starting point
cp results/multiome/nf_atac_results/work/71/e612117171963be950d9dd11336788/Sample_5124-NM-1-hg38_peaks.gappedPeak results/multiome/features

awk 'BEGIN {OFS="\t"} {print $1,$7,$8, $4}' results/multiome/features/Sample_5124-NM-1-hg38_peaks.gappedPeak > results/multiome/features/Sample_5124-NM-1-hg38_peaks.tmp.narrowPeak


### Regions for browser shots
RFX6 chr6:116,963,203-116,976,677

Intergenic ATAC peaks with and without 5GEX signal:
chr10:57Mb
SLC30A8 - multi:11,677,623-11,768,071

``` -->
