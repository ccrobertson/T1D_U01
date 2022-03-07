#!/bin/bash
#SBATCH --time="48:00:00"
#SBATCH --mem=40000M
#SBATCH --output=hg38_cvb4_star_index.log
#SBATCH --ntasks=10

## Note on memory usage
#I was able to generate the hg38+cvb4 hybrid STAR index
#and the max memory use was reported at 37.7GB, so 40GB allocation is probably just the right amount

##Note on calling singularity
#Since the environmental variable SINGULARITY_CACHEDIR is defined in my bashrc (set to /lab/work/ccrober/singularity-cache),
#using the "library://porchard/..." syntax below works.
#However, it causes the image to be stored in /lab/work/ccrober/singularity-cache/cache/library
#instead of using the already cached image  /lab/work/ccrober/singularity-cache/porchard-default-star-2.7.9a.img
#I haven't figured out how to use the already cached images generated when running the NF pipelines

module load golang/1.13.5 singularity
singularity exec library://porchard/default/star:2.7.9a \
 STAR --runMode genomeGenerate  \
 --genomeDir ${resources}/hg38/hg38_cvb4 \
 --genomeFastaFiles ${resources}/hg38/hg38_cvb4.fa \
 --sjdbGTFfile ${resources}/hg38/gencode.v39.annotation.CVB4.gtf \
--sjdbOverhang 150 --runThreadN 10


#--sjdbOverhang specifies the length of the genomic sequence around the annotated junction
#to be used in constructing the splice junctions database. Ideally, this length should be equal
#to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
#2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
#ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as
#well as the ideal value.
