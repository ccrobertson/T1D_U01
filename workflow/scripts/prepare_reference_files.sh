
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TSS file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prepare TSS file (formatted same as Peter's hg19 file: /home/porchard/github/ataqv/data/tss/hg19.tss.refseq.bed.gz)
## seems like there is more carefully generated TSS reference here: http://reftss.clst.riken.jp (DOI: 10.1016/j.jmb.2019.04.045) but their server seems to be down
cd $resources/hg38
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.sql #for column headers
mv refGene.txt.gz hg38.refGene.txt.gz
Rscript ${scripts}/make_tss.R
gzip hg38.refGene.tss.bed


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Hybrid reference genome (hg38 + CVB4)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### obtain CVB4 genome (GenBank AF311939.1)
#install NCBI E-utilities following these [instructions](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${PATH}:${HOME}/edirect
esearch -db nucleotide -query "AF311939.1" | efetch -format fasta > ${resources}/CVB4/AF311939.1.fasta


# ### interpreting FASTA headers
# # hg39 FASTA header
# >chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38
# AC: Locus identifier / Accession number at Genbank & NCBI
# gi: GI number (sometimes written in lower case, "gi") is simply a series of digits that are assigned consecutively to each sequence record processed by NCBI. The GI number bears no resemblance to the Accession number of the sequence record.
# LN: Length of the contig
# rl: role in the assembly, (chromosome, unplaced, EBV decoy)
# M5: MD5 checksum (file integrity)
# AS: Assemby version
#
# #cvb4 FASTA header
# >AF311939.1 Human coxsackievirus B4 strain E2 variant, complete genome


### append to hg38.fa as a separate chromosome
awk '$1 ~ /^>/ {print $1" rl:Chromosome"} $1 !~ /^>/' ${resources}/CVB4/AF311939.1.fasta > ${resources}/CVB4/AF311939.1_stripped_header.fasta
cat ${resources}/hg38/hg38.fa ${resources}/CVB4/AF311939.1_stripped_header.fasta > ${resources}/hg38/hg38_cvb4.fa


# ### index cvb4 fasta (just for testing)
# cd ${resources}/CVB4
# samtools faidx AF311939.1_stripped_header.fasta
# java -jar /lab/work/ccrober/sw/picard/picard.jar CreateSequenceDictionary R=AF311939.1_stripped_header.fa O=AF311939.1_stripped_header.dict


### index hybrid fasta
cd ${resources}/hg38
sbatch --time="48:00:00" --mem=20000M --output=hg38_cvb4_fai_index.log --wrap="samtools faidx hg38_cvb4.fa"
java -jar /lab/work/ccrober/sw/picard/picard.jar \
    CreateSequenceDictionary \
    R=hg38_cvb4.fa \
    O=hg38_cvb4.dict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## BWA index
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arushi set up hg38 index files for RNA seq with ERCC: /lab/work/arushiv/hg38_references/README.md
# Description of the hg38+ERCC reference created for TOPMed is here: https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
# This may be a useful roadmap for creating the hg38+coxsackievirus reference

## Prepare BWA index
cd ${resources}/hg38
sbatch --time="48:00:00" --mem=20000M --output=hg38_bwa_index.log --wrap="bwa index hg38.fa"

## Prepare BWA index with hg38+CVB4 reference
cd ${resources}/hg38
sbatch --time="48:00:00" --mem=20000M --output=hg38_cvb4_bwa_index.log --wrap="bwa index hg38_cvb4.fa"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prepare STAR index with hg38+CVB4 reference and GTF files
# Note make sure index is built using same version of STAR as is used in snRNA-seq pipeline singularity container
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### obtain gencode annotation
cd ${resources}/hg38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
bgzip -d gencode.v39.annotation.gtf.gz


### get CVB4 annotation file via GenBank browser.
#Search for AF311939.1 (https://www.ncbi.nlm.nih.gov/nuccore/AF311939.1#feature_AF311939.1)
#Get the GFF3 file using: Sent to -> Complete record -> File -> Format GFF3
#Edit third column from "region" to "gene"
cd ${resources}/CVB4
awk '$1 !~ /^#/' sequence.gff3 | awk -v OFS='\t' '$3=="region" {$3="gene"; print $0} $3=="CDS" {print}' > AF311939.1.gtf

### append to CVB4 GTF to gencode V38 GTF file
cd ${resources}/hg38
cat gencode.v39.annotation.gtf ${resources}/CVB4/AF311939.1.gtf > gencode.v39.annotation.CVB4.gtf

### create updated chrom sizes file
cd ${resources}/hg38
cat hg38.chrom.sizes <(awk -v OFS='\t' '$3=="gene" {print $1, $5}' ${resources}/CVB4/AF311939.1.gtf) > hg38_cvb4.chrom.sizes


### generate index
#cd /lab/work/ccrober/singularity-images
#singularity pull library://porchard/default/star:2.7.9a
cd ${resources}/hg38
sbatch ${scripts}/generate_star_index.slurm
