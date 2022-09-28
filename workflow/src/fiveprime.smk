#!/usr/bin/env python3

from os.path import join
import os
import pandas as pd
from functools import partial
import math
import glob


name = config["name"]
_data = partial(os.path.join, "data")
_results = partial(os.path.join, "results", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


configfile: _data("nandini_run_design.json")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def iterate_batches_in_run(run):
    batches = list(config["runs"][run]["batches_in_run"].keys())
    return batches

def get_uniq_sample_ids():
    s = []
    runs = config["runs"].keys()
    for r in runs:
        batches = iterate_batches_in_run(r)
        samples = ["Sample_"+r+"_"+b for b in batches]
        s = s + samples
    return s

def fastq_dir_by_run(run):
    dir = config['runs'][run]['fastq_dir']
    return dir

def fastqs_by_sample(run, batch):
    f = open(config['runs'][run]['bcl2fq_report'])
    for line in f:
        line = line.replace('"','')
        vals = line.split(',')
        Sample_ID = config["runs"][run]["batches_in_run"][batch]
        if vals[1] == Sample_ID:
            barcode = vals[3]
            dir = config['runs'][run]['fastq_dir']
            R1 = glob.glob(dir + "/*" + barcode + "*" + "R1_001.fastq.gz")[0]
            R2 = glob.glob(dir + "/*" + barcode + "*" + "R2_001.fastq.gz")[0]
    return [R1, R2]

def fastq_prefix_by_sample(run, batch):
    fullpath = fastqs_by_sample(run, batch)[0]
    file = fullpath.split('/')[len(fullpath.split('/'))-1]
    prefix = re.sub("_S[0-9]+_R1_001.fastq.gz","", file)
    return prefix


#print(fastqs_by_sample("5125-NM-1","NM-1"))
#print(" ".join(fastqs_by_sample("5125-NM-1","NM-1")))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 5GEX runs to process from fiveprime.yaml
runs = config["runs"].keys()
samples = get_uniq_sample_ids()
#print(samples)

### batch information from nandini_run_design.json
batches = config["batches"].keys()

### restrict wildcards to existing runs, batches, and samples
wildcard_constraints:
   run="|".join(runs),
   sample="|".join(samples),
   batch="|".join(batches),


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Run
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule all:
    input:
        #~~~~~~~~ cellranger
        #expand(_results("cellranger/{sample}.done.flag"), sample="Sample_5125-NM-2-8_NM-2"),
        #expand(_results("cellranger/{sample}.done.flag"), sample=samples),
        #~~~~~~~~ sort cellranger bam
        #expand(_results("cellranger_bam_sorted/{sample}/possorted_genome_bam_autosomes.bam"), sample=samples),
        #expand(_results("cellranger_bam_sorted/{sample}/possorted_genome_bam_autosomes.bam"), sample="Sample_5125-NM-2-8_NM-2"),
        #~~~~~~~~ starsolo
        #expand(_results("starsolo/{sample}/starsolo.Aligned.sortedByCoord.out.bam"), sample=samples),
        #~~~~~~~~ gzip cellranger
        #expand(_results("cellranger/{sample}/outs/raw_feature_bc_matrix/genes.tsv"), sample="Sample_5125-NM-1_NM-1"),
        #~~~~~~~~ gzip starsolo
        #expand(_results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/genes.tsv"),sample="Sample_5125-NM-1_NM-1", feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        #~~~~~~~~ qc
        #expand(_results("qc/{sample}.qc.txt"), sample="Sample_5125-NM-1_NM-1"),
        expand(_results("qc/{sample}.qc.metrics.png"), sample=samples),
        #~~~~~~~~ seurat_round3
        expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample=samples),


rule cellranger:
    input:
        fastqs = lambda wildcards: fastqs_by_sample(wildcards.run, wildcards.batch),
    output:
        flag = _results("cellranger/Sample_{run}_{batch}.done.flag"),
    params:
        id = "Sample_{run}_{batch}",
        fastq_dir = lambda wildcards: fastq_dir_by_run(wildcards.run),
        fastq_prefix = lambda wildcards: fastq_prefix_by_sample(wildcards.run, wildcards.batch),
        outdir = _results("cellranger"),
    shell:
        """
        fastq_dir=$PWD/{params.fastq_dir}
        flag=$PWD/{output.flag}
        cd {params.outdir}
        cellranger count --id={params.id} \
           --transcriptome=/lab/work/ccrober/sw/cellranger/refdata-gex-GRCh38-2020-A \
           --fastqs=$fastq_dir \
           --sample={params.fastq_prefix} \
           --chemistry=SC5P-PE \
           --localcores=10 \
           --localmem=50
        touch $flag
        """

rule gzip_cellranger:
    input:
        features_gz = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/features.tsv.gz"),
        barcodes_gz = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz"),
        matrix_gz = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz"),
    output:
        features = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/features.tsv"),
        barcodes = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv"),
        matrix = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/matrix.mtx"),
        genes = _results("cellranger/{sample}/outs/raw_feature_bc_matrix/genes.tsv"),
    shell:
        """
        gzip -d -c {input.matrix_gz} > {output.matrix}
        gzip -d -c {input.barcodes_gz} > {output.barcodes}
        gzip -d -c {input.features_gz} > {output.features}
        ln -s --relative --force {output.features} {output.genes}
        """

rule cellranger_bam_sorted:
    input:
        bam = _results("cellranger/{sample}/outs/possorted_genome_bam.bam"),
    output:
        sam =  _results("cellranger_bam_sorted/{sample}/possorted_genome_bam_autosomes.sam"),
        bam =  _results("cellranger_bam_sorted/{sample}/possorted_genome_bam_autosomes.bam"),
        bai = _results("cellranger_bam_sorted/{sample}/possorted_genome_bam_autosomes.bam.bai"),
    params:
        autosomes = expand("chr{chr}", chr=range(1,23)),
    container:
        "workflow/envs/arushi_general.simg",
    shell:
        """
        samtools view -H {input.bam} | grep -v "@SQ" > {output.sam}
        samtools view -H {input.bam} | grep "@SQ" | workflow/scripts/fix_cellranger_bam_header.sh >> {output.sam}
        samtools view -@ 10 {input.bam} {params.autosomes} >>  {output.sam}
        samtools view -b -@ 10 {output.sam} > {output.bam}
        samtools index -b {output.bam}
        """


rule starsolo:
    input:
        fastqs = lambda wildcards: fastqs_by_sample(wildcards.run, wildcards.batch),
    output:
        bam = _results("starsolo/Sample_{run}_{batch}/starsolo.Aligned.sortedByCoord.out.bam"),
        mtx = expand(_results("starsolo/Sample_{{run}}_{{batch}}/starsolo.Solo.out/{feature_type}/raw/matrix.mtx"),feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        barcodes = expand(_results("starsolo/Sample_{{run}}_{{batch}}/starsolo.Solo.out/{feature_type}/raw/barcodes.tsv"),feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        features = expand(_results("starsolo/Sample_{{run}}_{{batch}}/starsolo.Solo.out/{feature_type}/raw/features.tsv"),feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
    params:
        prefix = _results("starsolo/Sample_{run}_{batch}/starsolo."),
        genomeDir = _resources("hg38/hg38_cvb4"),
        sjdbGTFfile = _resources("hg38/gencode.v39.annotation.CVB4.gtf"),
        soloCBwhitelist = _resources("barcode_whitelist_5GEX.txt"),
        outSAMattributes = "NH HI nM AS CR CY CB UR UY UB sM GX GN",
    shell:
        """
        singularity exec --bind /lab workflow/envs/star_2.7.10a.sif STAR --soloBarcodeReadLength 0 \
              --runThreadN 10 \
              --outFileNamePrefix {params.prefix} \
              --genomeLoad NoSharedMemory \
              --limitBAMsortRAM 60000000000 \
              --runRNGseed 8675309 \
              --readFilesCommand gunzip -c \
              --outSAMattributes {params.outSAMattributes} \
              --outSAMtype BAM SortedByCoordinate \
              --genomeDir {params.genomeDir} \
              --outSAMunmapped Within KeepPairs \
              --sjdbGTFfile {params.sjdbGTFfile} \
              --soloFeatures Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS \
              --soloMultiMappers Uniform PropUnique EM Rescue \
              --soloUMIfiltering MultiGeneUMI \
              --soloCellFilter None \
              --soloCBmatchWLtype 1MM_multi_pseudocounts \
              --soloCBwhitelist {params.soloCBwhitelist} \
              --soloBarcodeMate 1 \
              --clip5pNbases 39 0 \
              --soloType CB_UMI_Simple \
              --soloCBstart 1 \
              --soloCBlen 16 \
              --soloUMIstart 17 \
              --soloUMIlen 10 \
              --readFilesIn {input.fastqs}
        """


rule gzip_starsolo:
    input:
        mtx = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/matrix.mtx"),
        barcodes = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/barcodes.tsv"),
        features = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/features.tsv"),
    output:
        mtx_gz = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/matrix.mtx.gz"),
        barcodes_gz = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/barcodes.tsv.gz"),
        features_gz = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/features.tsv.gz"),
        genes = _results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/genes.tsv"),
    shell:
        """
        gzip -c {input.mtx} > {output.mtx_gz}
        gzip -c {input.barcodes} > {output.barcodes_gz}
        gzip -c {input.features} > {output.features_gz}
        ln -s --relative --force {input.features} {output.genes}
        """


rule qc:
    input:
        bam = _results("starsolo/{sample}/starsolo.Aligned.sortedByCoord.out.bam"),
        mtx = _results("starsolo/{sample}/starsolo.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
        barcodes = _results("starsolo/{sample}/starsolo.Solo.out/GeneFull_ExonOverIntron/raw/barcodes.tsv"),
    output:
        qc = _results("qc/{sample}.qc.txt"),
    shell:
        """
        singularity exec --bind /lab workflow/envs/porchard_snatac_general_20220107.sif python workflow/scripts/qc-from-starsolo.py {input.bam} {input.mtx} {input.barcodes} > {output.qc}
        """


rule qc_plots:
    input:
        metrics = _results("qc/{sample}.qc.txt"),
    output:
        _results("qc/{sample}.qc.metrics.png"),
    params:
        prefix = _results("qc/{sample}.qc."),
    shell:
        """
        singularity exec --bind /lab workflow/envs/porchard_snatac_general_20220107.sif python workflow/scripts/plot-qc-metrics.py --prefix {params.prefix} {input.metrics}
        """


rule droplet_utils:
    input:
        _results("starsolo/{sample}/starsolo.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
    output:
        _results("droplet_utils/{sample}/barcode_rank.png"),
        _results("droplet_utils/{sample}/barcodes_empty.txt"),
        _results("droplet_utils/{sample}/barcodes_nuclei.txt"),
        _results("droplet_utils/{sample}/emptydrops.rds"),
    params:
        input_10x_dir = _results("starsolo/{sample}/starsolo.Solo.out/GeneFull_ExonOverIntron/raw"),
        outdir = _results("droplet_utils/{sample}"),
    shell:
        """
        Rscript workflow/scripts/run_dropletutils.R \
            --input_10x_dir {params.input_10x_dir} \
            --outdir {params.outdir}
        """

rule prepare_rna_counts:
    input:
        _results("starsolo/{sample}/starsolo.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
        barcodes_nuclei = _results("droplet_utils/{sample}/barcodes_nuclei.txt"),
        barcodes_empty = _results("droplet_utils/{sample}/barcodes_empty.txt"),
    output:
        counts_nuclei = _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_protein_coding/{sample}/counts_empty.rds"),
    params:
        input_10x_dir = _results("starsolo/{sample}/starsolo.Solo.out/GeneFull_ExonOverIntron/raw"),
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/prepare_rna_counts.R \
            --input_10x_dir {params.input_10x_dir} \
            --barcodes_nuclei {input.barcodes_nuclei} \
            --barcodes_empty {input.barcodes_empty} \
            --counts_nuclei {output.counts_nuclei} \
            --counts_empty {output.counts_empty}
        """


rule seurat_prelim:
    input:
        counts = _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
    output:
        _results("seurat_prelim/{sample}/seurat_obj.rds"),
        _results("seurat_prelim/{sample}/seurat_clusters.csv"),
    params:
        outdir = _results("seurat_prelim/{sample}"),
        resolution = 0.1,
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/run_seurat.R \
            --counts {input.counts} \
            --resolution {params.resolution} \
            --outdir {params.outdir}
        """

# rule doublet_filtering:
#     input:
#         seurat_obj = _results("seurat_prelim/{sample}/seurat_obj.rds"),
#         doublets = "results/demultiplex/demuxlet-unfiltered/{sample}/doublets.txt",
#         map = _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
#         counts =  _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
#     output:
#         counts = _results("counts_protein_coding/{sample}/counts_nuclei_no_doublets.rds"),
#         png = _results("seurat_prelim/{sample}/seurat_prelim_umap_doublets.png"),
#     conda:
#         "Renv"
#     shell:
#         """
#         Rscript workflow/scripts/doublet_filtering.R \
#             --seurat_obj {input.seurat_obj} \
#             --doublets {input.doublets} \
#             --map {input.map} \
#             --input_counts {input.counts} \
#             --output_counts {output.counts} \
#             --plotfile {output.png}
#         """


rule decontx_prelim:
    input:
        #counts_nuclei = _results("counts_protein_coding/{sample}/counts_nuclei_no_doublets.rds"),
        counts_nuclei = _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_protein_coding/{sample}/counts_empty.rds"),
        clusters = _results("seurat_prelim/{sample}/seurat_clusters.csv"),
    output:
        _results("decontx_prelim/{sample}/counts_low_contamination_raw.rds"),
    params:
        outdir = _results("decontx_prelim/{sample}"),
        max_contamination = 0.3,
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/run_decontx.R \
            --counts_nuclei {input.counts_nuclei} \
            --counts_empty {input.counts_empty} \
            --clusters {input.clusters} \
            --max_contamination {params.max_contamination} \
            --outdir {params.outdir}
        """

rule seurat_round2:
    input:
        counts = _results("decontx_prelim/{sample}/counts_low_contamination_raw.rds"),
    output:
        _results("seurat_round2/{sample}/seurat_obj.rds"),
        _results("seurat_round2/{sample}/seurat_clusters.csv"),
    params:
        outdir = _results("seurat_round2/{sample}"),
        resolution = 0.1,
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/run_seurat.R \
            --counts {input.counts} \
            --resolution {params.resolution} \
            --outdir {params.outdir}
        """


rule decontx_round2:
    input:
        counts_nuclei = _results("decontx_prelim/{sample}/counts_low_contamination_raw.rds"),
        counts_empty = _results("counts_protein_coding/{sample}/counts_empty.rds"),
        clusters = _results("seurat_round2/{sample}/seurat_clusters.csv"),
    output:
        results = _results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"),
    conda:
        "Renv"
    params:
        outdir = _results("decontx_round2/{sample}"),
        max_contamination = 0.3,
    shell:
        """
        Rscript workflow/scripts/run_decontx.R \
            --counts_nuclei {input.counts_nuclei} \
            --counts_empty {input.counts_empty} \
            --clusters {input.clusters} \
            --max_contamination {params.max_contamination} \
            --outdir {params.outdir}
        """


rule seurat_round3:
    input:
        counts = _results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"),
    output:
        _results("seurat_round3/{sample}/seurat_obj.rds"),
        _results("seurat_round3/{sample}/seurat_clusters.csv"),
    params:
        outdir = _results("seurat_round3/{sample}"),
        resolution = 0.1,
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/run_seurat.R \
            --counts {input.counts} \
            --resolution {params.resolution} \
            --outdir {params.outdir}
        """
