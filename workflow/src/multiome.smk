#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import pandas as pd
import math

name = config["name"]

_data = partial(os.path.join, "data")
_results = partial(os.path.join, "results", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


### Samples to process (multiome.yaml)
samples = config["samples"]
wildcard_constraints:
   sample="|".join(samples)

IONICE = 'ionice -c2 -n7'

rule all:
    input:
        #~~~~~~~~ gzip starsolo
        #expand(_results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/matrix.mtx.gz"),sample="Sample_test-hg38", feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        #expand(_results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/matrix.mtx.gz"), sample=samples, feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        #~~~~~~~~ droplet_utils
        #expand(_results("droplet_utils/{sample}/barcodes_nuclei.txt"), sample=samples),
        #~~~~~~~~ cross modality barcode filtering
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.tsv"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.tsv"), sample=samples),
        #~~~~~~~~ amulet
        #expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample=samples),
        #~~~~~~~~ prepare_rna_data
        #expand(_results("counts_by_sample_gex/{sample}/counts_nuclei.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("counts_by_sample_gex/{sample}/counts_nuclei.rds"), sample=samples),
        #~~~~~~~~ preclustering and decontamination
        #expand(_results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"), sample=samples),
        #~~~~~~~~ post decontx clustering
        #expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample="Sample_5124-NM-2-hg38"),
        expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample=samples),


rule barcode_map:
    input:
        _resources("barcode_whitelist_multiome_GEX.txt"),
        _resources("barcode_whitelist_multiome_ATAC.txt"),
    output:
        _resources("multiome_barcode_map.rds"),
        _resources("multiome_barcode_map.tsv"),
    shell:
        """
        Rscript workflow/scripts/make_multiome_barcode_map.R
        """

#Some downstream programs require gzipped output, while others need it unzipped
#so we will keep it available in both forms
rule gzip_starsolo:
    input:
        mtx = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/matrix.mtx"),
        barcodes = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/barcodes.tsv"),
        features = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/features.tsv"),
    output:
        mtx_gz = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/matrix.mtx.gz"),
        barcodes_gz = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/barcodes.tsv.gz"),
        features_gz = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/features.tsv.gz"),
        genes = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/{feature_type}/raw/genes.tsv"),
    shell:
        """
        gzip -c {input.mtx} > {output.mtx_gz}
        gzip -c {input.barcodes} > {output.barcodes_gz}
        gzip -c {input.features} > {output.features_gz}
        ln -s --relative --force {input.features} {output.genes}
        """


rule droplet_utils:
    input:
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
    output:
        _results("droplet_utils/{sample}/barcode_rank.png"),
        _results("droplet_utils/{sample}/barcodes_empty.txt"),
        _results("droplet_utils/{sample}/barcodes_nuclei.txt"),
        _results("droplet_utils/{sample}/emptydrops.rds"),
    params:
        input_10x_dir = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw"),
        outdir = _results("droplet_utils/{sample}"),
    shell:
        """
        Rscript workflow/scripts/run_dropletutils.R \
            --input_10x_dir {params.input_10x_dir} \
            --outdir {params.outdir}
        """

rule cross_modality_qc:
    input:
        gex_qc = _results("nf_gex_results/qc/{sample}-hg38.qc.txt"),
        atac_qc = _results("nf_atac_results/ataqv/single-nucleus/{sample}-hg38.txt"),
        barcode_map = _resources("multiome_barcode_map.rds"),
        droplet_utils_nuclei = _results("droplet_utils/{sample}/barcodes_nuclei.txt"),
        droplet_utils_empty = _results("droplet_utils/{sample}/barcodes_empty.txt"),
    output:
        _results("cross_modality_qc/{sample}/barcodes_nuclei.tsv"),
        _results("cross_modality_qc/{sample}/barcodes_empty.tsv"),
        _results("cross_modality_qc/{sample}/cross_modality_qc.tsv"),
        _results("cross_modality_qc/{sample}/qc_density_droplet_utils_facet_nuclei.png"),
        _results("cross_modality_qc/{sample}/qc_density_outlier_filtering.png"),
    params:
        outdir = _results("cross_modality_qc/{sample}"),
        sample = "{sample}",
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/cross_modality_qc.R \
            --gex_qc {input.gex_qc} \
            --atac_qc {input.atac_qc} \
            --barcode_map {input.barcode_map} \
            --droplet_utils_nuclei {input.droplet_utils_nuclei} \
            --droplet_utils_empty {input.droplet_utils_empty} \
            --outdir {params.outdir} \
            --sample {params.sample}
        """


# rule amulet:
#     input:
#         bam = _results("nf_atac_results/prune/{sample}-hg38.pruned.bam"),
#         qc = _results("cross_modality_qc/{sample}/cross_modality_qc.tsv"),
#         barcode_map = _resources("multiome_barcode_map.rds"),
#         autosomes = config["autosomes"],
#         blacklist = config["blacklist"],
#     output:
#         #_results("amulet/{sample}/MultipletBarcodes_01.txt"),
#         _results("amulet/{sample}/barcodes.csv"),
#     params:
#         outdir = _results("amulet/{sample}"),
#     conda:
#         "amulet"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         Rscript workflow/scripts/prep_for_amulet.R \
#                 --qc {input.qc} \
#                 --barcode_map {input.barcode_map}
#                 --outdir {params.outdir}
#         #AMULET.sh --bambc CB --bcidx 0 --cellidx 1 --iscellidx 2 {input.bam} {params.outdir}/barcodes.csv {input.autosomes} {input.blacklist} {params.outdir} $HOME/bin
#         """


#### ADD GEX DOUBLET DETECTION
# rule doublet_finder:
#     input:
#     output:
#         _results("doublet_finder/{sample}/doublet_barcodes.txt"),
#     shell:
#         """
#         Rscript workflow/scripts/run_doublet_finder.R \
#             --input
#
#         """



rule prepare_rna_counts:
    input:
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv"),
        barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei_gex.txt"),
        barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty_gex.txt"),
    output:
        counts_nuclei = _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
    params:
        input_10x_dir = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw"),
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
        counts = _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
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

#### ADD GEX DOUBLET DETECTION
rule doublet_filtering:
    input:
        seurat_obj = _results("seurat_prelim/{sample}/seurat_obj.rds"),
        doublets = _results("amulet/{sample}/MultipletBarcodes_01.txt"),
        barcode_map = _resources("multiome_barcode_map.rds"),
        counts =  _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
    output:
        counts = _results("counts_by_sample_gex/{sample}/counts_nuclei_no_doublets.rds"),
        png = _results("seurat_prelim/{sample}/seurat_prelim_umap_doublets.png"),
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/doublet_filtering.R \
            --seurat_obj {input.seurat_obj} \
            --doublets {input.doublets} \
            --barcode_map {input.barcode_map} \
            --input_counts {input.counts} \
            --output_counts {output.counts} \
            --plotfile {output.png}
        """


rule decontx_prelim:
    input:
        counts_nuclei = _results("counts_by_sample_gex/{sample}/counts_nuclei_no_doublets.rds"),
        counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
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
        counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Bam filtering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule filter_bams_for_nuclei:
    input:
        gex_bam = _results("nf_gex_results/prune/{sample}-hg38.before-dedup.bam"),
        atac_bam = _results("nf_atac_results/prune/{sample}-hg38.pruned.bam"),
        barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
    output:
        gex_bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex.bam"),
        atac_bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_atac.bam"),
    params:
        outdir = _results("bam_pass_qc_barcodes/{sample}"),
    container:
        "workflow/envs/arushi_general.simg"
    shell:
        """
        awk 'NR>1 {{print $1}}' {input.barcodes_nuclei} > {params.outdir}/pass_qc_barcodes_gex.txt
        awk 'NR>1 {{print $2}}' {input.barcodes_nuclei} > {params.outdir}/pass_qc_barcodes_atac.txt
        samtools view -h -b -@ 10 -D CB:{params.outdir}/pass_qc_barcodes_gex.txt {input.gex_bam} > {output.gex_bam}
        samtools view -h -b -@ 10 -D CB:{params.outdir}/pass_qc_barcodes_atac.txt {input.atac_bam} > {output.atac_bam}
        samtools index {output.gex_bam}
        samtools index {output.atac_bam}
        """


rule define_ambient_gene_regions:
    input:
        counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
        counts_nuclei = _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
        gtf = _resources("hg38/gencode.v39.annotation.gtf"),
    output:
        regions = _results("bam_pass_qc_barcodes/{sample}/ambient_gene_regions_{p}.bed"),
    params:
        outdir = _results("bam_pass_qc_barcodes/{sample}"),
        p = "{p}",
    conda:
        "Renv"
    wildcard_constraints:
        p = "auto|0.01|0.3"
    shell:
        """
        Rscript workflow/scripts/define_ambient_gene_regions.R \
            --counts_empty {input.counts_empty} \
            --counts_nuclei {input.counts_nuclei} \
            --gtf {input.gtf} \
            --p {params.p} \
            --outdir {params.outdir}
        """

rule filter_bams_for_ambient_genes:
    input:
        #gex_bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex.bam"),
        gex_bam = _results("nf_gex_results/prune/{sample}-hg38.before-dedup.bam"),
        regions = _results("bam_pass_qc_barcodes/{sample}/ambient_gene_regions_{p}.bed"),
    output:
        gex_bam = _results("bam_pass_qc_barcodes/{sample}/prune_barcodes_gex_exclude_ambient_genes_{p}.bam"),
        gex_bam_discarded = _results("bam_pass_qc_barcodes/{sample}/prune_barcodes_gex_discarded_{p}.bam"),
    container:
        "workflow/envs/arushi_general.simg"
    wildcard_constraints:
        p = "auto|0.01|0.3"
    shell:
        """
        samtools view -h -b -@ 10 -L {input.regions} -U {output.gex_bam} {input.gex_bam} > {output.gex_bam_discarded}
        samtools index {output.gex_bam}
        """
