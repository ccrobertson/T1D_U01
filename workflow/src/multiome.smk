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
        #~~~~~~~~ run dropkick
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"), sample="Sample_5124-NM-2-hg38", method="multiotsu", min_genes=[100]),
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"), sample=samples, method="multiotsu", min_genes=[100]),
        #~~~~~~~~ dropkick_plots
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick_score.png"), sample="Sample_5124-NM-2-hg38", method="multiotsu", min_genes=[100]),
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick_score.png"), sample=samples, method="multiotsu", min_genes=[100]),
        #~~~~~~~~ droplet_utils
        expand(_results("droplet_utils/{sample}/barcodes_nuclei.txt"), sample=samples),
        #~~~~~~~~ cross modality prefiltering
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.txt"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.txt"), sample=samples),
        #~~~~~~~~ amulet
        #expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample=samples),
        #~~~~~~~~ cross modality barcode filtering
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered_flagoutliers.h5ad"), sample=samples),
        #expand(_results("cross_modality_qc/{sample}/qc_grid_prefiltered_density_by_QuantileFilter_label.png"), sample=samples),
        #expand(_results("cross_modality_qc/{sample}/qc_grid_prefiltered_density_QuantileFilter_keep.png"), sample=samples),
        #expand(_results("cross_modality_qc/{sample}/barcodes_nuclei.txt"), sample=samples),
        #~~~~~~~~ prepare_rna_data
        #expand(_results("counts_by_sample_gex/{sample}/counts_nuclei.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("counts_by_sample_gex/{sample}/counts_nuclei.rds"), sample=samples),
        #~~~~~~~~ preclustering and decontamination
        #expand(_results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"), sample=samples),
        #~~~~~~~~ post decontx clustering
        #expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample=samples),


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
    shell:
        """
        gzip -c {input.mtx} > {output.mtx_gz}
        gzip -c {input.barcodes} > {output.barcodes_gz}
        gzip -c {input.features} > {output.features_gz}
        """


rule dropkick:
    input:
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx.gz"),
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/barcodes.tsv.gz"),
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv.gz"),
    output:
        _results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"),
        _results("dropkick_{method}_{min_genes}/{sample}/dropkick.csv"),
    conda:
        "dropkick"
        # NOTE: Snakemake cannot build dropkick environment for envs/dropkick.yaml
        # because dropkick is not available via pip
        #library://porchard/default/dropkick:20220225
    params:
        indir = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw"),
        outdir = _results("dropkick_{method}_{min_genes}/{sample}"),
        cpus = 2,
        min_genes = "{min_genes}",
        method = "{method}",
    shell:
        """
        python -u workflow/scripts/run_dropkick.py --indir {params.indir} --outdir {params.outdir} --cpus {params.cpus} --min_genes {params.min_genes} --method {params.method}
        """


rule dropkick_plots:
    input:
        ad = _results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"),
    output:
        _results("dropkick_{method}_{min_genes}/{sample}/dropkick_score.png"),
        _results("dropkick_{method}_{min_genes}/{sample}/dropkick_coef.png"),
    params:
        outdir = _results("dropkick_{method}_{min_genes}/{sample}"),
    conda:
        "dropkick"
    shell:
        """
        python -u workflow/scripts/run_dropkick_plots.py --ad {input.ad} --outdir {params.outdir}
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

rule cross_modality_prefiltering:
    input:
        gex = _results("nf_gex_results/qc/{sample}-hg38.qc.txt"),
        atac = _results("nf_atac_results/ataqv/single-nucleus/{sample}-hg38.txt"),
        dropkick = _results("dropkick_multiotsu_100/{sample}/dropkick.csv"),
        barcode_whitelist_gex = _resources("barcode_whitelist_multiome_GEX.txt"),
        barcode_whitelist_atac = _resources("barcode_whitelist_multiome_ATAC.txt"),
    output:
        _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
        _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.txt"),
        _results("cross_modality_qc/{sample}/barcodes_prefiltered.txt"),
        _results("cross_modality_qc/{sample}/barcodes_prefiltered_for_amulet.csv"),
        _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
        _results("cross_modality_qc/{sample}/qc_reads_unfiltered.png"),
        _results("cross_modality_qc/{sample}/qc_reads_prefiltered.png"),
        _results("cross_modality_qc/{sample}/qc_grid_unfiltered.png"),
        _results("cross_modality_qc/{sample}/qc_grid_prefiltered.png"),
    params:
        outdir = _results("cross_modality_qc/{sample}"),
        min_reads = 100,
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/cross_modality_prefiltering.R \
            --input_gex {input.gex} \
            --input_atac {input.atac} \
            --dropkick {input.dropkick} \
            --outdir {params.outdir} \
            --sample {wildcards.sample} \
            --min_reads {params.min_reads} \
            --barcode_whitelist_gex {input.barcode_whitelist_gex} \
            --barcode_whitelist_atac {input.barcode_whitelist_atac}
        """


rule doublet_detection:
    input:
        bam = _results("nf_atac_results/prune/{sample}-hg38.pruned.bam"),
        csv = _results("cross_modality_qc/{sample}/barcodes_prefiltered_for_amulet.csv"),
        autosomes = config["autosomes"],
        blacklist = config["blacklist"],
    output:
        #_results("amulet/{sample}/Overlaps.txt"),
        _results("amulet/{sample}/MultipletBarcodes_01.txt"),
    params:
        outdir = _results("amulet/{sample}"),
    conda:
        "amulet"
    shell:
        """
        mkdir -p {params.outdir}
        AMULET.sh --bambc CB --bcidx 0 --cellidx 1 --iscellidx 2 {input.bam} {input.csv} {input.autosomes} {input.blacklist} {params.outdir} $HOME/bin
        """

#NOTE: must generate H5AD file using this singularity container so that versions of scanpy are compatible
rule anndata_for_qc:
    input:
        mtx = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx.gz"),
        barcodes = _results("cross_modality_qc/{sample}/barcodes_prefiltered.txt"),
        qcstats = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.txt"),
    output:
        h5ad_file = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.h5ad"),
    params:
        starsolodir = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw"),
    conda:
        "dropkick"
    shell:
        """
        python workflow/scripts/prepare_anndata_qc.py --starsolodir {params.starsolodir} --barcodes {input.barcodes} --qcstats {input.qcstats} --h5ad_file {output.h5ad_file}
        """

rule outlier_detection:
    input:
        h5ad = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.h5ad"),
    output:
        h5ad = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered_flagoutliers.h5ad"),
        csv = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered_flagoutliers.csv"),
    params:
        metadata_columns = 'umis,hqaa,frac_mt_gex,frac_mt_atac,tss_enrichment,max_fraction_reads_from_single_autosome,dropkick_score',
        contamination = 0.1,
    conda:
         "dropkick"
    shell:
        """
        python workflow/scripts/outlier_detection.py \
                --input_h5ad {input.h5ad} \
                --metadata_columns {params.metadata_columns} \
                --contamination {params.contamination} \
                --output_h5ad {output.h5ad} \
                --output_csv {output.csv}
        """


rule cross_modality_visualization:
    input:
        csv = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered_flagoutliers.csv"),
    output:
        _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
        _results("cross_modality_qc/{sample}/qc_grid_prefiltered_density_by_dropkick_label.png"),
        _results("cross_modality_qc/{sample}/qc_grid_prefiltered_density_by_IsolationForest_label.png"),
        _results("cross_modality_qc/{sample}/qc_grid_prefiltered_density_by_QuantileFilter_label.png"),
        _results("cross_modality_qc/{sample}/qc_grid_prefiltered_density_QuantileFilter_keep.png"),
    params:
        outdir = _results("cross_modality_qc/{sample}"),
        sample = "{sample}",
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/cross_modality_visualization.R \
            --input_csv {input.csv} \
            --outdir {params.outdir} \
            --sample {params.sample}
        """

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

#### UPDATE BARCODE INPUT FILES HERE TO BE ONLY GEX BARCODES
rule prepare_rna_counts:
    input:
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv"),
        barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei_gex.txt"),
        barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty_gex.txt"),
    output:
        _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw/genes.tsv"),
        counts_nuclei = _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
    params:
        input_10x_dir = _results("nf_gex_results/starsolo/{sample}-hg38/{sample}-hg38.Solo.out/GeneFull_ExonOverIntron/raw"),
    conda:
        "Renv"
    shell:
        """
        ln -s --relative --force {params.input_10x_dir}/features.tsv {params.input_10x_dir}/genes.tsv
        Rscript workflow/scripts/prepare_rna_counts.R \
            --input_10x_dir {params.input_10x_dir} \
            --barcodes_nuclei {input.barcodes_nuclei} \
            --barcodes_empty {input.barcodes_empty} \
            --counts_nuclei {output.counts_nuclei} \
            --counts_empty {output.counts_empty}
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

rule doublet_filtering:
    input:
        seurat_obj = _results("seurat_prelim/{sample}/seurat_obj.rds"),
        doublets = _results("amulet/{sample}/MultipletBarcodes_01.txt"),
        map = _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
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
            --map {input.map} \
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
