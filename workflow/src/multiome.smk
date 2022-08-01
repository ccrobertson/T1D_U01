#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import pandas as pd
import math

_results = partial(os.path.join, "results/multiome")
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


configfile: _results("Sample_islet_list_for_multiomics_Batches_long_format_with_libname.json")
samples = config["samples"]
ambient_gene_thresholds = [0.01, 0.05, 0.1, 0.2, 0.3]


IONICE = 'ionice -c2 -n7'

def iterate_subjects_by_batch(samplename):
     return config["batches"][samplename].keys()

def count_subjects_by_batch(samplename):
     return len(config["batches"][samplename].keys())

def iterate_demuxlet_samples(demuxlet_report):
    d = pd.read_csv(demuxlet_report)
    samples = d["Sample_ID"].tolist()
    return samples

def calculate_max_missing_by_batch(samplename, genotyped_donors):
    n_subjects = len(config["batches"][samplename].keys())
    n_donors = pd.read_table(genotyped_donors, header=None).shape[0]
    max_missing_rate = math.ceil((n_subjects-1)/(n_donors+n_subjects)*100, )/100
    return max_missing_rate



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
        # post decontx clustering
        #expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample=samples),
        #~~~~~~~~ souporcell
        #expand(_results("souporcell_gex_{p}/{sample}/cluster_genotypes_reformatted.vcf.gz"), sample="Sample_5124-NM-2-hg38",p=ambient_gene_thresholds),
        #expand(_results("souporcell_atac/{sample}/cluster_genotypes_reformatted.vcf.gz"), sample="Sample_5124-NM-2-hg38"),
        #~~~~~~~~ demultiplexing
        #expand(_results("souporcell_gex_{p}/{sample}/corrplot_kinship_cat.png"), sample=samples, p=ambient_gene_thresholds),
        #expand(_results("souporcell_atac/{sample}/corrplot_kinship_cat.png"), sample=samples),
        expand(_results("demultiplex_gex_{p}/corrplot_r_factor_by_donor.png"), sample=samples, p=ambient_gene_thresholds),
        expand(_results("demultiplex_atac/corrplot_r_factor_by_donor.png"), sample=samples),
        #~~~~~~~ get tracks
        #expand(_results("tracks_by_sample/gex.{sample}.TPM.bw"), sample=samples),


rule symlinks_atac:
    input:
        demuxlet_report_atac = config["demuxlet_report_atac"],
    output:
        json = _results("nf_atac_config.json"),
    params:
        indir = config["fastq_dir_atac"],
        outdir = _results("fastq_atac"),
        samples = iterate_demuxlet_samples(config["demuxlet_report_atac"]),
    shell:
        """
        createSymLinksATAC () {{
          local target_dir=$1
          local link_dir=$2
          local sample=$3
          ln -s --force $target_dir/$sample/*R1_001.fastq.gz $link_dir/${{sample}}_R1_001.fastq.gz
          ln -s --force $target_dir/$sample/*R2_001.fastq.gz $link_dir/${{sample}}_R2_001.fastq.gz
          ln -s --force $target_dir/$sample/*R3_001.fastq.gz $link_dir/${{sample}}_R3_001.fastq.gz
        }}

        mkdir -p {params.outdir}
        for sample in {params.samples}; do
            echo $sample
            createSymLinksATAC {params.indir} {params.outdir} $sample
        done

        python workflow/scripts/build_multiome_json.py --demuxlet_report {input.demuxlet_report} --fastq_dir {params.outdir} --modality ATAC > {output.json}
        """


rule symlinks_gex:
    input:
        demuxlet_report = config["demuxlet_report_gex"],
    output:
        json = _results("nf_gex_config.json"),
    params:
        indir = config["fastq_dir_gex"],
        outdir = _results("fastq_gex"),
        samples = iterate_demuxlet_samples(config["demuxlet_report_gex"]),
    shell:
        """
        createSymLinksGEX () {{
          local target_dir=$1
          local link_dir=$2
          local sample=$3
          ln -s --force $target_dir/$sample/*R1_001.fastq.gz $link_dir/${{sample}}_R1_001.fastq.gz
          ln -s --force $target_dir/$sample/*R2_001.fastq.gz $link_dir/${{sample}}_R2_001.fastq.gz
        }}

        mkdir -p {params.outdir}
        for sample in {params.samples}; do
            echo $sample
            createSymLinksGEX {params.indir} {params.outdir} $sample
        done

        python workflow/scripts/build_multiome_json.py --demuxlet_report {input.demuxlet_report} --fastq_dir {params.outdir} --modality GEX > {output.json}
        """

## INSERT RULE FOR RUNNING ATAC NEXTFLOW
#rule nf_atac:

## INSERT RULE FOR RUNNING GEX NEXTFLOW
#rule_nf_gex:

## NEED TO UPDATE THIS TO PROVIDE SEPARATE BIGWIGS FOR EACH STRAND
rule bam_to_bigwig:
    input:
        bam = _results("nf_gex_results/prune/{sample}.before-dedup.bam"),
        blacklist = config["blacklist"],
    output:
        bw =  _results("tracks_by_sample/gex.{sample}.TPM.bw"),
    conda:
        "general"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} --normalizeUsing BPM --numberOfProcessors 8 --blackListFileName {input.blacklist}
        """


#Some downstream programs require gzipped output, while others need it unzipped
#so we will keep it available in both forms
rule gzip_starsolo:
    input:
        mtx = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/matrix.mtx"),
        barcodes = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/barcodes.tsv"),
        features = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/features.tsv"),
    output:
        mtx_gz = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/matrix.mtx.gz"),
        barcodes_gz = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/barcodes.tsv.gz"),
        features_gz = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/features.tsv.gz"),
    shell:
        """
        gzip -c {input.mtx} > {output.mtx_gz}
        gzip -c {input.barcodes} > {output.barcodes_gz}
        gzip -c {input.features} > {output.features_gz}
        """


rule dropkick:
    input:
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx.gz"),
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/barcodes.tsv.gz"),
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv.gz"),
    output:
        _results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"),
        _results("dropkick_{method}_{min_genes}/{sample}/dropkick.csv"),
    conda:
        "dropkick"
        # NOTE: Snakemake cannot build dropkick environment for envs/dropkick.yaml
        # because dropkick is not available via pip
        #library://porchard/default/dropkick:20220225
    params:
        indir = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw"),
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

rule cross_modality_prefiltering:
    input:
        gex = _results("nf_gex_results/qc/{sample}.qc.txt"),
        atac = _results("nf_atac_results/ataqv/single-nucleus/{sample}.txt"),
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
        bam = _results("nf_atac_results/prune/{sample}.pruned.bam"),
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
        mtx = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx.gz"),
        barcodes = _results("cross_modality_qc/{sample}/barcodes_prefiltered.txt"),
        qcstats = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.txt"),
    output:
        h5ad_file = _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.h5ad"),
    params:
        starsolodir = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw"),
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



rule prepare_rna_counts:
    input:
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv"),
        barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
        barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
    output:
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/genes.tsv"),
        counts_nuclei = _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
    params:
        input_10x_dir = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw"),
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


rule filter_bams_for_nuclei:
    input:
        gex_bam = _results("nf_gex_results/prune/{sample}.before-dedup.bam"),
        atac_bam = _results("nf_atac_results/prune/{sample}.pruned.bam"),
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
        gtf = _resources("hg38/gencode.v39.annotation.gtf"),
    output:
        regions = _results("bam_pass_qc_barcodes/{sample}/ambient_gene_regions_{p}.bed"),
    params:
        outdir = _results("bam_pass_qc_barcodes/{sample}"),
        p = "{p}",
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/define_ambient_gene_regions.R \
            --counts_empty {input.counts_empty} \
            --gtf {input.gtf} \
            --p {params.p} \
            --outdir {params.outdir}
        """

rule filter_bams_for_ambient_genes:
    input:
        gex_bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex.bam"),
        regions = _results("bam_pass_qc_barcodes/{sample}/ambient_gene_regions_{p}.bed"),
    output:
        gex_bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex_exclude_ambient_genes_{p}.bam"),
        gex_bam_discarded = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex_discarded_{p}.bam"),
    container:
        "workflow/envs/arushi_general.simg"
    shell:
        """
        samtools view -h -b -@ 10 -L {input.regions} -U {output.gex_bam} {input.gex_bam} > {output.gex_bam_discarded}
        samtools index {output.gex_bam}
        """

rule souporcell_gex:
    input:
        bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex_exclude_ambient_genes_{p}.bam"),
        barcodes=_resources("barcode_whitelist_multiome_GEX_cp.txt"),
        fasta=_resources("hg38/hg38_cvb4.fa"),
        common_variants=_resources("common_variants_grch38_fixed.vcf"),
        donor_genotypes = config["donor_genotypes"],
    output:
        #clusters = _results("souporcell_gex_{p}/{sample}/clusters.tsv"),
        cluster_genotypes = _results("souporcell_gex_{p}/{sample}/cluster_genotypes.vcf"),
        #cluster_genotypes_reformatted = _results("souporcell_gex_{p}/{sample}/cluster_genotypes_reformatted.vcf.gz"),
    params:
        outdir = _results("souporcell_gex_{p}/{sample}"),
        k = lambda wildcards: count_subjects_by_batch(wildcards.sample),
        threads=10,
        sample = "{sample}",
    shell:
        """
        rm -f {params.outdir}/clustering.done
        rm -f {params.outdir}/consensus.done
        rm -f {params.outdir}/troublet.done

        singularity exec workflow/envs/souporcell_latest.sif souporcell_pipeline.py \
            -i {input.bam} \
            -b {input.barcodes} \
            -f {input.fasta} \
            -t {params.threads} \
            --cluster {params.k} \
            --common_variants {input.common_variants} \
            --skip_remap True \
            -o {params.outdir}
        """

## NOTE: --ignore True and --no_umi option tells souporcell to ignore UMIs
## (since this is ATAC-seq, there are no UMIs)
rule souporcell_atac:
    input:
        bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_atac.bam"),
        barcodes = _resources("barcode_whitelist_multiome_ATAC_cp.txt"),
        fasta = _resources("hg38/hg38_cvb4.fa"),
        common_variants = _resources("common_variants_grch38_fixed.vcf"),
        donor_genotypes = config["donor_genotypes"],
    output:
        #clusters = _results("souporcell_atac/{sample}/clusters.tsv"),
        cluster_genotypes = _results("souporcell_atac/{sample}/cluster_genotypes.vcf"),
        #cluster_genotypes_reformatted = _results("souporcell_atac/{sample}/cluster_genotypes_reformatted.vcf.gz"),
    params:
        outdir = _results("souporcell_atac/{sample}"),
        k = lambda wildcards: count_subjects_by_batch(wildcards.sample),
        threads=10,
        sample = "{sample}",
    shell:
        """
        rm -f {params.outdir}/clustering.done
        rm -f {params.outdir}/consensus.done
        rm -f {params.outdir}/troublet.done

        singularity exec workflow/envs/souporcell_latest.sif souporcell_pipeline.py \
            -i {input.bam} \
            -b {input.barcodes} \
            -f {input.fasta} \
            --ignore True \
            --no_umi True \
            -t {params.threads} \
            --cluster {params.k} \
            --common_variants {input.common_variants} \
            --skip_remap True \
            -o {params.outdir}
        """

#NEXT - explore a range of k
#https://github.com/wheaton5/souporcell/issues/7

#I think this may be failing because one of the subjects in our pool is missing from the genotype VCF
#https://github.com/wheaton5/souporcell/issues/77
#https://github.com/wheaton5/souporcell/issues/141
# rule souporcell_known_genotypes:
#     input:
#         bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex_exclude_ambient_genes.bam"),
#         barcodes = _resources("barcode_whitelist_multiome_GEX_cp.txt"),
#         fasta = _resources("hg38/hg38_cvb4.fa"),
#         common_variants = _resources("common_variants_grch38_fixed.vcf"),
#         donor_genotypes = config["donor_genotypes"],
#     output:
#         clusters = _results("souporcell_known_genotypes/{sample}/clusters.tsv"),
#         cluster_genotypes = _results("souporcell_known_genotypes/{sample}/cluster_genotypes.vcf"),
#     params:
#         outdir = _results("souporcell_known_genotypes/{sample}"),
#         k = lambda wildcards: count_subjects_by_batch(wildcards.sample),
#         threads=10,
#         sample = "{sample}",
#     shell:
#         """
#         singularity exec workflow/envs/souporcell_latest.sif souporcell_pipeline.py \
#             -i {input.bam} \
#             -b {input.barcodes} \
#             -f {input.fasta} \
#             -t {params.threads} \
#             --cluster {params.k} \
#             --known_genotypes {input.donor_genotypes} \
#             --known_genotypes_sample_names HPAP093 HPAP105 ICRH139 ICRH142 HPAP107 \
#             --skip_remap True \
#             -o {params.outdir}
#         """


rule souporcell_fix_vcf:
    input:
        cluster_genotypes = _results("souporcell_{input}/{sample}/cluster_genotypes.vcf"),
    output:
        cluster_genotypes_reformatted = _results("souporcell_{input}/{sample}/cluster_genotypes_reformatted.vcf.gz"),
    params:
        sample = "{sample}",
    conda:
        "genetics"
    shell:
        """
        awk -v OFS='\t' -v sample={params.sample} '$1!~/^#CHROM/ {{print $0}} $1~/^#CHROM/ {{for(i=10; i<=NF; ++i) $i=sample"_souporcell_"$i; print $0 }}' {input.cluster_genotypes} | grep -v BACKGROUND | bgzip > {output.cluster_genotypes_reformatted}
        tabix -p vcf {output.cluster_genotypes_reformatted}
        """

rule demultiplex_by_sample:
    input:
        souporcell_vcf = _results("souporcell_{input}/{sample}/cluster_genotypes_reformatted.vcf.gz"),
        imputed_vcf = "results/imputation/2022_02_16_T1D_genotypes/imputation_results/chrALL.donors_only.maf_gt_0.01__rsq_gt_0.95.dose.vcf.gz",
    output:
        kin0 =  _results("souporcell_{input}/{sample}/souporcell_clusters_and_donors_filtered.kin0"),
    params:
        prefix = _results("souporcell_{input}/{sample}/souporcell_clusters_and_donors"),
        prefix_filtered = _results("souporcell_{input}/{sample}/souporcell_clusters_and_donors_filtered"),
        max_missing = lambda wildcards: calculate_max_missing_by_batch(wildcards.sample, config["genotyped_donors"])
    conda:
        "genetics"
    shell:
        """
        #merge souporcell and donor vcfs
        bcftools merge -O z -o {params.prefix}.vcf.gz {input.souporcell_vcf} {input.imputed_vcf}
        tabix -p vcf {params.prefix}.vcf.gz

        #convert merged vcf to plink
        plink --vcf {params.prefix}.vcf.gz --double-id --make-bed --out {params.prefix}

        #filter merged vcf for missingness
        plink --bfile {params.prefix} --geno {params.max_missing} --make-bed --out {params.prefix_filtered}

        #run king relationship inference
        plink2 --bfile {params.prefix_filtered} --make-king-table --out {params.prefix_filtered}
        """


rule demultiplex_by_sample_plots:
    input:
        kin0 = _results("souporcell_{input}/{sample}/souporcell_clusters_and_donors_filtered.kin0"),
    output:
        png = _results("souporcell_{input}/{sample}/corrplot_kinship_cat.png"),
    params:
        outdir = _results("souporcell_{input}/{sample}"),
        donors = "HPAP105,HPAP093,ICRH139,ICRH142,HPAP107",
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/demultiplexing_plots.R --kin0 {input.kin0} --outdir {params.outdir} --donors {params.donors}
        """


rule demultiplex:
    input:
        souporcell_vcfs = expand(_results("souporcell_{{input}}/{sample}/cluster_genotypes_reformatted.vcf.gz"), sample=samples),
        imputed_vcf = "results/imputation/2022_02_16_T1D_genotypes/imputation_results/chrALL.donors_only.maf_gt_0.01__rsq_gt_0.95.dose.vcf.gz",
    output:
        kin0 =  _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat =  _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered.king"),
        kinIds =  _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered.kin.id"),
    params:
        prefix = _results("demultiplex_{input}/souporcell_clusters_and_donors"),
        prefix_filtered = _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered"),
        max_missing = 0.5
    shell:
        """
        #merge souporcell and donor vcfs
        bcftools merge -O z -o {params.prefix}.vcf.gz {input.souporcell_vcfs} {input.imputed_vcf}
        tabix -p vcf {params.prefix}.vcf.gz

        #convert merged vcf to plink
        plink --vcf {params.prefix}.vcf.gz --double-id --make-bed --out {params.prefix}

        #filter merged vcf for missingness
        plink --bfile {params.prefix} --geno {params.max_missing} --make-bed --out {params.prefix_filtered}

        #run king relationship inference
        plink2 --bfile {params.prefix_filtered} --make-king square --make-king-table --out {params.prefix_filtered}

        #recode as alt allele count matrix
        plink --bfile {params.prefix_filtered} --recodeA --out {params.prefix_filtered}
        """


rule demultiplex_plots:
    input:
        kin0 = _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("demultiplex_{input}/souporcell_clusters_and_donors_filtered.king.id"),
    output:
        _results("demultiplex_{input}/corrplot_kinship_factor_by_donor.png"),
        _results("demultiplex_{input}/corrplot_r_factor_by_donor.png"),
    params:
        outdir = _results("demultiplex_{input}"),
        donors = "HPAP105,HPAP093,ICRH139,ICRH142,HPAP107",
    shell:
        """
        Rscript workflow/scripts/demultiplexing_scratch.R \
            --kin0 {input.kin0} \
            --kinMat {input.kinMat} \
            --kinIds {input.kinIds} \
            --donors "HPAP105,HPAP093,ICRH139,ICRH142,HPAP107" \
            --outdir {params.outdir}
        """
