#!/usr/bin/env python3

from os.path import join
import os
from functools import partial


_results = partial(os.path.join, config["results"])
#_resources = partial(os.path.join, config["resources"])
_logs = partial(_results, "logs")
samples = config["samples"]

IONICE = 'ionice -c2 -n7'

rule all:
    input:
        #~~~~~~~~ gzip starsolo
        #expand(_results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/matrix.mtx.gz"),sample="Sample_test-hg38", feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        #expand(_results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/{feature_type}/raw/matrix.mtx.gz"), sample=samples, feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        #~~~~~~~~ run dropkick
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"), sample="Sample_5124-NM-2-hg38", method=["otsu", "multiotsu"], min_genes=[30, 50]),
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick.h5ad"), sample=samples, method=["otsu", "multiotsu"], min_genes=[30, 50]),
        #~~~~~~~~ dropkick_plots
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick_score.png"), sample="Sample_5124-NM-2-hg38", method=["otsu", "multiotsu"], min_genes=[30, 50]),
        #expand(_results("dropkick_{method}_{min_genes}/{sample}/dropkick_score.png"), sample=samples, method=["otsu", "multiotsu"], min_genes=[30, 50]),
        #~~~~~~~~ select_quality_nuclei
        expand(_results("cross_modality_qc/{sample}/cross_modality_qc.txt"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.txt"), sample=samples),
        #~~~~~~~~ amulet
        expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample=samples),
        #~~~~~~~~ prepare_rna_data
        expand(_results("decontx/{sample}/counts_nuclei.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("decontx/{sample}/counts_nuclei.rds"), sample=samples),
        #~~~~~~~~ decontx
        expand(_results("decontx/{sample}/counts_nuclei_decontaminated.rds"), sample="Sample_5124-NM-2-hg38"),


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


rule select_quality_nuclei:
    input:
        gex = _results("nf_gex_results/qc/{sample}.qc.txt"),
        atac = _results("nf_atac_results/ataqv/single-nucleus/{sample}.txt"),
        dropkick = _results("dropkick_otsu_50/{sample}/dropkick.csv"),
    output:
        # metrics = _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
        # barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
        # barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
        _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
        _results("cross_modality_qc/{sample}/cross_modality_qc_plot_status_dropkick.png"),
        _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
        _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
        amulet = _results("amulet/{sample}/amulet_input.csv"),
    params:
        outdir = _results("cross_modality_qc/{sample}"),
        barcode_filtering_config = "workflow/scripts/cross_modality_qc.yaml",
    shell:
        """
        Rscript workflow/scripts/cross_modality_qc.R \
            --input_gex {input.gex} \
            --input_atac {input.atac} \
            --input_dropkick {input.dropkick} \
            --outdir {params.outdir} \
            --output_amulet {output.amulet} \
            --sample {wildcards.sample} \
            --barcode_filtering_config {params.barcode_filtering_config}
        """


rule doublet_detection:
    input:
        bam = _results("nf_atac_results/prune/{sample}.pruned.bam"),
        csv = _results("amulet/{sample}/amulet_input.csv"),
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


rule prepare_rna_counts:
    input:
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv"),
        barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
        barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
        barcodes_doublets = _results("amulet/{sample}/MultipletBarcodes_01.txt"),
    output:
        _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/genes.tsv"),
        counts_nuclei = _results("decontx/{sample}/counts_nuclei.rds"),
        counts_empty = _results("decontx/{sample}/counts_empty.rds"),
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
            --barcodes_doublets {input.barcodes_doublets} \
            --counts_nuclei {output.counts_nuclei} \
            --counts_empty {output.counts_empty}
        """


rule decontx:
    input:
        counts_nuclei = _results("decontx/{sample}/counts_nuclei.rds"),
        counts_empty = _results("decontx/{sample}/counts_empty.rds"),
    output:
        _results("decontx/{sample}/counts_nuclei_decontaminated.rds"),
    conda:
        "Renv"
    params:
        outdir = _results("decontx/{sample}"),
    shell:
        """
        Rscript workflow/scripts/run_decontx.R \
            --counts_nuclei {input.counts_nuclei} \
            --counts_empty {input.counts_empty} \
            --outdir {params.outdir}
        """

rule liger:
    input:
        _results("decontx/{sample}/counts_nuclei_decontaminated.rds"),
    output:
    shell:
        """
        Rscript run_liger.R
        """









# rule convert_to_h5ad:
#     input:
#         dir = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/Gene/raw")
#     output:
#         h5ad = _results("dropkick/{sample}/gex.h5ad"),
#         h5Seurat = _results("dropkick/{sample}/gex.h5Seurat"),
#     conda:
#         "envs/Renv.yaml"
#     shell:
#         """
#         Rscript workflow/scripts/convert_to_h5ad.R --datadir {input.dir} --outad {output.h5ad} --outSeurat {output.h5Seurat}
#         """
# # NOTE: Snakemake cannot build dropkick environment for envs/dropkick.yaml
# # because dropkick is not available via pip
# rule dropkick:
#     input:
#         counts = _results("dropkick/{sample}/gex.h5ad")
#     output:
#         _results("dropkick/{sample}/dropkick.h5ad"),
#         _results("dropkick/{sample}/dropkick.csv"),
#     conda:
#         "dropkick"
#     params:
#         outdir = _results("dropkick/{sample}"),
#         cpus = 2,
#         min_genes = 50,
#     shell:
#         """
#         python -u workflow/scripts/run_dropkick.py --counts {input.counts} --outdir {params.outdir} --cpus {params.cpus} --min_genes {params.min_genes}
#         """
