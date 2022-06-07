#!/usr/bin/env python3

from os.path import join
import os
from functools import partial


_results = partial(os.path.join, "results/multiome")
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


configfile: _results("Sample_islet_list_for_multiomics_Batches_long_format_with_libname.json")
samples = config["samples"]

IONICE = 'ionice -c2 -n7'

def iterate_subjects_by_batch(samplename):
     return config["batches"][samplename].keys()

def count_subjects_by_batch(samplename):
     return len(config["batches"][samplename].keys())

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
        #~~~~~~~~ select_quality_nuclei
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.txt"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("cross_modality_qc/{sample}/cross_modality_qc.txt"), sample=samples),
        #~~~~~~~~ amulet
        #expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample="Sample_5124-NM-2-hg38"),
        expand(_results("amulet/{sample}/MultipletBarcodes_01.txt"), sample=samples),
        #~~~~~~~~ prepare_rna_data
        #expand(_results("decontx/{sample}/counts_nuclei.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("decontx/{sample}/counts_nuclei.rds"), sample=samples),
        #~~~~~~~~ decontx
        #expand(_results("decontx/{sample}/counts_nuclei_decontaminated.rds"), sample="Sample_5124-NM-2-hg38"),
        #expand(_results("decontx/{sample}/counts_nuclei_decontaminated.rds"), sample=samples),
        #~~~~~~~~ souporcell
        #expand(_results("souporcell/{sample}/clusters.tsv"), sample=samples),
        #~~~~~~~~ seurat
        #expand(_results("seurat_prelim/{sample}/seurat_obj.rds"), sample=samples),
        #expand(_results("seurat_round2/{sample}/seurat_obj.rds"), sample=samples)

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
    shell:
        """
        Rscript workflow/scripts/cross_modality_qc2.R \
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


# rule anndata:
#     input:
#         _results("cross_modality_qc/{sample}/cross_modality_qc_prefiltered.txt")
#     output:
#     shell:
#
#
# rule wtsi_hgi_qc:
#     input:
#     output:
#     shell:
#         """
#         singularity exec workflow/envs/nf_qc_cluster.sif workflow/scripts/filter_outlier_cells.py
#         """

# rule cross_modality_qc_plots:
#     input:
#         gex = _results("nf_gex_results/qc/{sample}.qc.txt"),
#         atac = _results("nf_atac_results/ataqv/single-nucleus/{sample}.txt"),
#         #dropkick = _results("dropkick_otsu_50/{sample}/dropkick.csv"),
#     output:
#         # metrics = _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
#         # barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
#         # barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
#         _results("cross_modality_qc/{sample}/cross_modality_qc.txt"),
#         #_results("cross_modality_qc/{sample}/cross_modality_qc_plot_status_dropkick.png"),
#         _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
#         _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
#         #amulet = _results("amulet/{sample}/amulet_input.csv"),
#     params:
#         outdir = _results("cross_modality_qc/{sample}"),
#     shell:
#         """
#         Rscript workflow/scripts/cross_modality_qc.R \
#             --input_gex {input.gex} \
#             --input_atac {input.atac} \
#             --outdir {params.outdir} \
#             --sample {wildcards.sample} \
#             --min_genes {params.min_genes}
#         """





# rule prepare_rna_counts:
#     input:
#         _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
#         _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv"),
#         barcodes_nuclei = _results("cross_modality_qc/{sample}/barcodes_nuclei.txt"),
#         barcodes_empty = _results("cross_modality_qc/{sample}/barcodes_empty.txt"),
#         #barcodes_doublets = _results("amulet/{sample}/MultipletBarcodes_01.txt"),
#     output:
#         _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/genes.tsv"),
#         counts_nuclei = _results("decontx/{sample}/counts_nuclei.rds"),
#         counts_empty = _results("decontx/{sample}/counts_empty.rds"),
#     params:
#         input_10x_dir = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw"),
#     conda:
#         "Renv"
#     shell:
#         """
#         ln -s --relative --force {params.input_10x_dir}/features.tsv {params.input_10x_dir}/genes.tsv
#         Rscript workflow/scripts/prepare_rna_counts.R \
#             --input_10x_dir {params.input_10x_dir} \
#             --barcodes_nuclei {input.barcodes_nuclei} \
#             --barcodes_empty {input.barcodes_empty} \
#             --counts_nuclei {output.counts_nuclei} \
#             --counts_empty {output.counts_empty}
#         """
#
#
# rule prelim_clustering:
#     input:
#         counts = _results("decontx/{sample}/counts_nuclei.rds"),
#     output:
#         clusters = _results("seurat_prelim/{sample}/seurat_obj.rds"),
#     params:
#         outdir = _results("seurat_prelim/{sample}"),
#     shell:
#         """
#         Rscript workflow/scripts/run_seurat_prelim.R --counts {input.counts} --outdir {params.outdir}
#         """
#
# rule decontx:
#     input:
#         counts_nuclei = _results("decontx/{sample}/counts_nuclei.rds"),
#         counts_empty = _results("decontx/{sample}/counts_empty.rds"),
#     output:
#         _results("decontx/{sample}/counts_nuclei_decontaminated.rds"),
#     conda:
#         "Renv"
#     params:
#         outdir = _results("decontx/{sample}"),
#     shell:
#         """
#         Rscript workflow/scripts/run_decontx.R \
#             --counts_nuclei {input.counts_nuclei} \
#             --counts_empty {input.counts_empty} \
#             --outdir {params.outdir}
#         """
#
#
# rule second_clustering:
#     input:
#         counts = _results("decontx/{sample}/counts_nuclei_decontaminated.rds"),
#     output:
#         clusters = _results("seurat_round2/{sample}/seurat_obj.rds"),
#     params:
#         outdir = _results("seurat_round2/{sample}"),
#     shell:
#         """
#         Rscript workflow/scripts/run_seurat_prelim.R --counts {input.counts} --outdir {params.outdir}
#         """
#
# # rule liger:
# #     input:
# #         counts = _results("decontx/{sample}/counts_nuclei_decontaminated.rds"),
# #     output:
# #     shell:
# #         """
# #         Rscript run_liger.R \
# #             --counts {input.counts} \
# #         """
#
#
# rule souporcell:
#     input:
#         bam=_results("nf_gex_results/prune/{sample}.before-dedup.bam"),
#         barcodes=_resources("barcode_whitelist_multiome_GEX_cp.txt"),
#         fasta=_resources("hg38/hg38_cvb4.fa"),
#         common_variants=_resources("common_variants_grch38_fixed.vcf"),
#     output:
#         _results("souporcell/{sample}/clusters.tsv"),
#     params:
#         outdir = _results("souporcell/{sample}"),
#         k = lambda wildcards: count_subjects_by_batch(wildcards.sample),
#         threads=10,
#     shell:
#         """
#           singularity exec workflow/envs/souporcell_latest.sif workflow/scripts/souporcell_pipeline.py \
#             -i {input.bam} \
#             -b {input.barcodes} \
#             -f {input.fasta} \
#             -t {params.threads} \
#             --cluster {params.k} \
#             --common_variants {input.common_variants} \
#             --skip_remap True \
#             -o {params.outdir}
#         """
