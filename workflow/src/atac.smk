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



rule all:
    input:
        #~~~~~~~~ gzip starsolo
        #expand(_results("starsolo/{sample}/starsolo.Solo.out/{feature_type}/raw/genes.tsv"),sample="Sample_5125-NM-1_NM-1", feature_type=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron']),
        #~~~~~~~~ qc
        #expand(_results("qc/{sample}.qc.txt"), sample="Sample_5125-NM-1_NM-1"),
        expand(_results("qc/{sample}.qc.metrics.png"), sample=samples),
        #~~~~~~~~ seurat_round3
        expand(_results("seurat_round3/{sample}/seurat_obj.rds"), sample=samples),


## Approach for generating gene count matrix form snATAC from LIGER authors
# http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html
rule prep_count_matrix:


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
