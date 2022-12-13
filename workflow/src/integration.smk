#!/usr/bin/env python3

from os.path import join
import os
import pandas as pd
from functools import partial
import math
import glob


_results = partial(os.path.join, "results")
_data = partial(os.path.join, "data")
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


def get_samples_for_set(set):
    samples = config["sets"][set]
    return samples

rule all:
    input:
        #~~~~~ get barcode map
        expand("results/master_barcode_map_{set}.rds", set="multiome_5GEX"),
        #expand("results/master_barcode_map_{set}.rds", set="3GEX_all"),
        #~~~~~ get barcode lists
        #expand(_results("atac_barcodes_sample_cluster/{set}/atac_barcodes_{sample}_{cluster}.txt"), sample=get_samples_for_set("multiome_5GEX"), cluster=get_clusters_for_set("multiome_5GEX")),
        #expand(_results("atac_barcodes_sample_cluster/all/atac_barcodes_{sample}_{cluster}.txt"), sample=samples_all, cluster=cell_types_all),
        #~~~~~ filter bams
        #expand(_results("atac_bam_sample_cluster/noNM3/atac_{sample}_{cluster}.bam"),sample=samples_noNM3, cluster=cell_types_noNM3),
        #expand(_results("atac_bam_sample_cluster/all/atac_{sample}_{cluster}.bam"),sample=samples_all, cluster=cell_types_all),
        #~~~~~ merge bams
        #expand(_results("atac_bam_cluster/noNM3/atac_{cluster}.bam"), cluster=cell_types_noNM3),
        #expand(_results("atac_bam_cluster/all/atac_{cluster}.bam"), cluster=cell_types_all),


### NOTE: NEED TO INCLUDE BOTH MULTIOME BARCODES ATAC+GEX
rule master_barcode_map:
    input:
        barcode_to_cell_type = "results/liger/{set}/barcode_to_cluster_to_cell_type.csv",
        run_design = _data("batch_design.tsv"),
        demuxlet_files = lambda wildcards: expand("results/demultiplex/demuxlet-unfiltered/{sample}/demuxlet.best", sample = get_samples_for_set(wildcards.set))
    output:
        "results/master_barcode_maps/master_barcode_map_{set}.rds"
    params:
        prefix = "results/master_barcode_maps/master_barcode_map_{set}",
        PSNG_threshold = 0.8,
    shell:
        """
        Rscript workflow/scripts/master_barcode_map.R \
            --barcode_to_cell_type {input.barcode_to_cell_type} \
            --run_design {input.run_design} \
            --PSNG_threshold {params.PSNG_threshold} \
            --prefix {params.prefix} \
            {input.demuxlet_files}
        """

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   For each modality - get bams by cluster/condition
#
# NOTE: THIS IS A WORK IN PROGRESS -- CURRENTLY NOT FUNCTIONAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule barcodes_by_sample_by_cluster:
    input:
        barcode_to_cell_type = "results/liger/{set}/barcode_to_cluster_to_cell_type.csv",
    output:
        barcodes = _results("atac_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt"),
    params:
        sample = "{sample}",
        cluster = "{cluster}",
    shell:
        """
        sample=$(echo "{params.sample}" | sed 's/-/_/g')
        awk -v s="$sample" -v c="{params.cluster}" 'BEGIN {{FS=","}} $1==s && $4==c {{ print $2"-1"}}' {input.barcode_to_cell_type} > {output.barcodes}
        """

rule bam_by_sample_by_cluster:
    input:
        bam = lambda wildcards: get_bam(wildcards.set, wildcards.sample),
        barcodes = _results("atac_barcodes_sample_cluster/{set}/atac_barcodes_{sample}_{cluster}.txt"),
    output:
        bam = _results("atac_bam_sample_cluster/{set}/atac_{sample}_{cluster}.bam"),
        bai = _results("atac_bam_sample_cluster/{set}/atac_{sample}_{cluster}.bam.bai"),
    container:
        "workflow/envs/arushi_general.simg"
    shell:
        """
        samtools view -h -b -@ 10 -D CB:{input.barcodes} {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule merge_bams:
    input:
        bams = lambda wildcards: expand(_results("atac_bam_sample_cluster/{{set}}/atac_{sample}_{{cluster}}.bam"), sample=get_samples_for_set(wildcards.set))
    output:
        bam = _results("atac_bam_cluster/{set}/atac_{cluster}.bam"),
        bai = _results("atac_bam_cluster/{set}/atac_{cluster}.bam.bai"),
    container:
        "workflow/envs/arushi_general.simg"
    shell:
        """
    	samtools merge -@ 15 -O BAM {output.bam} {input.bams}
        samtools index {output.bam}
        """
