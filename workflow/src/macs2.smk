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
configfile: "workflow/src/markers.yaml"

def get_bam(set, sample):
    bam = config["sets"][set]["samples"][sample]
    return bam

def get_samples_for_set(set):
    samples = list(config["sets"][set]["samples"].keys())
    return samples

def get_clusters_for_set(set):
    d = pd.read_csv("results/liger/3GEX_"+set+"/cluster_to_cell_type.csv")
    return d['cell_type'].unique().tolist()

def barcode_files_for_set(set):
    files = expand(
        _results("atac_barcodes_sample_cluster"+set+"/atac_barcodes_{sample}_{cluster}.txt"),
            sample=get_samples_for_set(set),
            cluster=get_clusters_for_set(set)
        )
    return files

cell_types_all = get_clusters_for_set("all")
cell_types_noNM3 = get_clusters_for_set("noNM3")
samples_all = get_samples_for_set("all")
samples_noNM3 = get_samples_for_set("noNM3")


rule all:
    input:
        #~~~~~ get barcode lists
        #expand(_results("atac_barcodes_sample_cluster/noNM3/atac_barcodes_{sample}_{cluster}.txt"), sample=samples_noNM3, cluster=cell_types_noNM3),
        #expand(_results("atac_barcodes_sample_cluster/all/atac_barcodes_{sample}_{cluster}.txt"), sample=samples_all, cluster=cell_types_all),
        #~~~~~ filter bams
        #expand(_results("atac_bam_sample_cluster/noNM3/atac_{sample}_{cluster}.bam"),sample=samples_noNM3, cluster=cell_types_noNM3),
        #expand(_results("atac_bam_sample_cluster/all/atac_{sample}_{cluster}.bam"),sample=samples_all, cluster=cell_types_all),
        #~~~~~ merge bams
        #expand(_results("atac_bam_cluster/noNM3/atac_{cluster}.bam"), cluster=cell_types_noNM3),
        #expand(_results("atac_bam_cluster/all/atac_{cluster}.bam"), cluster=cell_types_all),
        #~~~~~ call peaks
        expand(_results("peakcalls/noNM3/{cluster}_summits_ext500_noblacklist.bed"), cluster=cell_types_noNM3),
        #expand(_results("peakcalls/{set}/{cluster}_summits.bed"), cluster=get_clusters_for_set("all"), set="all"),
        #expand(_results("peakcalls/{set}/{cluster}_summits.bed"), cluster=get_clusters_for_set("noNM3"), set="noNM3"),

#want a file with the following columns:
#   sample barcode donor condition cluster cell_type
# rule barcode_map:
#     input:
#         demuxlet_best
#         liger_clusters = "results/liger/liger_clusters.tsv")
#     output:
#         """
#         """


rule barcodes_by_sample_by_cluster_noNM3:
    input:
        barcode_to_cell_type = "results/liger/3GEX_noNM3/barcode_to_cluster_to_cell_type.csv",
        barcode_map = "resources/multiome_barcode_map.rds",
    output:
        expand(_results("atac_barcodes_sample_cluster/noNM3/atac_barcodes_{sample}_{cluster}.txt"), sample=samples_noNM3, cluster=cell_types_noNM3),
        #lambda wildcards: barcode_files_for_set(wildcards.set),
    params:
        samples = samples_noNM3,
        cell_types = cell_types_noNM3,
        outdir = _results("atac_barcodes_sample_cluster/noNM3"),
    shell:
        """
        Rscript workflow/scripts/barcodes_by_sample_by_cluster.R \
            --barcode_to_cell_type {input.barcode_to_cell_type} \
            --barcode_map {input.barcode_map} \
            --outdir {params.outdir}
        """

rule barcodes_by_sample_by_cluster_all:
    input:
        barcode_to_cell_type = "results/liger/3GEX_all/barcode_to_cluster_to_cell_type.csv",
        barcode_map = "resources/multiome_barcode_map.rds",
    output:
        expand(_results("atac_barcodes_sample_cluster/all/atac_barcodes_{sample}_{cluster}.txt"), sample=samples_all, cluster=cell_types_all),
        #lambda wildcards: barcode_files_for_set(wildcards.set),
    params:
        samples = samples_noNM3,
        cell_types = cell_types_noNM3,
        outdir = _results("atac_barcodes_sample_cluster/all"),
    shell:
        """
        Rscript workflow/scripts/barcodes_by_sample_by_cluster.R \
            --barcode_to_cell_type {input.barcode_to_cell_type} \
            --barcode_map {input.barcode_map} \
            --outdir {params.outdir}
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

rule bamtobed:
    input:
        bam = _results("atac_bam_cluster/{set}/atac_{cluster}.bam"),
    output:
        bed = _results("atac_bam_cluster/{set}/atac_{cluster}.bed"),
    container:
        "workflow/envs/porchard_snatac_general_20220107.sif"
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output.bed}
        """

rule run_macs2:
    input:
        bed = _results("atac_bam_cluster/{set}/atac_{cluster}.bed"),
    output:
        bdg = _results("peakcalls/{set}/{cluster}_treat_pileup.bdg"),
        peaks = _results("peakcalls/{set}/{cluster}_peaks.narrowPeak"),
        summits = _results("peakcalls/{set}/{cluster}_summits.bed"),
    params:
        outdir = _results("peakcalls/{set}"),
        cluster = "{cluster}",
    container:
        "workflow/envs/porchard_snatac_general_20220107.sif"
    shell:
        """
        macs2 callpeak -t {input.bed} \
            -f BED \
            -g hs \
            --outdir {params.outdir} \
            -n {params.cluster} \
            --bdg --SPMR \
            --nomodel \
            --shift -100 --extsize 200 \
            --seed 762873 \
            --call-summits \
            --keep-dup all
        """



rule remove_blacklist:
    input:
        summits = _results("peakcalls/{set}/{cluster}_summits.bed"),
    output:
        summits_extended1 = _results("peakcalls/{set}/{cluster}_summits_ext150.bed"),
        summits_extended2 = _results("peakcalls/{set}/{cluster}_summits_ext500.bed"),
        summits_extended_noblacklist1 = _results("peakcalls/{set}/{cluster}_summits_ext150_noblacklist.bed"),
        summits_extended_noblacklist2 = _results("peakcalls/{set}/{cluster}_summits_ext500_noblacklist.bed"),
    params:
        blacklist = config["blacklist"]
    shell:
        """
        #extend summits by 150bp
        awk 'BEGIN {{ OFS = "\t"}} {{if ($2-150 > 0) {{$2 -= 150}} else {{$2 = 0}}; $3 += 150; $4 = $1"_"$2"_"$3; print $0}}' {input.summits} > {output.summits_extended1}

        #extend summits by 500bp
        awk 'BEGIN {{ OFS = "\t"}} {{if ($2-500 > 0) {{$2 -= 500}} else {{$2 = 0}}; $3 += 500; $4 = $1"_"$2"_"$3; print $0}}' {input.summits} > {output.summits_extended2}

        #remove blacklisted regions
        bedtools intersect -a {output.summits_extended1} -b {params.blacklist} -v > {output.summits_extended_noblacklist1}
        bedtools intersect -a {output.summits_extended2} -b {params.blacklist} -v > {output.summits_extended_noblacklist2}
        """
