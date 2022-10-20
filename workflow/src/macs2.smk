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


configfile: _data("batch_design.json")
configfile: "workflow/src/markers.yaml"

def get_bam(set, sample):
    bam = config["sets"][set]["samples"][sample]
    return bam

def get_samples_for_set(set):
    samples = list(config["sets"][set]["samples"].keys())
    return samples

def get_clusters_for_set(set):
    d = pd.read_csv("results/liger/"+set+"/cluster_to_cell_type.csv")
    return d['cell_type'].unique().tolist()


rule all:
    input:
        #~~~~~ get barcode lists
        expand(_results("atac_barcodes_sample_cluster/{set}/atac_barcodes_{sample}_{cluster}.txt"), sample=get_samples_for_set("multiome_5GEX"), cluster=get_clusters_for_set("multiome_5GEX")),
        #expand(_results("atac_barcodes_sample_cluster/all/atac_barcodes_{sample}_{cluster}.txt"), sample=samples_all, cluster=cell_types_all),
        #~~~~~ filter bams
        #expand(_results("atac_bam_sample_cluster/noNM3/atac_{sample}_{cluster}.bam"),sample=samples_noNM3, cluster=cell_types_noNM3),
        #expand(_results("atac_bam_sample_cluster/all/atac_{sample}_{cluster}.bam"),sample=samples_all, cluster=cell_types_all),
        #~~~~~ merge bams
        #expand(_results("atac_bam_cluster/noNM3/atac_{cluster}.bam"), cluster=cell_types_noNM3),
        #expand(_results("atac_bam_cluster/all/atac_{cluster}.bam"), cluster=cell_types_all),
        #~~~~~ call peaks
        expand(_results("peakcalls/{set}/{cluster}_summits_ext500_noblacklist.bed"), cluster=cell_types_noNM3),
        #expand(_results("peakcalls/{set}/{cluster}_summits.bed"), cluster=get_clusters_for_set("all"), set="all"),
        #expand(_results("peakcalls/{set}/{cluster}_summits.bed"), cluster=get_clusters_for_set("noNM3"), set="noNM3"),



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
