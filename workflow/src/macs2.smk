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
        #~~~~~ call peaks
        #expand(_results("peakcalls/{set}/{cluster}_summits_ext500_noblacklist.bed"), cluster=cell_types_noNM3),
        #expand(_results("peakcalls/{set}/{cluster}_summits.bed"), cluster=get_clusters_for_set("all"), set="all"),
        #expand(_results("peakcalls/{set}/{cluster}_summits.bed"), cluster=get_clusters_for_set("noNM3"), set="noNM3"),
        #~~~~~ bigwigs
        expand(_results("peakcalls/noNM3/{cluster}_treat_pileup.bw"), cluster=['acinar','alpha','beta','delta','ductal','endothelial','gamma','immune','stellate'])



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

rule bdg_to_bw:
    input:
        bdg = _results("peakcalls/{set}/{cluster}_treat_pileup.bdg"),
    output:
        bdg2 = _results("peakcalls/{set}/{cluster}_treat_pileup_clipped.bdg"),
        bw = _results("peakcalls/{set}/{cluster}_treat_pileup.bw"),
    params:
        chr_sizes = _resources("hg38/hg38_cvb4.chrom.sizes"),
    shell:
        """
        bedClip {input.bdg} {params.chr_sizes} {output.bdg2}
        bedGraphToBigWig {output.bdg2} {params.chr_sizes} {output.bw}
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
