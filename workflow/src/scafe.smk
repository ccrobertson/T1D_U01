#!/usr/bin/env python3

from os.path import join
import os
import pandas as pd
from functools import partial
import math
import glob


name = config["name"]
_results = partial(os.path.join, "results", name)
_data = partial(os.path.join, "data")
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


def get_samples_for_set(set):
    samples = list(config["sets"][set]["samples"].keys())
    return samples

def get_clusters_for_set(set):
    d = pd.read_csv("results/liger/"+set+"/cluster_to_cell_type.csv")
    return d['cell_type'].unique().tolist()


configfile: _data("batch_design.json")

IONICE = 'ionice -c2 -n7'

rule all:
    input:
        expand(_results("cellranger_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt"), sample=get_samples_for_set("multiome_5GEX"), cluster=get_clusters_for_set("multiome_5GEX"), set="multiome_5GEX"),
        #expand(_results("scafe_by_cluster/{set}/{cluster}/ctss_to_bigwig/{cluster}/wig/{cluster}.cpm.fwd.bw"), sample=get_samples_for_set("multiome_5GEX"), cluster=get_clusters_for_set("multiome_5GEX"), set="multiome_5GEX"),
        expand(_results("scafe_by_sample_by_cluster/{set}/{sample}_{cluster}/ctss_to_bigwig/{sample}_{cluster}/bed/{sample}_{cluster}.unencoded_G.collapse.ctss.bed.gz"), sample=get_samples_for_set("multiome_5GEX"), cluster=get_clusters_for_set("multiome_5GEX"), set="multiome_5GEX"),


rule bgzip:
    input:
        barcodes ="results/fiveprime/barcode_filtering/{sample}/barcodes_nuclei.txt",
    output:
        barcodes_gz = "results/fiveprime/barcode_filtering/{sample}/barcodes_nuclei.txt.gz",
    conda:
        "genetics"
    shell:
        """
        bgzip -c {input.barcodes} > {output.barcodes_gz}
        """

rule scafe_by_sample:
    input:
        #bam = "results/fiveprime/starsolo/{sample}/starsolo.Solo.out/starsolo.Aligned.sortedByCoord.out.bam",
        bam = _results("cellranger/{sample}/outs/possorted_genome_bam.bam"),
        barcodes ="results/fiveprime/barcode_filtering/{sample}/barcodes_nuclei.txt.gz",
    output:
        cpm_fw = _results("scafe_by_sample/{sample}/ctss_to_bigwig/{sample}/wig/{sample}.cpm.fwd.bw"),
        cpm_rev = _results("scafe_by_sample/{sample}/ctss_to_bigwig/{sample}/wig/{sample}.cpm.rev.bw"),
    params:
        tag = "{sample}",
        outdir = _results("scafe_by_sample/{sample}")
    shell:
        """
        scafe.workflow.sc.solo \
          --run_bam_path={input.bam} \
          --run_cellbarcode_path={input.barcodes} \
          --genome=hg38.gencode_v32 \
          --run_tag={params.tag} \
          --run_outDir={params.outdir} \
          --max_thread=10 \
          --overwrite=yes
        """


rule cellranger_barcodes_by_sample_by_cluster:
    input:
        barcode_to_cell_type = "results/liger/{set}/barcode_to_cluster_to_cell_type.csv",
    output:
        barcodes = _results("cellranger_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt"),
    params:
        sample = "{sample}",
        cluster = "{cluster}",
    shell:
        """
        sample=$(echo "{params.sample}" | sed 's/-/_/g')
        awk -v s="$sample" -v c="{params.cluster}" 'BEGIN {{FS=","}} $1==s && $4==c {{ print $2"-1"}}' {input.barcode_to_cell_type} > {output.barcodes}
        """

rule cellranger_bgzip:
    input:
        barcodes = _results("cellranger_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt"),
    output:
        barcodes_gz = _results("cellranger_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt.gz"),
    conda:
        "genetics"
    shell:
        """
        bgzip -c {input.barcodes} > {output.barcodes_gz}
        """

rule cellranger_bam_by_sample_by_cluster:
    input:
        bam = "results/fiveprime/cellranger/{sample}/outs/possorted_genome_bam.bam",
        barcodes = _results("cellranger_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt"),
    output:
        bam = _results("cellranger_bam_sample_cluster/{set}/{sample}_{cluster}.bam"),
        bai = _results("cellranger_bam_sample_cluster/{set}/{sample}_{cluster}.bam.bai"),
    shell:
        """
        subset-bam --bam {input.bam} --cell-barcodes {input.barcodes} --out-bam {output.bam} --cores 10
        samtools index {output.bam}
        """

# rule merge_bams:
#     input:
#         bams = lambda wildcards: expand(_results("cellranger_bam_sample_cluster/{{set}}/{sample}_{{cluster}}.bam"), sample=get_samples_for_set(wildcards.set)),
#     output:
#         bam = _results("cellranger_bam_cluster/{set}/{cluster}.bam"),
#         bai = _results("cellranger_bam_cluster/{set}/{cluster}.bam.bai"),
#     container:
#         "workflow/envs/arushi_general.simg"
#     shell:
#         """
#     	samtools merge -@ 15 -O BAM {output.bam} {input.bams}
#         samtools index {output.bam}
#         """



# WARNING!!! when rerunning SCAFE must first delete old files
# as it will not restart a step if it detects output
### ADD A CHECK STEP HERE: if exists(outdir); then Error
rule scafe_by_sample_by_cluster:
    input:
        bam = _results("cellranger_bam_sample_cluster/{set}/{sample}_{cluster}.bam"),
        barcodes = _results("cellranger_barcodes_sample_cluster/{set}/barcodes_{sample}_{cluster}.txt"),
    output:
        bed = _results("scafe_by_sample_by_cluster/{set}/{sample}_{cluster}/ctss_to_bigwig/{sample}_{cluster}/bed/{sample}_{cluster}.unencoded_G.collapse.ctss.bed.gz"),
    params:
        tag = "{sample}_{cluster}",
        outdir = _results("scafe_by_sample_by_cluster/{set}/{sample}_{cluster}")
    shell:
        """
        scafe.workflow.sc.solo \
          --run_bam_path={input.bam} \
          --run_cellbarcode_path={input.barcodes} \
          --genome=hg38.gencode_v32 \
          --run_tag={params.tag} \
          --run_outDir={params.outdir} \
          --max_thread=10 \
          --overwrite=yes
        """


# NOTE: model lib_list_path.txt after
# /lab/work/ccrober/sw/SCAFE/demo/input/cm.aggregate/lib_list_path.txt
#
# Columns in scafe_by_cluster/{set}/{cluster}_lib_list_path.txt will be:
#   col1: {sample}_{cluster}
#   col2: scafe_by_sample/{sample}_{cluster}/remove_strand_invader/{sample}_{cluster}/bed/{sample}_{cluster}.pass.ctss.bed.gz
#   col3: scafe_by_sample/{sample}_{cluster}/bam_to_ctss/{sample}_{cluster}/bed/{sample}_{cluster}.unencoded_G.collapse.ctss.bed.gz
# rule scafe_aggregate:
#     input:
#         beds = lambda wildcards:
#             expand(_results("scafe_by_sample_by_cluster/{{set}}/{sample}_{{cluster}}/ctss_to_bigwig/{sample}_{{cluster}}/bed/{sample}_{{cluster}}.unencoded_G.collapse.ctss.bed.gz"),
#             sample=get_samples_for_set(wildcards.set)),
#         lib_list_path = _results("scafe_by_cluster/{set}/{cluster}_lib_list_path.txt"),
#     output:
#         cpm_fw = _results("scafe_by_cluster/{set}/{cluster}/ctss_to_bigwig/{cluster}/wig/{cluster}.cpm.fwd.bw"),
#         cpm_rev = _results("scafe_by_cluster/{set}/{cluster}/ctss_to_bigwig/{cluster}/wig/{cluster}.cpm.rev.bw"),
#     params:
#         tag = "{cluster}",
#         outdir = _results("scafe_by_cluster/{set}/{cluster}"),
#     shell:
#         """
#         scafe.workflow.cm.aggregate \
#           --lib_list_path {input.lib_list_path} \
#           --genome=hg38.gencode_v32 \
#           --run_tag={params.tag} \
#           --run_outDir={params.outdir} \
#           --max_thread=10 \
#           --overwrite=yes
#         """
#
# rule unzip:
#     input:
#         bed = _results("scafe_by_cluster/{set}/{cluster}/annotate/{cluster}/bed/{cluster}.CRE.annot.bed.gz"),
#     output:
#         bed = _results("scafe_by_cluster/{set}/{cluster}/annotate/{cluster}/bed/{cluster}.CRE.annot.bed"),
#     conda:
#         "genetics"
#     shell:
#         """
#         bgzip -d -c {input.bed} > {output.bed}
#         """

#
# rule scafe_by_cluster:
#     input:
#         bam = _results("bam_cluster/{cluster}.bam"),
#         barcodes = _results("barcodes_cluster/barcodes_Sample_5125-NM-1_NM-1_{cluster}.txt"),
#     output:
#         cpm_fw = _results("scafe_by_cluster/{cluster}/ctss_to_bigwig/{cluster}/wig/{cluster}.cpm.fwd.bw"),
#         cpm_rev = _results("scafe_by_cluster/{cluster}/ctss_to_bigwig/{cluster}/wig/{cluster}.cpm.rev.bw"),
#     params:
#         tag = "{cluster}",
#         outdir = _results("scafe_by_cluster/{cluster}"),
#     shell:
#         """
#         scafe.workflow.sc.solo \
#           --run_bam_path={input.bam} \
#           --run_cellbarcode_path={input.barcodes} \
#           --genome=hg38.gencode_v32 \
#           --run_tag={params.tag} \
#           --run_outDir={params.outdir}
#         """

# OLD; using starsolo
# rule starsolo_barcodes_by_sample_by_cluster:
#     input:
#         barcode_to_cell_type = "results/liger/multimodal_noNM3/barcode_to_cluster_to_cell_type.csv",
#     output:
#         barcodes = _results("starsolo_barcodes_sample_cluster/barcodes_{sample}_{cluster}.txt"),
#     params:
#         sample = "{sample}",
#         cluster = "{cluster}",
#     shell:
#         """
#         sample=$(echo "{params.sample}" | sed 's/-/_/g')
#         awk -v s="$sample" -v c="{params.cluster}" 'BEGIN {{FS=","}} $1==s && $4==c {{ print $2}}' {input.barcode_to_cell_type} > {output.barcodes}
#         """
#
# rule bam_by_sample_by_cluster:
#     input:
#         bam = "results/fiveprime/starsolo/{sample}/starsolo.Aligned.sortedByCoord.out.bam",
#         barcodes = _results("starsolo_barcodes_sample_cluster/barcodes_{sample}_{cluster}.txt"),
#     output:
#         bam = _results("starsolo_bam_sample_cluster/{sample}_{cluster}.bam"),
#     container:
#         "workflow/envs/arushi_general.simg"
#     shell:
#         """
#         samtools view -h -b -@ 10 -D CB:{input.barcodes} {input.bam} > {output.bam}
#         samtools index {output.bam}
#         """


# OLD; Using cellranger output
# rule scafe:
#     input:
#         bam = _results("cellranger/{sample}/outs/possorted_genome_bam.bam"),
#         barcodes = _results("cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),
#     output:
#         cpm_fw = _results("scafe/{sample}/ctss_to_bigwig/{sample}/wig/{sample}.cpm.fwd.bw"),
#         cpm_rev = _results("scafe/{sample}/ctss_to_bigwig/{sample}/wig/{sample}.cpm.rev.bw"),
#     params:
#         sample = "{sample}",
#         outdir = _results("scafe/{sample}")
#     shell:
#         """
#         scafe.workflow.sc.solo \
#           --run_bam_path={input.bam} \
#           --run_cellbarcode_path={input.barcodes} \
#           --genome=hg38.gencode_v32 \
#           --run_tag={params.sample} \
#           --run_outDir={params.outdir}
#         """
