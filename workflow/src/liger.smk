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

clusters = ["ductal", "acinar", "alpha", "beta", "delta", "gamma", "stellate", "endothelial", "immune"]


def get_samples_for_set(set):
    samples = config["liger_sets"][set]
    return samples

def get_count_files_for_set(set):
    samples = config["liger_sets"][set]
    counts_files = [config["counts"][s] for s in samples]
    return counts_files


rule all:
    input:
        expand(_results("{set}/umap_INS.png"), set=["multimodal_all", "multimodal_noNM3", "3GEX_all", "3GEX_noNM3"]),
        expand(_results("{set}/liger_clusters.tsv"), set=["multimodal_all", "multimodal_noNM3", "3GEX_all", "3GEX_noNM3"]),
        expand(_results("{set}/cluster_to_cell_type.csv"), set=["multimodal_all", "multimodal_noNM3", "3GEX_all", "3GEX_noNM3"]),
        expand(_results("{set}/barcode_to_cluster_to_cell_type.csv"), set=["multimodal_all", "multimodal_noNM3", "3GEX_all", "3GEX_noNM3"]),


rule liger_iNMF:
    input:
        count_matrices = lambda wildcards: get_count_files_for_set(wildcards.set),
    output:
        _results("{set}/liger_obj.rds"),
    params:
        outdir = _results("{set}"),
    shell:
        """
        Rscript workflow/scripts/run_liger_iNMF.R --outdir {params.outdir} {input.count_matrices}
        """

rule liger_clustering:
    input:
        liger_obj = _results("{set}/liger_obj.rds"),
    output:
        _results("{set}/liger_obj_clusters.rds"),
    params:
        outdir = _results("{set}"),
    shell:
        """
        Rscript workflow/scripts/run_liger_clustering.R --liger_obj {input.liger_obj} --outdir {params.outdir}
        """

rule liger_plots:
    input:
        liger_obj = _results("{set}/liger_obj_clusters.rds"),
    output:
        _results("{set}/umap_INS.png"),
    params:
        outdir = _results("{set}"),
    shell:
        """
        Rscript workflow/scripts/run_liger_plots.R --liger_obj {input.liger_obj} --outdir {params.outdir}
        """

rule export_clusters:
    input:
        liger_obj = _results("{set}/liger_obj_clusters.rds"),
    output:
        tsv = _results("{set}/liger_clusters.tsv"),
    shell:
        """
        Rscript workflow/scripts/get_liger_clusters.R {input.liger_obj} {output.tsv}
        """

rule cluster_to_cell_type:
    input:
        barcode_to_cluster = _results("{set}/liger_clusters.tsv"),
        deg = _results("{set}/liger_deg_by_cluster.rds"),
    output:
        _results("{set}/cluster_to_cell_type.csv"),
        _results("{set}/barcode_to_cluster_to_cell_type.csv"),
    params:
        outdir = _results("{set}"),
        markers = "workflow/src/markers.yaml",
    shell:
        """
        Rscript workflow/scripts/cluster_to_cell_type.R \
            --barcode_to_cluster {input.barcode_to_cluster} \
            --deg {input.deg} \
            --markers {params.markers} \
            --outdir {params.outdir}
        """
