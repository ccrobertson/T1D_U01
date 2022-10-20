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
        expand("results/master_barcode_map_{set}.rds", set="multiome_5GEX"),
        #expand("results/master_barcode_map_{set}.rds", set="3GEX_all"),

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


# """
# Rscript workflow/scripts/master_barcode_map.R \
#     --barcode_to_cell_type {input.barcode_to_cell_type} \
#     --run_design {input.run_design} \
#     --PSNG_threshold {params.PSNG_threshold} \
#     --outdir {params.outdir} \
#     {input.demuxlet_files}
# """
