#!/usr/bin/env python3

#fgwas user manual:
#https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf

from os.path import join
import os
from functools import partial
import re
import glob

_results = partial(os.path.join, "results")


configfile: "workflow/src/reformat_ricardo_bed_for_fgwas.yaml"

rule all:
    input:
        expand(_results("fgwas/annotations/ricardo_annots_hg19_reformatted/{file}.bed3"), file=config["ricardo_files"])

rule reformat_ricardo_annotations:
    input:
        bed = _results("fgwas/annotations/ricardo_annots_hg19/{file}"),
    output:
        bed = _results("fgwas/annotations/ricardo_annots_hg19_reformatted/{file}.bed3"),
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{print $1,$2,$3}}' {input.bed} > {output.bed}
        """
