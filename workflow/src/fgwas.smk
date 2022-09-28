#!/usr/bin/env python3

#fgwas user manual:
#https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf

from os.path import join
import os
from functools import partial
import re
import glob

name = config["name"]

_results = partial(os.path.join, "results", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")

configfile: "workflow/src/fgwas_annotations.yaml"
annotations = config["fgwas_annotations"].keys()

def get_annot_file(annot):
    return config["fgwas_annotations"][annot]



rule all:
    input:
        #~~~ symlinks
        expand(_results("annotations/hg38/{annot}.bed"), annot=annotations),
        #~~~ liftover
        expand(_results("annotations/hg19/{annot}.bed"), annot=annotations),
        #~~~ fgwas_prep
        expand(_results("fgwas_run/{gwas}-{build}/fgwas_input.txt.gz"), gwas="T1D", build="hg19"),
        #~~~ fgwas
        expand(_results("fgwas_run/{gwas}-{build}/{annot}/fgwas.params"), gwas="T1D", build="hg19", annot=annotations),



rule symlinks:
    input:
        bed = lambda wildcards: get_annot_file(wildcards.annot),
    output:
        bed = _results("annotations/hg38/{annot}.bed"),
    shell:
        """
        ln -s --force --relative {input.bed} {output.bed}
        """

rule liftover_annotations:
    input:
        bed = _results("annotations/hg38/{annot}.bed"),
    output:
        bed = _results("annotations/hg19/{annot}.bed"),
    params:
        chain = _resources("hg38ToHg19.over.chain.gz"),
    shell:
        """
        liftOver -bedPlus=6 {input.bed} resources/hg38ToHg19.over.chain.gz {output.bed} {output.bed}.unmapped
        """


rule fgwas_prep:
    input:
        gwas_file = _results("gwas_stats/{gwas}-{build}/formatted.txt.gz"),
        annots = expand(_results("annotations/{{build}}/{annot}.bed"), annot=annotations),
    output:
        fgwas_file = _results("fgwas_run/{gwas}-{build}/fgwas_input.txt.gz"),
    conda:
        "genetics"
    shell:
        """
        Rscript workflow/scripts/fgwas_prep.R \
            --gwas_file {input.gwas_file} \
            --outfile {output.fgwas_file} \
            {input.annots}
        """


rule fgwas_run:
    input:
        fgwas_input = _results("fgwas_run/{gwas}-{build}/fgwas_input.txt.gz"),
    output:
        _results("fgwas_run/{gwas}-{build}/{annot}/fgwas.llk"),
        _results("fgwas_run/{gwas}-{build}/{annot}/fgwas.params"),
    params:
        prefix = _results("fgwas_run/{gwas}-{build}/{annot}/fgwas"),
        annot = "{annot}",
    shell:
        """
        fgwas -i {input.fgwas_input} -w {params.annot} -o {params.prefix}
        """

# rule format_for_fgwas:
#     output:
#         "data/2020_chiou_t1d_gwas/all_snps.maf0.05.fgwas.txt.gz"
#     shell:
#         """
#         Rscript scripts/fgwas_prep_chiou_2020.R
#         """
#
# rule get_fine_input:
#     output:
#         "data/2020_chiou_t1d_gwas/for_fgwas/chiou_2020.lead_plus_5K_adj_snps.txt.gz"
#     shell:
#         """
#         Rscript scripts/fgwas_make_chunks.Rmd
#         """
#
# rule get_new_chunks:
#     output:
#         "data/2020_chiou_t1d_gwas/for_fgwas/chunks.bed"
#     shell:
#         """
#         Rscript scripts/fgwas_make_chunks.Rmd
#         """
#
# rule make_input:
#     input:
#         # i = lambda wildcards: get_file(wildcards.feature),
#         i = lambda wildcards: get_file_from_wildcards(
#             wildcards.cluster, wildcards.cond, wildcards.fc
#         ),
#         snps = rules.format_for_fgwas.output,
#     output:
#         # join(RESULTS, "input_files", "{feature}.txt.gz"),
#         join(RESULTS, "input_files", "{cluster}__{cond}__{fc}.txt.gz"),
#     params:
#         # handle = "{feature}"
#         handle = "{cluster}__{cond}__{fc}"
#     shell:
#         """
#         {IONICE} Rscript scripts/fgwas_get_input_chiou_2020_all_snps.R \
#             {input.i} {output} {params}
#         """
#
# rule fgwas:
#     input:
#         rules.make_input.output,
#     output:
#         join(MODEL_DIR, "fgwas.llk"),
#     params:
#         # col = "{feature}",
#         col = "{cluster}__{cond}__{fc}",
#         handle = join(MODEL_DIR, "fgwas"),
#     shell:
#         """
#         {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print
#         """
#
# rule fgwas_new_chunks:
#     input:
#         snps = rules.make_input.output,
#         bed = rules.get_new_chunks.output,
#     output:
#         join(BED_DIR, "fgwas.llk"),
#     params:
#         # col = "{feature}",
#         col = "{cluster}__{cond}__{fc}",
#         handle = join(BED_DIR, "fgwas"),
#     shell:
#         """
#         {IONICE} fgwas -i {input.snps} -bed {input.bed} -w {params.col} \
#             -o {params.handle} -print
#         """
#
# rule make_input_fine:
#     input:
#         i = lambda wildcards: get_file(wildcards.feature),
#         snps = rules.format_for_fgwas.output,
#     output:
#         join(RESULTS, "input_files_fine", "{feature}.txt.gz"),
#     params:
#         handle = "{feature}"
#     shell:
#         """
#         {IONICE} Rscript scripts/fgwas_get_input_chiou_2020_fine.R \
#             {input.i} {output} {params}
#         """
#
# rule fgwas_fine:
#     input:
#         rules.make_input_fine.output,
#     output:
#         join(FINE_DIR, "fgwas.llk"),
#     params:
#         col = "{feature}",
#         handle = join(FINE_DIR, "fgwas"),
#     shell:
#         """
#         {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print -fine
#         """
#
# rule make_big_input:
#     """Generate one file with all the output models"""
#     input:
#         expand(
#             # join(RESULTS, "input_files", "{feature}.txt.gz"),
#             # feature=features
#             join(RESULTS, "input_files", "{cluster}__{cond}__{fc}.txt.gz"),
#             cluster=clusters, cond=conds, fc=fcs
#         ),
#     output:
#         join(RESULTS, "input_files", "combined.txt.gz"),
#     shell:
#         """
#         {IONICE} Rscript scripts/fgwas_combine_annotations.R {input} {output}
#         """
#
# rule fgwas_combo:
#     input:
#         rules.make_big_input.output,
#     output:
#         join(MODEL2_DIR, "fgwas.llk"),
#     params:
#         # col = get_full_model(),
#         col = "{combo}",
#         handle = join(MODEL2_DIR, "fgwas"),
#     shell:
#         """
#         {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print
#         """
#
# def get_annotations_list(dir):
#     files = glob.glob(os.path.join(dir, "*.bed"))
#     annots = [re.sub(".bed","", re.sub("/","", re.sub(dir, "", x))) for x in files]
#     return(annots)
#
# IONICE = "ionice -c2 -n7"
#
# rule all:
#     input:
#         ### Liftover annotations
#         expand(_results("annotations/whole_islet-hg19/{annot}.bed"), annot = get_annotations_list(_results("annotations/whole_islet-hg38"))),
#         ### Run fgwas
#         _results("fgwas_run/T1D-whole_islet-hg19/fgwas_input.txt.gz"),
#         expand(_results("fgwas_run/T1D-whole_islet-hg19/{annot}/fgwas.llk"), annot = ['cres_overlapping_intergenic_atac_features','atac_summits_all', 'atac_summits_all+cres_overlapping_intergenic_atac_features']),
