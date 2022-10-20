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
sets = config["fgwas_sets"].keys()
annots = list(config["fgwas_annotation_sources"]["hg19"].keys()) + list(config["fgwas_annotation_sources"]["hg38"].keys())

wildcard_constraints:
   set="|".join(sets),
   annot="|".join(annots),



def get_annots_in_set(set):
    return config["fgwas_sets"][set]

def get_annot_source(build, annot):
    return config["fgwas_annotation_sources"][build][annot]

def get_gwas_file(gwas, build):
    return config["gwas_stats"][gwas][build]

rule all:
    input:
        #~~~ symlinks
        #expand(_results("annotations/hg19/{annot}.bed"), annot=config["fgwas_annotation_sources"]["hg19"].keys()),
        #expand(_results("annotations/hg38/{annot}.bed"), annot=config["fgwas_annotation_sources"]["hg38"].keys()),
        #~~~ liftover
        #expand(_results("annotations/hg19/{annot}.bed"), annot=config["fgwas_annotation_sources"]["hg38"].keys()),
        #expand(_results("annotations/hg38/{annot}.bed"), annot=config["fgwas_annotation_sources"]["hg19"].keys()),
        #~~~ remove mhc from annotations
        #expand(_results("annotations/{build}/{annot}_noMHC.bed"), build="hg19", annot=get_annots_in_set("tCRE")),
        #expand(_results("annotations/{build}/{annot}_noMHC.bed"), build="hg19", annot=get_annots_in_set("ATAC_1kb")),
        #expand(_results("annotations/{build}/{annot}_noMHC.bed"), build="hg19", annot=get_annots_in_set("ATAC_300b")),
        #expand(_results("annotations/{build}/{annot}_noMHC.bed"), build="hg19", annot=get_annots_in_set("Ricardo_clusters")),
        #expand(_results("annotations/{build}/{annot}_noMHC.bed"), build="hg38", annot=get_annots_in_set("Ricardo_clusters")),
        #~~~~~ format gwas
        #expand(_results("gwas_stats/{gwas}-{build}/formatted.txt.gz"), gwas=["T1D_Crouch", "T1D_Chiou"], build="hg38"),
        #~~~~~ remove mhc from gwas
        #expand(_results("gwas_stats/{gwas}-{build}/formatted_noMHC.txt.gz"), gwas="T1D_Crouch", build="hg38"),
        #~~~ fgwas_prep
        #expand(_results("fgwas_run/{gwas}-{build}-{set}/fgwas_input_{set}.txt.gz"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", set=["ATAC_1kb","ATAC_300b","Ricardo_clusters","tCRE"]),
        #expand(_results("fgwas_run/{gwas}-{build}-{set}/fgwas_input_{set}.txt.gz"), gwas=["T1D","T2D"], build="hg19", set=["ATAC_1kb","ATAC_300b","Ricardo_clusters","tCRE"]),
        #~~~ fgwas
        #expand(_results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas.params"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", annot=get_annots_in_set("tCRE"), set="tCRE"),
        #expand(_results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas.params"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", annot=get_annots_in_set("ATAC_300b"), set="ATAC_300b"),
        #expand(_results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas.params"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", annot=get_annots_in_set("ATAC_1kb"), set="ATAC_1kb"),
        #expand(_results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas.params"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", annot=get_annots_in_set("Ricardo_clusters"), set="Ricardo_clusters"),
        #~~~ forest plot
        expand(_results("fgwas_run/{gwas}-{build}-{set}.png"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", set=["ATAC_1kb","ATAC_300b","Ricardo_clusters","tCRE"]),
        expand(_results("fgwas_run/{gwas}-{build}-{set}_cond.png"), gwas=["T1D_Crouch","T1D_Chiou"], build="hg38", set=["ATAC_1kb","ATAC_300b","Ricardo_clusters","tCRE"])

# NOTE: SNAKEMAKE CAN'T FIGURE OUT WHETHER TO RUN SYMLINKS OR LIFTOVER FIRST
# MAY HAVE TO SET UP SYMLINKS AS A PRELIM STEP OUTSIDE OF PIPELINE
# rule symlinks:
#     input:
#         bed = lambda wildcards: get_annot_source(wildcards.build, wildcards.annot),
#     output:
#         bed = _results("annotations/{build}/{annot}.bed"),
#     shell:
#         """
#         ln -s --force --relative {input.bed} {output.bed}
#         """


# rule liftover_to_hg19:
#     input:
#         bed = _results("annotations/hg38/{annot}.bed"),
#     output:
#         bed = _results("annotations/hg19/{annot}.bed"),
#     params:
#         chain = _resources("hg38ToHg19.over.chain.gz"),
#     shell:
#         """
#         liftOver -bedPlus=6 {input.bed} {params.chain} {output.bed} {output.bed}.unmapped
#         """
#
# rule liftover_to_hg38:
#     input:
#         bed = _results("annotations/hg19/{annot}.bed"),
#     output:
#         bed = _results("annotations/hg38/{annot}.bed"),
#     params:
#         chain = _resources("hg19ToHg38.over.chain.gz"),
#     shell:
#         """
#         liftOver {input.bed} {params.chain} {output.bed} {output.bed}.unmapped
#         """


rule remove_mhc_annot:
    input:
        bed = _results("annotations/{build}/{annot}.bed"),
    output:
        bed = _results("annotations/{build}/{annot}_noMHC.bed"),
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} $1!="chr6" {{print}} $1=="chr6" && $2<20000000 {{print}} $1=="chr6" && $2>40000000 {{print}}' {input.bed} > {output.bed}
        """

rule format_gwas:
    input:
        gwas_file = lambda wildcards: get_gwas_file(wildcards.gwas, wildcards.build),
    output:
        gwas_file = _results("gwas_stats/{gwas}-{build}/formatted.txt.gz"),
    shell:
        """
        Rscript workflow/scripts/format_gwas_stats_for_fgwas.R \
            --gwas_file {input.gwas_file} \
            --outfile {output.gwas_file}
        """

rule remove_mhc_gwas:
    input:
        gwas_file = _results("gwas_stats/{gwas}-{build}/formatted.txt.gz"),
    output:
        gwas_file = _results("gwas_stats/{gwas}-{build}/formatted_noMHC.txt.gz"),
    conda:
        "genetics"
    shell:
        """
        zcat {input.gwas_file} | awk 'BEGIN {{OFS="\t"}} $2!="chr6" {{print}} $2=="chr6" && $3<20000000 {{print}} $2=="chr6" && $3>40000000 {{print}}' | bgzip > {output.gwas_file}
        """

rule fgwas_prep:
    input:
        gwas_file = _results("gwas_stats/{gwas}-{build}/formatted_noMHC.txt.gz"),
        annots = lambda wildcards: expand(_results("annotations/{{build}}/{annot}_noMHC.bed"), annot=get_annots_in_set(wildcards.set))
    output:
        fgwas_file = _results("fgwas_run/{gwas}-{build}-{set}/fgwas_input_{set}.txt.gz"),
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
        fgwas_input = _results("fgwas_run/{gwas}-{build}-{set}/fgwas_input_{set}.txt.gz"),
    output:
        _results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas.llk"),
        _results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas.params"),
    params:
        prefix = _results("fgwas_run/{gwas}-{build}-{set}/{annot}/fgwas"),
        annot = "{annot}_noMHC",
    shell:
        """
        fgwas -i {input.fgwas_input} -w {params.annot} -o {params.prefix}
        """

rule aggregate:
    input:
        results = lambda wildcards: expand(_results("fgwas_run/{{gwas}}-{{build}}-{{set}}/{annot}/fgwas.params"), annot=get_annots_in_set(wildcards.set)),
    output:
        txt = _results("fgwas_run/{gwas}-{build}-{set}.txt"),
        plotfile = _results("fgwas_run/{gwas}-{build}-{set}.png"),
    shell:
        """
        cat {input.results} | grep -v parameter | grep -v pi_region > {output.txt}
        Rscript workflow/scripts/forest_plot.R \
            --input {output.txt} \
            --plotfile {output.plotfile}
        """


rule fgwas_run_cond:
    input:
        fgwas_input = _results("fgwas_run/{gwas}-{build}-{set}/fgwas_input_{set}.txt.gz"),
    output:
        _results("fgwas_run/{gwas}-{build}-{set}/{annot}_cond/fgwas.llk"),
        _results("fgwas_run/{gwas}-{build}-{set}/{annot}_cond/fgwas.params"),
    params:
        prefix = _results("fgwas_run/{gwas}-{build}-{set}/{annot}_cond/fgwas"),
        annot = "{annot}_noMHC",
        cond_annot = "wholeIslet_{set}_noMHC",
    shell:
        """
        fgwas -i {input.fgwas_input} -w {params.annot} -cond {params.cond_annot} -o {params.prefix}
        """


rule aggregate_cond:
    input:
        results = lambda wildcards: expand(_results("fgwas_run/{{gwas}}-{{build}}-{{set}}/{annot}_cond/fgwas.params"), annot=get_annots_in_set(wildcards.set)),
    output:
        txt = _results("fgwas_run/{gwas}-{build}-{set}_cond.txt"),
        plotfile = _results("fgwas_run/{gwas}-{build}-{set}_cond.png"),
    shell:
        """
        cat {input.results} | grep -v parameter | grep -v pi_region > {output.txt}
        Rscript workflow/scripts/forest_plot.R \
            --input {output.txt} \
            --plotfile {output.plotfile}
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
