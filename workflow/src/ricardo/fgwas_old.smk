from os.path import join
from itertools import combinations

#fgwas user manual:
#https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf

def get_file(feature):
    """Gets the associated bed file for each feature"""
    f = config["features"][feature][0]
    f = join(BASE_DIR, f)
    return f


def get_file_from_wildcards(cluster, cond, fc):
    """Gets the associated bed file for each wildcard combo"""
    f = DA_STRING.format(cluster=cluster, cond=cond, fc=fc)
    return f


def get_combinations(features):
    """
    Makes all model combinations - right now only making the full model
    """
    l = len(features)
    models = []
    # for r in range(2, l+1):
    for r in range(l, l+1):
        iters = combinations(features, r)
        model_partial = ["+".join(i) for i in iters]
        models.append(model_partial)
    models = [item for sublist in models for item in sublist]
    return models


def get_full_model():
    """
    Making the full model
    """
    return "+".join(features_fullmodel)


IONICE = "ionice -c2 -n7"

# Dirs
LIGER = config["liger_handle"]
BASE_DIR = join(config["results"], LIGER)

# RESULTS = join(BASE_DIR, "fgwas_summits_all_snps")
# MODEL_DIR = join(RESULTS, "output_models", "{feature}")
# FINE_DIR = join(RESULTS, "output_models_fine", "{feature}")
# BED_DIR = join(RESULTS, "output_new_chunks", "{feature}")

RESULTS = join(BASE_DIR, "fgwas", "fgwas_summits_all_fcs_deseq2_min3")
# RESULTS = join(BASE_DIR, "fgwas_henry_random_model_FDR_notDAR")
MODEL_DIR = join(RESULTS, "output_models", "{cluster}__{cond}__{fc}")
BED_DIR = join(RESULTS, "output_new_chunks", "{cluster}__{cond}__{fc}")
FINE_DIR = join(RESULTS, "output_models_fine", "{cluster}__{cond}__{fc}")


MODEL2_DIR = join(RESULTS, "combined_models", "{combo}")

DA_STRING = join(BASE_DIR, config["summits_da_string"])

# Wildcards
conds = ["AAB"]
# clusters = ["GCG", "INS", "KRT19", "PDGFRA", "PPY", "PRSS1", "RGS5",
#             "SDS", "SST", "VWF"]
clusters = ["INS"]
fcs = ["0.8", "0.85", "0.9", "0.95", "0.99"]
# fcs = ["0.01", "0.05", "0.1"]


features = [x for x in config["features"].keys()]
features_fullmodel = [x for x in config["features"].keys()]

sig_features = ["immune", "stellateQ", "acinar", "ductal", "beta", "alpha"]
combos = "+".join([x + "_summits" for x in sig_features])

wildcard_constraints:
    cluster = "\w+",
    fc = "0\.[0-9]+",
    cond = "\w+",


rule all:
    input:
        expand(
            join(BED_DIR, "fgwas.llk"),
            cluster=clusters, cond=conds, fc=fcs
        ),
        # expand(join(MODEL_DIR, "fgwas.llk"), feature=features),
        # expand(join(FINE_DIR, "fgwas.llk"), feature=features),
        # expand(join(BED_DIR, "fgwas.llk"), feature=features),
        # expand(join(MODEL2_DIR, "fgwas.llk"), combo=combos),
        # expand(join(MODEL_DIR, "highest_ppa.bed"), feature=features),
        # # expand(
        # #     join(MODEL_DIR, "highest_ppa.bed"),
        # #     cluster=clusters, cond=conds, fc=fcs
        # # ),
        # expand(join(MODEL2_DIR, "highest_ppa.bed"), combo=combos),


rule format_for_fgwas:
    output:
        "data/2020_chiou_t1d_gwas/all_snps.maf0.05.fgwas.txt.gz"
    shell:
        """
        Rscript scripts/fgwas_prep_chiou_2020.R
        """

rule get_fine_input:
    output:
        "data/2020_chiou_t1d_gwas/for_fgwas/chiou_2020.lead_plus_5K_adj_snps.txt.gz"
    shell:
        """
        Rscript scripts/fgwas_make_chunks.Rmd
        """

rule get_new_chunks:
    output:
        "data/2020_chiou_t1d_gwas/for_fgwas/chunks.bed"
    shell:
        """
        Rscript scripts/fgwas_make_chunks.Rmd
        """

rule make_input:
    input:
        # i = lambda wildcards: get_file(wildcards.feature),
        i = lambda wildcards: get_file_from_wildcards(
            wildcards.cluster, wildcards.cond, wildcards.fc
        ),
        snps = rules.format_for_fgwas.output,
    output:
        # join(RESULTS, "input_files", "{feature}.txt.gz"),
        join(RESULTS, "input_files", "{cluster}__{cond}__{fc}.txt.gz"),
    params:
        # handle = "{feature}"
        handle = "{cluster}__{cond}__{fc}"
    shell:
        """
        {IONICE} Rscript scripts/fgwas_get_input_chiou_2020_all_snps.R \
            {input.i} {output} {params}
        """

rule fgwas:
    input:
        rules.make_input.output,
    output:
        join(MODEL_DIR, "fgwas.llk"),
    params:
        # col = "{feature}",
        col = "{cluster}__{cond}__{fc}",
        handle = join(MODEL_DIR, "fgwas"),
    shell:
        """
        {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print
        """

rule fgwas_new_chunks:
    input:
        snps = rules.make_input.output,
        bed = rules.get_new_chunks.output,
    output:
        join(BED_DIR, "fgwas.llk"),
    params:
        # col = "{feature}",
        col = "{cluster}__{cond}__{fc}",
        handle = join(BED_DIR, "fgwas"),
    shell:
        """
        {IONICE} fgwas -i {input.snps} -bed {input.bed} -w {params.col} \
            -o {params.handle} -print
        """

rule make_input_fine:
    input:
        i = lambda wildcards: get_file(wildcards.feature),
        snps = rules.format_for_fgwas.output,
    output:
        join(RESULTS, "input_files_fine", "{feature}.txt.gz"),
    params:
        handle = "{feature}"
    shell:
        """
        {IONICE} Rscript scripts/fgwas_get_input_chiou_2020_fine.R \
            {input.i} {output} {params}
        """

rule fgwas_fine:
    input:
        rules.make_input_fine.output,
    output:
        join(FINE_DIR, "fgwas.llk"),
    params:
        col = "{feature}",
        handle = join(FINE_DIR, "fgwas"),
    shell:
        """
        {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print -fine
        """

rule make_big_input:
    """Generate one file with all the output models"""
    input:
        expand(
            # join(RESULTS, "input_files", "{feature}.txt.gz"),
            # feature=features
            join(RESULTS, "input_files", "{cluster}__{cond}__{fc}.txt.gz"),
            cluster=clusters, cond=conds, fc=fcs
        ),
    output:
        join(RESULTS, "input_files", "combined.txt.gz"),
    shell:
        """
        {IONICE} Rscript scripts/fgwas_combine_annotations.R {input} {output}
        """

rule fgwas_combo:
    input:
        rules.make_big_input.output,
    output:
        join(MODEL2_DIR, "fgwas.llk"),
    params:
        # col = get_full_model(),
        col = "{combo}",
        handle = join(MODEL2_DIR, "fgwas"),
    shell:
        """
        {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print
        """

# rule make_bed_files:
#     input:
#         rules.fgwas.output,
#     output:
#         join(MODEL_DIR, "highest_ppa.bed"),
#     params:
#         name = "{feature}"
#         # name = "{cluster}__{cond}__{fc}"
#     shell:
#         """
#         {IONICE} Rscript scripts/fgwas_get_highest_ppa.R {input} {output} {params}
#         """

# rule make_bed_files_combo:
#     input:
#         rules.fgwas_combo.output,
#     output:
#         join(MODEL2_DIR, "highest_ppa.bed"),
#     params:
#         name = "{combo}"
#     shell:
#         """
#         {IONICE} Rscript scripts/fgwas_get_highest_ppa.R {input} \
#             {output} {params}
#         """
