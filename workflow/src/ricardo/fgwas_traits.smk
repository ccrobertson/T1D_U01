from os.path import join
from itertools import combinations


def get_file(wildcards):
    """Gets the associated bed file for each feature"""
    my_feature = wildcards.feature
    f = config["features"][my_feature][0]
    f = join(BASE_DIR, f)
    return f


def get_file_from_wildcards(wildcards):
    """Gets the associated bed file for each wildcard combo"""
    f = DA_STRING.format(
        cluster=wildcards.cluster, 
        cond=wildcards.cond, 
        fc=wildcards.fc, 
        direction=wildcards.direction
    )
    return f


IONICE = "ionice -c2 -n7"

# Dirs
LIGER = config["liger_handle"]
BASE_DIR = join(config["results"], LIGER)
SUMMSTATS_DIR = config["summary_stats_dir"]
DA_STRING = join(BASE_DIR, config["summits_da_string"])

RESULTS = join(BASE_DIR, "fgwas", "mast_random_effect_0nAGQ_model_no_min_new_traits")
# MODEL_DIR = join(RESULTS, "output_models", "{trait}", "{feature}")
MODEL_DIR = join(
    RESULTS, "output_models", 
    "{trait}", "{cluster}__{cond}__{direction}__{fc}"
)

# Wildcards
clusters = config["clusters"]
conds = config["conditions"]
traits = [x for x in config["GWAS"].keys()]
# fcs = ["0.01", "0.05", "0.1"]
fcs = ["0.8", "0.85", "0.9", "0.95", "0.99"]
directions = ["above", "below"]
# features = [x for x in config["features"].keys()]


rule all:
    input:
        expand(
            join(MODEL_DIR, "fgwas.llk"), 
            # feature=features
            cluster=clusters, cond=conds, fc=fcs, trait=traits, 
            direction=directions
        ),
        

rule format_for_fgwas:
    output:
        expand(join(SUMMSTATS_DIR, "{trait}.txt.gz"), trait=traits)
    shell:
        """
        Rscript scripts/fgwas/prep_FGluadjBMI_MAGIC.R  
        Rscript scripts/fgwas/prep_FInsadjBMI_MAGIC.R  
        Rscript scripts/fgwas/prep_T2D_DIAMANTE.R
        """

rule make_input:
    input:
        snps = join(SUMMSTATS_DIR, "{trait}.txt.gz"),
        # i = get_file,
        i = get_file_from_wildcards,
    output:
        # join(RESULTS, "input_files", "{trait}", "{feature}.txt.gz"),
        join(
            RESULTS, "input_files", "{trait}", 
            "{cluster}__{cond}__{direction}__{fc}.txt.gz"
        ),
    params:
        # handle = "{feature}"
        handle = "{cluster}__{cond}__{direction}__{fc}"
    shell:
        """
        {IONICE} Rscript scripts/fgwas/prepare_fgwas_input.R \
            {input.snps} {input.i} {output} {params}
        """

rule fgwas:
    input:
        rules.make_input.output,
    output:
        join(MODEL_DIR, "fgwas.llk"),
    params:
        # col = "{feature}",
        col = "{cluster}__{cond}__{direction}__{fc}",
        handle = join(MODEL_DIR, "fgwas"),
    shell:
        """
        {IONICE} fgwas -i {input} -w {params.col} -o {params.handle} -print
        """