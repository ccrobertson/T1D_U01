from os.path import join

IONICE = config["ionice"]

# Directories
MAINDIR = config["maindir"]
GARFIELD_SW = config["garfield"]
GARFIELD_DATA = config["garfield_data"]
FGWAS_DIR = config["fgwas_sumary_stats_dir"]

RESULTS = join(MAINDIR, "garfield", "mast_random_effect_0nAGQ_model_no_min")
SUMMSTATS_DIR = join(GARFIELD_DATA, "pval", "{trait}")
ANNOT_DIR = join(RESULTS, "annotation", "{trait}")
PREP_CHR_DIR = join(RESULTS, "prep-chr", "{trait}")
OUTDIR = join(RESULTS, "output", "{trait}")

# Wildcards
chroms = ["chr{}".format(x) for x in range(1, 23)]
traits = config["GWAS"]
clusters = config["clusters"]
conds = config["conditions"]
fcs = ["0.8", "0.85", "0.9", "0.95", "0.99"]
# fcs = ["0.8", "0.99"]
directions = ["above", "below"]

# List of input features
summits = expand(
    join(MAINDIR, config["summit_string"]), 
    cluster=clusters
)
dars = expand(
    join(MAINDIR, config["da_string"]), 
    cluster=clusters, cond=conds, direction=directions, fc=fcs
)


rule all:
    input:
        # expand(join(SUMMSTATS_DIR, "{chrom}"), trait=traits, chrom=chroms),
        expand(join(ANNOT_DIR, "link_file.txt"), trait=traits),
        expand(join(PREP_CHR_DIR, "{chrom}"), trait=traits, chrom=chroms),
        expand(join(PREP_CHR_DIR, "prepfile.txt"), trait=traits),
        expand(join(OUTDIR, "garfield.c0.out.txt"), trait=traits),


rule prep_summary_stats:
    output:
        expand(join(SUMMSTATS_DIR, "{chrom}"), trait=traits, chrom=chroms)
    shell:
        """
        # Rscript scripts/garfield/prep_T1D_chiou_2021.R
        # Rscript scripts/garfield/prep_FGluadjBMI_MAGIC.R  
        # Rscript scripts/garfield/prep_FInsadjBMI_MAGIC.R  
        # Rscript scripts/garfield/prep_T2D_DIAMANTE.R
        """

rule prep_annotations:
    """
    Needs around 50-60G memory
    """
    input:
        trait = lambda wildcards: expand(
            join(SUMMSTATS_DIR, "{chrom}"), trait=wildcards.trait, chrom=chroms
        ),
        # features = summits + dars,
        features = dars,
    output:
        chrom = expand(join(ANNOT_DIR, "{chrom}"), chrom=chroms, trait="{trait}"),
        link = join(ANNOT_DIR, "link_file.txt"),
    params:
        features = join(ANNOT_DIR, "features.txt"),
        sumstat = lambda wildcards: SUMMSTATS_DIR.format(trait=wildcards.trait),
        outdir = ANNOT_DIR,
    shell:
        """
        for i in {input.features}; do echo "${{i}}";done > {params.features}
        Rscript scripts/garfield/prepare_garfield_input.R {params.sumstat} \
            {params.outdir} {params.features}
        """

rule garfield_prep_chr:
    input:
        ptags = join(GARFIELD_DATA, "tags", "r01", "{chrom}"),
        ctags = join(GARFIELD_DATA, "tags", "r08", "{chrom}"),
        maftss = join(GARFIELD_DATA, "maftssd", "{chrom}"),
        pvalue = join(SUMMSTATS_DIR, "{chrom}"),
        annot = join(ANNOT_DIR, "{chrom}"),
    output:
        join(PREP_CHR_DIR, "{chrom}"),
    params:
        bin = join(GARFIELD_SW, "garfield-prep-chr"),
        chr = "{chrom}",
    shell:
        """
        {IONICE} {params.bin} -ptags {input.ptags} -ctags {input.ctags} \
            -maftss {input.maftss} -pval {input.pvalue} -ann {input.annot} \
            -chr {params.chr} -o {output}
        """

rule combined_prep_results:
    input:
        expand(
            join(PREP_CHR_DIR, "{chrom}"),
            trait="{trait}", chrom=chroms
        ),
    output:
        join(PREP_CHR_DIR, "prepfile.txt"),
    params:
        indir = PREP_CHR_DIR
    shell:
        """
        cat {params}/chr1 > {output} && \
        for i in {{2..22}}; do cat {params}/chr${{i}}; done >> {output}
        """

rule garfield_Meff_Padj:
    input:
        rules.combined_prep_results.output,
    output:
        join(OUTDIR, "meff_padj.out.txt"),
    params:
        R = join(GARFIELD_SW, "garfield-Meff-Padj.R")
    shell:
        """
        {IONICE} {params.R} -i {input} -o {output}
        """

rule garfield_test_all:
    input:
        link = rules.prep_annotations.output.link,
        prep = rules.combined_prep_results.output,
    output:
        join(OUTDIR, "garfield.c0.out.txt"),
    params:
        R = join(GARFIELD_SW, "garfield-test.R"),
        pt = "-pt 1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8",
        bins = "-b n5,m5,t5",
    shell:
        """
        {IONICE} Rscript {params.R} -i {input.prep} -o {output} \
            -l {input.link} {params.pt} {params.bins} -c 0
        """