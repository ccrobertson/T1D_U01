#!/usr/bin/env python3

from os.path import join


def list_gwas():
    for gwas in config["GWAS"].keys():
        yield gwas


def get_annotations_bed(wildcards):
    feature = wildcards.feature
    if any(x in feature for x in ["CVB4","Cytokine"]):
        outdir = DA_DIR_PERT
    else:
        outdir = DA_DIR_T1D
    out = join(outdir, "{feature}.bed")
    return out


# Config
unique_id = config["liger_obj"]
MAF = config["maf"]
SNPS = config["snps"]
OTHER_ANNOTATIONS = config["other_annots"]
T1D = config["GWAS"]["T1D"]
T1D2020 = config["GWAS"]["T1D2020"]
T1D2021 = config["GWAS"]["T1D2021"]
CLUSTERS_N = config["n_clusters"]
IONICE = "ionice -c2 -n7"

# Directories
BIM_DIR = config["plink_dir"]
BASELINE_DIR = config["baseline_dir"]
WEIGHTS_DIR = config["weights_dir"]
FRQ_DIR = config["freq_dir"]
WHM3 = config["w_hm3"]
BASE_DIR = join(config["basedir"], "clusters", unique_id) 

# Broadpeaks or summits?

# RESULTS_handle = "ldsc_dars"
# DA_DIR = join(BASE_DIR, config["da_peaks"])
RESULTS_handle = "ldsc_dars_henry_random_effect"
DA_DIR_PERT = join(BASE_DIR, config["da_summits_pert"])
DA_DIR_T1D = join(BASE_DIR, config["da_summits_t1d"])


RESULTS = join(BASE_DIR, RESULTS_handle, config["outdir"])
GWAS_DIR = join(BASE_DIR, RESULTS_handle, "gwas")
ANNOT_DIR = join(RESULTS, "annotations")
SUMSTATS_DIR = join(RESULTS, "sumstats")
H2_DIR = join(RESULTS, "partitioned_heritability")


# Wildcards
chroms = [x for x in range(1, 23)]
conds = config["conditions"]

# features, = glob_wildcards(
#     join(DA_DIR.format(cluster=3), "{feature}.bed")
# )
# gwas_list = list_gwas()

clusters = config["clusters"]

features = expand(
    "{cond}_dar_{direction}_{fc}pct",
    cond = conds,
    direction = ["above", "below"],
    fc = [0.8, 0.85, 0.9, 0.95, 0.99]
)

# features = expand(
#     "{cond}_dar_{direction}_{fc}FDR",
#     cond = ["T1D", "preT1D", "CVB4", "Cytokine"],
#     direction = ["sig", "nsig"],
#     fc = [0.1, 0.05, 0.01, 0.001]
# )
gwas_list = ["T1D2021"]
# clusters = [3]

print(DA_DIR_T1D)
print(DA_DIR_PERT)
print(features)



rule all:
    input:
        expand(
            join(H2_DIR, "{feature}__{gwas}.results"),
            cluster=clusters, feature=features, gwas=gwas_list
        )



rule exists_bim:
    output:
        join(BIM_DIR, "1000G.EUR.QC.{chrom}.bim"),

rule make_annotations:
    input:
        bim = rules.exists_bim.output,
        bed = get_annotations_bed,
    output:
        join(ANNOT_DIR, "{feature}__{chrom}.annot.gz")
    shell:
        """
        {IONICE} scripts/make_annot_custom.py --bimfile {input.bim} \
            --bed-files {input.bed} {OTHER_ANNOTATIONS} \
            --annot-file {output}
        """

rule prep_gwas_t1d:
    input:
        T1D,
    output:
        join(GWAS_DIR, "T1D.stats.txt")
    shell:
        """
        {IONICE} scripts/prep-aylward-ldsc.py {input} > {output}
        """

rule prep_gwas_t1d_2020:
    input:
        T1D2020,
    output:
        join(GWAS_DIR, "T1D2020.stats.txt")
    shell:
        """
        zcat {input} > {output}
        """

rule prep_gwas_t1d_2021:
    input:
        T1D2021,
    output:
        join(GWAS_DIR, "T1D2021.stats.txt")
    shell:
        """
        zcat {input} > {output}
        """

rule prep_gwas_others:
    input:
        t2d = config["GWAS"]["T2Dbmiadj"],
        fglu = config["GWAS"]["FGluadjBMI"],
        fins = config["GWAS"]["FInsadjBMI"],
    output:
        t2d = join(GWAS_DIR, "T2Dbmiadj.stats.txt"),
        fglu = join(GWAS_DIR, "FGluadjBMI.stats.txt"),
        fins = join(GWAS_DIR, "FInsadjBMI.stats.txt"),
    shell:
        """
        {IONICE} cp {input.t2d} {output.t2d}
        {IONICE} cp {input.fglu} {output.fglu}
        {IONICE} cp {input.fins} {output.fins}
        """

rule calculate_ld_scores:
    input:
        bim = rules.exists_bim.output,
        annot = rules.make_annotations.output,
        snps = SNPS,
    output:
        M = join(ANNOT_DIR, "{feature}__{chrom}.l2.M"),
        M_5_50 = join(ANNOT_DIR, "{feature}__{chrom}.l2.M_5_50"),
        ldscore = join(ANNOT_DIR, "{feature}__{chrom}.l2.ldscore.gz"),
    params:
        handle = join(ANNOT_DIR, "{feature}__{chrom}"),
        bim_handle = join(BIM_DIR, "1000G.EUR.QC.{chrom}")
    threads: 8
    shell:
        """
        export MKL_NUM_THREADS=1 && \
        export NUMEXPR_NUM_THREADS=1 && \
        export OPENBLAS_NUM_THREADS=1 && \
        export OMP_NUM_THREADS=1 && \
        export VECLIB_MAXIMUM_THREADS=1 && \
        {IONICE} ldsc.py --l2 --thin-annot --bfile {params.bim_handle} --maf {MAF} \
            --ld-wind-cm 1 --print-snps {input.snps}  --annot {input.annot} \
            --out {params.handle}
        """

rule sumstats:
    input:
        gwas = join(GWAS_DIR, "{gwas}.stats.txt"),
        snps = WHM3,
    output:
        join(SUMSTATS_DIR, "{gwas}.sumstats.gz")
    params:
        out_handle = join(SUMSTATS_DIR, "{gwas}")
    threads: 8
    shell:
        """
        export MKL_NUM_THREADS=1 && \
        export NUMEXPR_NUM_THREADS=1 && \
        export OPENBLAS_NUM_THREADS=1 && \
        export OMP_NUM_THREADS=1 && \
        export VECLIB_MAXIMUM_THREADS=1 && \
        {IONICE} munge_sumstats.py --sumstats {input.gwas} \
            --merge-alleles {input.snps} \
            --out {params.out_handle}
        """

rule process_ldsc:
    input:
        gwas = rules.sumstats.output,
        annot = lambda wildcards: expand(
            join(ANNOT_DIR.format(cluster=wildcards.cluster), 
                "{}".format(wildcards.feature) + "__{chrom}.l2.M"),
            chrom=chroms
        ),
    output:
        join(H2_DIR, "{feature}__{gwas}.results"),
    params:
        baseline_handle = join(BASELINE_DIR, "baselineLD."),
        annot_handle = join(ANNOT_DIR, "{feature}__"),
        w_handle = join(WEIGHTS_DIR, "weights.hm3_noMHC."),
        frq_handle = join(FRQ_DIR, "1000G.EUR.QC."),
        out_handle = join(H2_DIR, "{feature}__{gwas}"),
    threads: 8
    shell:
        """
        export MKL_NUM_THREADS=1 && \
        export NUMEXPR_NUM_THREADS=1 && \
        export OPENBLAS_NUM_THREADS=1 && \
        export OMP_NUM_THREADS=1 && \
        export VECLIB_MAXIMUM_THREADS=1 && 
        {IONICE} ldsc.py --h2 {input.gwas} --maf {MAF} \
            --ref-ld-chr {params.baseline_handle},{params.annot_handle} \
            --w-ld-chr {params.w_handle} \
            --overlap-annot --frqfile-chr {params.frq_handle} \
            --out {params.out_handle} --print-coefficients
        """