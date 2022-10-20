#!/usr/bin/env python3

from os.path import join


def list_gwas():
    for gwas in config["GWAS"].keys():
        yield gwas


# Config
unique_id = config["liger_obj"]
CLUSTERS_N = config["n_clusters"]
MAF = config["maf"]
SNPS = config["snps"]
EXT = config["summit_extension"]
OTHER_ANNOTATIONS = config["other_annots"]

T1D = config["GWAS"]["T1D"]
T1D2020 = config["GWAS"]["T1D2020"]
T1D2021 = config["GWAS"]["T1D2021"]


# Directories
BIM_DIR = config["plink_dir"]
BASELINE_DIR = config["baseline_dir"]
WEIGHTS_DIR = config["weights_dir"]
FRQ_DIR = config["freq_dir"]
WHM3 = config["w_hm3"]

BASE_DIR = join(config["basedir"], "clusters", unique_id) 
# RESULTS = join(BASE_DIR, "ldsc_maf0.01_min2samples")
RESULTS = join(BASE_DIR, "ldsc_maf0.01_summits" + EXT)
MACS2_DIR = join(BASE_DIR, "macs2")
GWAS_DIR = join(RESULTS, "gwas")
ANNOT_DIR = join(RESULTS, "annotations")
SUMSTATS_DIR = join(RESULTS, "sumstats")
H2_DIR = join(RESULTS, "partitioned_heritability")



# Wildcards
clusters = ["GCG", "INS", "KRT19", "PDGFRA", "PPY", "PRSS1", "RGS5", 
            "SDS", "SST", "VWF"]
chroms = [x for x in range(1, 23)]
gwas_list = list_gwas()


rule all:
    input:
        expand(
            join(H2_DIR, "cluster_{cluster}__{gwas}.results"),
            cluster=clusters, gwas=gwas_list
        )



rule exists_bim:
    output:
        join(BIM_DIR, "1000G.EUR.QC.{chrom}.bim"),

rule make_annotations:
    input:
        bim = rules.exists_bim.output,
        # bed = join(MACS2_DIR, "{cluster}_peaks.min2.bed")
        bed = join(MACS2_DIR, "{cluster}_summits.ext" + EXT + ".min2.bed"),
    output:
        join(ANNOT_DIR, "cluster_{cluster}__{chrom}.annot.gz")
    shell:
        """
        scripts/make_annot_custom.py --bimfile {input.bim} \
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
        scripts/prep-aylward-ldsc.py {input} > {output}
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
        cp {input.t2d} {output.t2d}
        cp {input.fglu} {output.fglu}
        cp {input.fins} {output.fins}
        """

rule calculate_ld_scores:
    input:
        bim = rules.exists_bim.output,
        annot = rules.make_annotations.output,
        snps = SNPS,
    output:
        M = join(ANNOT_DIR, "cluster_{cluster}__{chrom}.l2.M"),
        M_5_50 = join(ANNOT_DIR, "cluster_{cluster}__{chrom}.l2.M_5_50"),
        ldscore = join(ANNOT_DIR, "cluster_{cluster}__{chrom}.l2.ldscore.gz"),
    params:
        handle = join(ANNOT_DIR, "cluster_{cluster}__{chrom}"),
        bim_handle = join(BIM_DIR, "1000G.EUR.QC.{chrom}")
    shell:
        """
        export MKL_NUM_THREADS=1 && \
        export NUMEXPR_NUM_THREADS=1 && \
        export OPENBLAS_NUM_THREADS=1 && \
        export OMP_NUM_THREADS=1 && \
        export VECLIB_MAXIMUM_THREADS=1 && \
        ldsc.py --l2 --thin-annot --bfile {params.bim_handle} --maf {MAF} \
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
    shell:
        """
        export MKL_NUM_THREADS=1 && \
        export NUMEXPR_NUM_THREADS=1 && \
        export OPENBLAS_NUM_THREADS=1 && \
        export OMP_NUM_THREADS=1 && \
        export VECLIB_MAXIMUM_THREADS=1 && \
        munge_sumstats.py --sumstats {input.gwas} \
            --merge-alleles {input.snps} \
            --out {params.out_handle}
        """

rule process_ldsc:
    input:
        gwas = rules.sumstats.output,
        annot = lambda wildcards: expand(
            join(ANNOT_DIR, "cluster_{}".format(wildcards.cluster) + 
                            "__{chrom}.l2.M"),
            chrom=chroms
        ),
    output:
        join(H2_DIR, "cluster_{cluster}__{gwas}.results"),
    params:
        baseline_handle = join(BASELINE_DIR, "baselineLD."),
        clus_handle = join(ANNOT_DIR, "cluster_{cluster}__"),
        w_handle = join(WEIGHTS_DIR, "weights.hm3_noMHC."),
        frq_handle = join(FRQ_DIR, "1000G.EUR.QC."),
        out_handle = join(H2_DIR, "cluster_{cluster}__{gwas}"),
    shell:
        """
        export MKL_NUM_THREADS=1 && \
        export NUMEXPR_NUM_THREADS=1 && \
        export OPENBLAS_NUM_THREADS=1 && \
        export OMP_NUM_THREADS=1 && \
        export VECLIB_MAXIMUM_THREADS=1 && \
        ldsc.py --h2 {input.gwas} --maf {MAF} \
            --ref-ld-chr {params.baseline_handle},{params.clus_handle} \
            --w-ld-chr {params.w_handle} \
            --overlap-annot --frqfile-chr {params.frq_handle} \
            --out {params.out_handle} --print-coefficients
        """