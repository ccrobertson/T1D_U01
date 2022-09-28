#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import pandas as pd
import math


name = config["name"]

_results = partial(os.path.join, "results", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


configfile: "data/nandini_run_design.json"
configfile: "workflow/src/demultiplex_samples.yaml"
samples = config["samples"].keys()
batches = config["batches"].keys()

### restrict wildcards to existing runs, batches, and samples
wildcard_constraints:
   sample="|".join(samples),
   batch="|".join(batches),
   bamtype="unfiltered|auto|p01",

IONICE = 'ionice -c2 -n7'


def genotyped_subjects_by_batch_as_string(sample, genotyped_donors):
    samplename = config["samples"][sample]['batch']
    donors_in_batch = config["batches"][samplename].keys()
    donors_in_vcf = pd.read_table(genotyped_donors, header=None).iloc[:,0].tolist()
    overlap = ' '.join(sorted(list(set(donors_in_batch) & set(donors_in_vcf))))
    return overlap

def genotyped_subjects_by_batch_as_string2(sample, genotyped_donors):
    samplename = config["samples"][sample]['batch']
    donors_in_batch = config["batches"][samplename].keys()
    donors_in_vcf = pd.read_table(genotyped_donors, header=None).iloc[:,0].tolist()
    overlap = ','.join(sorted(list(set(donors_in_batch) & set(donors_in_vcf))))
    return overlap

def count_genotyped_subjects_by_batch(sample, genotyped_donors):
    samplename = config["samples"][sample]['batch']
    donors_in_batch = config["batches"][samplename].keys()
    donors_in_vcf = pd.read_table(genotyped_donors, header=None).iloc[:,0].tolist()
    overlap = list(set(donors_in_batch) & set(donors_in_vcf))
    return len(overlap)

def count_subjects_by_batch(samplename):
    return len(config["batches"][samplename].keys())

def calculate_max_missing_by_batch(sample, genotyped_donors):
    samplename = config["samples"][sample]['batch']
    n_subjects = len(config["batches"][samplename].keys())
    n_donors = pd.read_table(genotyped_donors, header=None).shape[0]
    max_missing_rate = math.ceil((n_subjects-1)/(n_donors+n_subjects)*100, )/100
    return max_missing_rate

def iterate_donors_by_sample(sample):
    batch = config["samples"][sample]['batch']
    return config["batches"][batch].keys()


rule all:
    input:
        #~~~~~~~~ demuxlet
        expand(_results("demuxlet-unfiltered/5125-{batch}/demuxlet.best"), batch="NM-1"),
        #expand(_results("demuxlet-unfiltered/5124-{batch}/demuxlet.best"), batch=batches),
        #expand(_results("demuxlet-auto/5124-{batch}/demuxlet.best"), batch=batches),
        #expand(_results("demuxlet-p01/5124-{batch}/demuxlet.best"), batch=batches),
        #~~~~~~~~ souporcell
        #expand(_results("check_by_sample-unfiltered/5125-{batch}/corrplot_r_by_donor.png"), batch="NM-1"),
        #expand(_results("check_by_sample-unfiltered/5124-{batch}/corrplot_r_by_donor.png"), batch=batches),
        #expand(_results("check_by_sample-auto/5124-{batch}/corrplot_r_by_donor.png"), batch=batches),
        #expand(_results("check_by_sample-p01/5124-{batch}/corrplot_r_by_donor.png"), batch=batches),
        #summary
        #_results("check-unfiltered/corrplot_r_by_donor.png"),



#DEMUXLET2: https://github.com/statgen/popscle
#NOTE: example for using demuxlet "v2":
#https://github.com/porchard/Multiome-QC-NextFlow/blob/master/qc.nf#L114
#popscle demuxlet --sam $bam
# --vcf $vcf
# --alpha 0
# --alpha 0.5
# --group-list $barcodes
# --field GT
# --out ${library}-${modality}
rule demuxlet:
    input:
        bam = lambda wildcards: config['libraries'][wildcards.sample]["bam_"+wildcards.bamtype],
        barcodes =
        vcf = config["donor_genotypes"],
    output:
        _results("demuxlet-{bamtype}/{sample}/demuxlet.best"),
    params:
        prefix = _results("demuxlet-{bamtype}/{sample}/demuxlet"),
        donors = lambda wildcards: genotyped_subjects_by_batch_as_string2(wildcards.sample, config["genotyped_donors"]),
    container:
        "workflow/envs/demuxlet_20220204.sif",
    shell:
        """
        echo {params.donors} | sed 's/,/\\n/g' > {params.prefix}.donor_list.txt
        popscle demuxlet --sam {input.bam} \
            --tag-group CB \
            --tag-UMI UB \
            --vcf {input.vcf} \
            --field GP \
            --sm-list {params.prefix}.donor_list.txt \
            --group-list {input.barcodes} \
            --out {params.prefix}
        """


rule define_doublets:
    input:
        demuxlet_best = _results("demuxlet-{bamtype}/{sample}/demuxlet.best"),
    output:
        doublets = _results("demuxlet-{bamtype}/{sample}/doublets.txt"),
    shell:
        """
        awk '$6~/^DBL/ {print $1}' {input.demuxlet_best} > {output.doublets}
        """


# Issue about exploring range of k
# https://github.com/wheaton5/souporcell/issues/7
#
# Issues about what happens when donors are missing from VCF
# https://github.com/wheaton5/souporcell/issues/77
# https://github.com/wheaton5/souporcell/issues/141
rule souporcell:
    input:
        bam = lambda wildcards: config['libraries'][wildcards.sample]["bam_"+wildcards.bamtype],
        barcodes = lambda wildcards: config['libraries'][wildcards.sample]['barcodes'],
        fasta = _resources("hg38/hg38_cvb4.fa"),
        common_variants = _resources("common_variants_grch38_fixed.vcf"),
        donor_genotypes = config["donor_genotypes_unzipped"],
    output:
        cluster_genotypes = _results("souporcell-{bamtype}/{sample}/cluster_genotypes.vcf"),
    params:
        outdir = _results("souporcell-{bamtype}/{sample}"),
        k = lambda wildcards: count_genotyped_subjects_by_batch(wildcards.sample, config["genotyped_donors"]),
        donor_list_as_string = lambda wildcards: genotyped_subjects_by_batch_as_string(wildcards.sample, config["genotyped_donors"]),
        threads=10,
        sample = "{sample}",
    shell:
        """
        singularity exec workflow/envs/souporcell_latest.sif souporcell_pipeline.py \
            -i {input.bam} \
            -b {input.barcodes} \
            -f {input.fasta} \
            -t {params.threads} \
            --cluster {params.k} \
            --known_genotypes {input.donor_genotypes} \
            --known_genotypes_sample_names {params.donor_list_as_string} \
            --skip_remap True \
            -o {params.outdir}
        """



rule souporcell_fix_vcf:
    input:
        cluster_genotypes = _results("souporcell-{bamtype}/{sample}/cluster_genotypes.vcf"),
    output:
        cluster_genotypes_reformatted = _results("souporcell-{bamtype}/{sample}/cluster_genotypes_reformatted.vcf.gz"),
    params:
        sample = "{sample}",
    conda:
        "genetics"
    shell:
        """
        awk -v OFS='\t' -v sample={params.sample} '$1!~/^#CHROM/ {{print $0}} $1~/^#CHROM/ {{for(i=10; i<=NF; ++i) $i=sample"_souporcell_"$i; print $0 }}' {input.cluster_genotypes} | grep -v BACKGROUND | bgzip > {output.cluster_genotypes_reformatted}
        tabix -p vcf {output.cluster_genotypes_reformatted}
        """

rule demultiplex_by_sample:
    input:
        souporcell_vcf = _results("souporcell-{bamtype}/{sample}/cluster_genotypes_reformatted.vcf.gz"),
        imputed_vcf = config["donor_genotypes"],
    output:
        kin0 =  _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.raw"),
    params:
        prefix = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors"),
        prefix_filtered = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered"),
        max_missing = lambda wildcards: calculate_max_missing_by_batch(wildcards.sample, config["genotyped_donors"])
    conda:
        "genetics"
    shell:
        """
        #merge souporcell and donor vcfs
        bcftools merge -O z -o {params.prefix}.vcf.gz {input.souporcell_vcf} {input.imputed_vcf}
        tabix -p vcf {params.prefix}.vcf.gz

        #convert merged vcf to plink
        plink --vcf {params.prefix}.vcf.gz --double-id --make-bed --out {params.prefix}

        #filter merged vcf for missingness
        plink --bfile {params.prefix} --geno {params.max_missing} --make-bed --out {params.prefix_filtered}

        #run king relationship inference
        plink2 --bfile {params.prefix_filtered} --make-king square --make-king-table --out {params.prefix_filtered}

        #recode as alt allele count matrix
        plink --bfile {params.prefix_filtered} --recodeA --out {params.prefix_filtered}
        """


rule demultiplex_by_sample_plots:
    input:
        kin0 = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.raw"),
    output:
        _results("check_by_sample-{bamtype}/{sample}/corrplot_r_by_donor.png"),
    params:
        outdir = _results("check_by_sample-{bamtype}/{sample}"),
        donors = "HPAP105,HPAP093,ICRH139,ICRH142,ICRH143,HPAP107",
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/demultiplexing_plots.R \
            --kin0 {input.kin0} \
            --kinMat {input.kinMat} \
            --kinIds {input.kinIds} \
            --donors "HPAP105,HPAP093,ICRH139,ICRH142,ICRH143,HPAP107" \
            --genoraw {input.genoraw} \
            --outdir {params.outdir}
        """


rule demultiplex:
    input:
        souporcell_vcfs = expand(_results("souporcell-{{bamtype}}/{sample}/cluster_genotypes_reformatted.vcf.gz"), sample=samples),
        imputed_vcf = config["donor_genotypes"],
    output:
        kin0 =  _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat =  _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.king"),
        kinIds =  _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.raw"),
    params:
        prefix = _results("check-{bamtype}/souporcell_clusters_and_donors"),
        prefix_filtered = _results("check-{bamtype}/souporcell_clusters_and_donors_filtered"),
        max_missing = 0.5,
    conda:
        "genetics"
    shell:
        """
        #merge souporcell and donor vcfs
        bcftools merge -O z -o {params.prefix}.vcf.gz {input.souporcell_vcfs} {input.imputed_vcf}
        tabix -p vcf {params.prefix}.vcf.gz

        #convert merged vcf to plink
        plink --vcf {params.prefix}.vcf.gz --double-id --make-bed --out {params.prefix}

        #filter merged vcf for missingness
        plink --bfile {params.prefix} --geno {params.max_missing} --make-bed --out {params.prefix_filtered}

        #run king relationship inference
        plink2 --bfile {params.prefix_filtered} --make-king square --make-king-table --out {params.prefix_filtered}

        #recode as alt allele count matrix
        plink --bfile {params.prefix_filtered} --recodeA --out {params.prefix_filtered}
        """

rule demultiplex_plots:
    input:
        kin0 = _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("check-{bamtype}/souporcell_clusters_and_donors_filtered.raw"),
    output:
        _results("check-{bamtype}/corrplot_r_by_donor.png"),
    params:
        outdir = _results("check-{bamtype}"),
        donors = "HPAP105,HPAP093,ICRH139,ICRH142,ICRH143,HPAP107",
    shell:
        """
        Rscript workflow/scripts/demultiplexing_plots.R \
            --kin0 {input.kin0} \
            --kinMat {input.kinMat} \
            --kinIds {input.kinIds} \
            --donors {params.donors} \
            --genoraw {input.genoraw} \
            --outdir {params.outdir}
        """
