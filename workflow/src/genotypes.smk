#!/usr/bin/env python3


from os.path import join
import os
from functools import partial
import requests
import json


name = config["name"]
autosomes = list(range(1, 23))

_results = partial(os.path.join, "results", name)
_logs = partial(_results, "logs")

def _source(source, file):
    f = os.path.join(config['sources'][source], file)
    return(f)

def get_donors_string():
    return(','.join(config["donors"]))

def get_sources():
    return(config["sources"].keys())

def n_sources():
    return(len(config["sources"].keys()))


rule all:
     input:
        #expand(_results("imputation_results/{source}/chr{chr}.donors_only.dose.vcf.gz"), chr=20, source=['nih_20220216','6719-NM']),
        #expand(_results("imputation_results_combined/chr20.donors_only.{filter}.dose.vcf.gz"), filter="maf_gt_1pct__rsq_gt_pt95"),
        expand(_results("imputation_results_combined/all_chromosomes.donors_only.{filter}.bed"), filter=['maf_gt_1pct__rsq_gt_pt95','maf_gt_1pct__rsq_gt_pt30']),


### NOTE: rsq filter did not work in TOPMED imputation
rule filter_variants_round1:
    input:
        vcf = lambda wildcards: _source(wildcards.source,"chr{chr}.dose.vcf.gz"),
    output:
        vcf = _results("imputation_results/{source}/chr{chr}.maf_gt_1pct__rsq_gt_pt30.dose.vcf.gz"),
    conda:
        "genetics"
    shell:
        """
        bcftools view -i 'MAF>.01 & R2>0.3' --threads 10 -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """

rule extract_donors:
    input:
        vcf = _results("imputation_results/{source}/chr{chr}.maf_gt_1pct__rsq_gt_pt30.dose.vcf.gz"),
    output:
        vcf = _results("imputation_results/{source}/chr{chr}.donors_only.maf_gt_1pct__rsq_gt_pt30.dose.vcf.gz"),
    params:
        donors = get_donors_string(),
    conda:
        "genetics"
    shell:
        """
        bcftools view --samples {params.donors} --force-samples --threads 10 -Oz -o {output.vcf} {input.vcf}
	    tabix -p vcf {output.vcf}
        """


### stringent rsq for identity matching
rule filter_variants_round2:
    input:
        vcf = _results("imputation_results/{source}/chr{chr}.donors_only.maf_gt_1pct__rsq_gt_pt30.dose.vcf.gz"),
    output:
        vcf = _results("imputation_results/{source}/chr{chr}.donors_only.maf_gt_1pct__rsq_gt_pt95.dose.vcf.gz"),
    conda:
        "genetics"
    shell:
        """
        bcftools view -i 'R2>0.95' --threads 10 -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """

rule merge_vcfs:
    input:
        vcfs = expand(_results("imputation_results/{source}/chr{{chr}}.donors_only.{{filter}}.dose.vcf.gz"), source=get_sources()),
    output:
        tmpfiles = expand(_results("imputation_results_combined/bcftools_isec.{{filter}}.chr{{chr}}/000{k}.vcf.gz"), k=list(range(0,n_sources()))),
        vcf = _results("imputation_results_combined/chr{chr}.donors_only.{filter}.dose.vcf.gz"),
    params:
        tmpdir = _results("imputation_results_combined/bcftools_isec.{filter}.chr{chr}"),
        nsources = n_sources(),
    conda:
        "genetics"
    shell:
        """
        bcftools isec -p {params.tmpdir} -n {params.nsources} -c none -Oz {input.vcfs}
        bcftools merge -Oz -o {output.vcf} {output.tmpfiles}
        """

rule combine_chromosomes:
    input:
        vcfs = expand(_results("imputation_results_combined/chr{chr}.donors_only.{{filter}}.dose.vcf.gz"), chr=autosomes),
    output:
        vcf = _results("imputation_results_combined/all_chromosomes.donors_only.{filter}.dose.vcf.gz"),
        tbi = _results("imputation_results_combined/all_chromosomes.donors_only.{filter}.dose.vcf.gz.tbi"),
    conda:
        "genetics"
    shell:
        """
        bcftools concat -O z -o {output.vcf} {input.vcfs}
        tabix -p vcf {output.vcf}
        """

rule convert_to_plink:
    input:
        vcf = _results("imputation_results_combined/all_chromosomes.donors_only.{filter}.dose.vcf.gz"),
    output:
        bed = _results("imputation_results_combined/all_chromosomes.donors_only.{filter}.bed"),
    params:
        prefix = _results("imputation_results_combined/all_chromosomes.donors_only.{filter}"),
    conda:
        "genetics"
    shell:
        """
        plink --vcf {input.vcf} --double-id --make-bed --out {params.prefix}
        """
