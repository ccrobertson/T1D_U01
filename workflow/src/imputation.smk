#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import requests
import json


name = config["name"]

_results = partial(os.path.join, "results/imputation", name)
_source = partial(os.path.join, config["source"])
_logs = partial(_results, "logs")

chromosomes = list(range(1, 23))


# def iterate_subjects_by_batch(samplename):
#      return config["batches"][samplename].keys()
#
def get_donorstring(donors):
    return(','.join(donors))

rule all:
    input:
        #_results("imputation_input/chrALL.vcf.gz"),
        _results("imputation_results/chrALL.donors_only.dose.vcf.gz"),
        _results("imputation_results/chrALL.donors_only.maf_gt_0.01.dose.vcf.gz"),
        _results("imputation_results/chrALL.maf_gt_0.01.dose.vcf"),


rule combine_chromosomes:
    input:
        vcf_by_chr = expand(_results("imputation_results/chr{chr}.dose.vcf.gz"), chr=chromosomes),
    output:
        vcf_combined = _results("imputation_results/chrALL.dose.vcf.gz"),
        tbi = _results("imputation_results/chrALL.dose.vcf.gz.tbi"),
    shell:
        """
        bcftools concat -O z -o {output.vcf_combined} {input.vcf_by_chr}
        tabix -p vcf {output.vcf_combined}
        """

rule remove_1kg:
    input:
        vcf_combined = _results("imputation_results/chrALL.dose.vcf.gz"),
        donors = _results("donors.txt"),
    output:
        vcf_donors = _results("imputation_results/chrALL.donors_only.dose.vcf.gz"),
        tbi = _results("imputation_results/chrALL.donors_only.dose.vcf.gz.tbi"),
    shell:
        """
        bcftools view --samples-file {input.donors} -O z -o {output.vcf_donors} {input.vcf_combined}
	    tabix -p vcf {output.vcf_donors}
        """

rule filter_variants:
    input:
        vcf = _results("imputation_results/chrALL.donors_only.dose.vcf.gz"),
    output:
        vcf_filtered = _results("imputation_results/chrALL.donors_only.maf_gt_0.01.dose.vcf"),
    conda:
        "genetics"
    shell:
        """
        vcftools --gzvcf {input.vcf} --maf 0.01 --recode --stdout > {output.vcf_filtered}
        bgzip -c {output.vcf_filtered} > {output.vcf_filtered}.gz
        tabix -p vcf {output.vcf_filtered}.gz
        """

rule filter_variants_1kg:
    input:
        vcf = _results("imputation_results/chrALL.dose.vcf.gz"),
    output:
        vcf_filtered = _results("imputation_results/chrALL.maf_gt_0.01.dose.vcf"),
    conda:
        "genetics"
    shell:
        """
        vcftools --gzvcf {input.vcf} --maf 0.01 --recode --stdout > {output.vcf_filtered}
        bgzip -c {output.vcf_filtered} > {output.vcf_filtered}.gz
        tabix -p vcf {output.vcf_filtered}.gz
        """


# imputation server url
# url = 'https://imputationserver.sph.umich.edu/api/v2'
#
# # add token to header (see Authentication)
# headers = {'X-Auth-Token' : token }
#
# # submit new job
# vcf = '/path/to/genome.vcf.gz';
# files = {'input-files' : open(vcf, 'rb')}
# r = requests.post(url + "/jobs/submit/minimac4", files=files, headers=headers)
# if r.status_code != 200:
#     raise Exception('POST /jobs/submit/minimac4 {}'.format(r.status_code))
#
# # print message
# print r.json()['message']
# print r.json()['id']
