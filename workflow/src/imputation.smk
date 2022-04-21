#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import requests
import json


_results = partial(os.path.join, config["results"])
_source = partial(os.path.join, config["source"])
#_resources = partial(os.path.join, config["resources"])
_logs = partial(_results, "logs")

chromosomes = list(range(1, 23))

def get_imputation_input(chr):
    vcf = _source("*chr{chr}.vcf.gz")

rule all:
    input:
        expand(_results("imputation_input/chr{chr}.vcf.gz"), chr=chromosomes),


rule symlinks:
    input:
        target = get_imputation_input(chr),
    output:
        linkname = _results("imputation_input/chr{chr}.vcf.gz"),
    shell:
        """
        ln -s {input.target} {output.linkname}
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
