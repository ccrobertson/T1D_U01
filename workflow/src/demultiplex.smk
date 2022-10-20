#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import pandas as pd
import math
import glob


name = config["name"]

_results = partial(os.path.join, "results", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


### Batch design and library key
configfile: "data/batch_design.json"
configfile: "data/library_key.json"

### Bam file and barcodes to use
configfile: "workflow/src/demultiplex_samples.yaml"

### Samples and
samples = config["samples"].keys()
batches = config["batches"].keys()
wildcard_constraints:
   sample="|".join(samples),
   batch="|".join(batches),
   bamtype="unfiltered|auto|p01|p10",



IONICE = 'ionice -c2 -n7'

def sample_to_batch(sample):
    return config["libraries"][sample]["Batch"]

def batch_to_donors(batch):
    return config["batches"][batch]["donors"].keys()

def sample_to_donors(sample):
    batch = sample_to_batch(sample)
    donors = batch_to_donors(batch)
    return donors

def count_donors(sample):
    return len(sample_to_donors(sample))

def genotyped_donors_by_sample_as_list(sample, genotyped_donors):
    donors_in_sample = sample_to_donors(sample)
    donors_in_vcf = pd.read_table(genotyped_donors, header=None).iloc[:,0].tolist()
    overlap = sorted(list(set(donors_in_sample) & set(donors_in_vcf)))
    return overlap

def genotyped_donors_by_sample_as_string(sample, genotyped_donors):
    overlap = ' '.join(genotyped_donors_by_sample_as_list(sample, genotyped_donors))
    return overlap

def genotyped_donors_by_sample_as_string2(sample, genotyped_donors):
    overlap = ','.join(genotyped_donors_by_sample_as_list(sample, genotyped_donors))
    return overlap

def count_genotyped_donors_by_sample(sample, genotyped_donors):
    return len(genotyped_donors_by_sample_as_list(sample, genotyped_donors))

def calculate_max_missing_by_sample(sample, genotyped_donors):
    n_donors_pool = count_genotyped_donors_by_sample(sample, genotyped_donors)
    n_donors_vcf = pd.read_table(genotyped_donors, header=None).shape[0]
    max_missing_rate = math.ceil((n_donors_pool-1)/(n_donors_vcf+n_donors_pool)*100, )/100
    return max_missing_rate


def get_bam(sample, bamtype):
    bam = config["samples"][sample]["bam_"+bamtype]
    return bam

def get_barcodes(sample):
    barcodes = config["samples"][sample]["barcodes"]
    return barcodes

def get_chunks(sample):
    chunk_files = glob.glob(_results("barcodes_chunked/"+sample+"/barcodes.*.txt"))
    chunk_suffixes = [re.sub(".txt", "", re.sub("barcodes.", "", os.path.basename(x))) for x in chunk_files]
    return chunk_suffixes


rule all:
    input:
        #~~~~~~~~ chunk barcodes
        expand(_results("barcodes_chunked/{sample}"), sample=samples),
        #~~~~~~~~ demuxlet
        expand(_results("demuxlet-unfiltered/{sample}/demuxlet.best"), sample=samples),
        #expand(_results("demuxlet-unfiltered/{sample}/demuxlet.best"), sample=samples),
        #expand(_results("demuxlet-unfiltered/5124-{batch}/demuxlet.best"), batch=batches),
        #expand(_results("demuxlet-auto/5124-{batch}/demuxlet.best"), batch=batches),
        #expand(_results("demuxlet-p01/5124-{batch}/demuxlet.best"), batch=batches),
        #~~~~~~~~ souporcell
        #expand(_results("souporcell_check_by_sample-unfiltered/{sample}/corrplot_r_by_donor.png"), sample=samples),
        #expand(_results("souporcell_check_by_sample-unfiltered/5124-{batch}/corrplot_r_by_donor.png"), batch=batches),
        #expand(_results("souporcell_check_by_sample-auto/5124-{batch}/corrplot_r_by_donor.png"), batch=batches),
        #expand(_results("souporcell_check_by_sample-p01/5124-{batch}/corrplot_r_by_donor.png"), batch=batches),
        #summary
        #_results("souporcell_check-unfiltered/corrplot_r_by_donor.png"),


rule chunk_barcodes:
    input:
        barcodes = lambda wildcards: get_barcodes(wildcards.sample),
    output:
        directory(_results("barcodes_chunked/{sample}")),
    params:
        prefix = _results("barcodes_chunked/{sample}/barcodes."),
    shell:
        """
        mkdir {output}
        split --lines 500 --suffix-length 10 --additional-suffix=.txt {input.barcodes} {params.prefix}
        """



#DEMUXLET2: https://github.com/statgen/popscle
#NOTE: example for using demuxlet "v2":
#https://github.com/porchard/Multiome-QC-NextFlow/blob/master/qc.nf#L114
# MAJOR BUGS WHEN I RUN THIS VERSION:
# SNG.POSTERIOR is always 1
# BEST.POSTERIOR is always 1 for doublets and 0 for singlets
#popscle demuxlet --sam $bam
# --vcf $vcf
# --alpha 0
# --alpha 0.5
# --group-list $barcodes
# --field GT
# --out ${library}-${modality}
rule demuxlet:
    input:
        bam = lambda wildcards: get_bam(wildcards.sample, wildcards.bamtype),
        barcodes = _results("barcodes_chunked/{sample}/barcodes.{chunk}.txt"),
        vcf = config["donor_genotypes"],
    output:
        _results("demuxlet-{bamtype}/{sample}/chunked/demuxlet.{chunk}.best"),
    params:
        prefix = _results("demuxlet-{bamtype}/{sample}/chunked/demuxlet.{chunk}"),
        donors = lambda wildcards: genotyped_donors_by_sample_as_string2(wildcards.sample, config["genotyped_donors"]),
    container:
        "workflow/envs/demuxlet_20220204.sif",
    shell:
        """
        echo {params.donors} | sed 's/,/\\n/g' > {params.prefix}.donors.txt
        popscle demuxlet --sam {input.bam} \
            --tag-group CB \
            --tag-UMI UB \
            --vcf {input.vcf} \
            --field GP \
            --sm-list {params.prefix}.donors.txt \
            --group-list {input.barcodes} \
            --out {params.prefix}
        """

rule concat_demuxlet:
    input:
        first = lambda wildcards: expand(_results("demuxlet-unfiltered/{{sample}}/chunked/demuxlet.{chunk}.best"), chunk=get_chunks(wildcards.sample)[0]),
        bests = lambda wildcards: expand(_results("demuxlet-unfiltered/{{sample}}/chunked/demuxlet.{chunk}.best"), chunk=get_chunks(wildcards.sample)),
    output:
        best = _results("demuxlet-unfiltered/{sample}/demuxlet.best"),
    params:
        indir = _results("demuxlet-unfiltered/{sample}/chunked/"),
    shell:
        """
        awk 'NR==1' {input.first} > {output.best}
        for f in {input.bests}; do
            awk 'NR>1' $f >> {output.best}
        done
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
        bam = lambda wildcards: get_bam(wildcards.sample, wildcards.bamtype),
        barcodes = lambda wildcards: get_barcodes(wildcards.sample),
        fasta = _resources("hg38/hg38_cvb4.fa"),
        common_variants = _resources("common_variants_grch38_fixed.vcf"),
        donor_genotypes = config["donor_genotypes_unzipped"],
    output:
        cluster_genotypes = _results("souporcell-{bamtype}/{sample}/cluster_genotypes.vcf"),
    params:
        outdir = _results("souporcell-{bamtype}/{sample}"),
        k = lambda wildcards: count_genotyped_donors_by_sample(wildcards.sample, config["genotyped_donors"]),
        donor_list_as_string = lambda wildcards: genotyped_donors_by_sample_as_string(wildcards.sample, config["genotyped_donors"]),
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



rule souporcell_check_by_sample:
    input:
        souporcell_vcf = _results("souporcell-{bamtype}/{sample}/cluster_genotypes_reformatted.vcf.gz"),
        imputed_vcf = config["donor_genotypes"],
    output:
        kin0 =  _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.raw"),
    params:
        prefix = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors"),
        prefix_filtered = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered"),
        max_missing = lambda wildcards: calculate_max_missing_by_sample(wildcards.sample, config["genotyped_donors"])
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


rule souporcell_check_by_sample_plots:
    input:
        kin0 = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("souporcell_check_by_sample-{bamtype}/{sample}/souporcell_clusters_and_donors_filtered.raw"),
    output:
        _results("souporcell_check_by_sample-{bamtype}/{sample}/corrplot_r_by_donor.png"),
    params:
        outdir = _results("souporcell_check_by_sample-{bamtype}/{sample}"),
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


rule souporcell_check:
    input:
        souporcell_vcfs = expand(_results("souporcell-{{bamtype}}/{sample}/cluster_genotypes_reformatted.vcf.gz"), sample=samples),
        imputed_vcf = config["donor_genotypes"],
    output:
        kin0 =  _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat =  _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.king"),
        kinIds =  _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.raw"),
    params:
        prefix = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors"),
        prefix_filtered = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered"),
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

rule souporcell_check_plots:
    input:
        kin0 = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("souporcell_check-{bamtype}/souporcell_clusters_and_donors_filtered.raw"),
    output:
        _results("souporcell_check-{bamtype}/corrplot_r_by_donor.png"),
    params:
        outdir = _results("souporcell_check-{bamtype}"),
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
