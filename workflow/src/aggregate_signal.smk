#!/usr/bin/env python3

from os.path import join
import os
from functools import partial
import pandas as pd
import math

#name = config["name"]

_results = partial(os.path.join, "results/multiome")
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")



rule all:
    input:
        #~~~~~~~ get tracks
        #expand(_results("tracks_by_sample/gex.{sample}.TPM.bw"), sample=samples),
        #expand(_results("tracks_by_sample/gex.{sample}.TPM.{strand}.bw"), sample=samples, strand=['forward','reverse']),
        #~~~~~ rerun macs2 to get summits
        #expand(_results("macs2/{sample}_summits.bed"), sample="Sample_5124-NM-1-hg38"),
        #expand(_results("ngsplot/{sample}/ngsplot_multiplot.png"), sample="Sample_5124-NM-1-hg38"),
        expand(_results("deeptools/{sample}/profile_plot_tss_vs_intergenic.png"), sample="Sample_5124-NM-1-hg38"),


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregate signal plots -- MULTIOME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### DEFINING INTERGENIC ATAC SUMMITS
rule bamtobed:
    input:
        bam = "results/multiome/nf_atac_results/prune/{sample}.pruned.bam"
    output:
        bed = _results("macs2/{sample}.pruned.bed"),
    container:
        "workflow/envs/porchard_snatac_general_20220107.sif"
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output.bed}
        """

rule run_macs2:
    input:
        bed = _results("macs2/{sample}.pruned.bed"),
    output:
        bdg = _results("macs2/{sample}_treat_pileup.bdg"),
        peaks = _results("macs2/{sample}_peaks.narrowPeak"),
        summits = _results("macs2/{sample}_summits.bed"),
    params:
        outdir = _results("macs2"),
        sample = "{sample}",
    container:
        "workflow/envs/porchard_snatac_general_20220107.sif"
    shell:
        """
        macs2 callpeak -t {input.bed} \
            -f BED \
            -g hs \
            --outdir {params.outdir} \
            -n {params.sample} \
            --bdg --SPMR \
            --nomodel \
            --shift -100 --extsize 200 \
            --seed 762873 \
            --call-summits \
            --keep-dup all
        """

rule define_features:
    input:
        summits = _results("macs2/Sample_5124-NM-1-hg38_summits.bed"),
        gtf = "resources/hg38/gencode.v39.annotation.gtf",
    output:
        _results("features/tss_nonoverlapping_genes.bed"),
        _results("features/atac_summits_intergenic.bed"),
        #_results("features/atac_summits_1kb_intergenic.bed"),
        #_results("features/atac_summits_1kb_all.bed"),
    params:
        outdir = _results("features"),
    shell:
        """
        Rscript workflow/scripts/define_features.R \
            --summits {input.summits} \
            --gtf {input.gtf} \
            --outdir {params.outdir}
        """


### PLOTS WITH BAM
rule ngs_plots:
    input:
        bam =  _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex.bam"),
        tss = _results("features/tss_nonoverlapping_genes.bed"),
        intergenic = _results("features/atac_summits_intergenic.bed"),
    output:
        config = _results("ngsplot/{sample}/ngsplot_config.txt"),
        png = _results("ngsplot/{sample}/ngsplot_multiplot.png"),
    params:
        outdir = _results("ngsplot/{sample}"),
    shell:
        """
        ## make config
        echo -e "{input.bam}\t{input.tss}\tTSS" > {output.config}
        echo -e "{input.bam}\t{input.intergenic}\tIntergenic_Summit" >> {output.config}

        ## run ngs
        export NGSPLOT=/lab/work/ccrober/sw/ngsplot
        ngs.plot.r -G hg38 -R bed -C {output.config} -O {params.outdir} -L 1000
        """


### PLOTS WITH STRANDED BIGWIGS
rule bam_to_bigwig_stranded:
    input:
        bam = _results("bam_pass_qc_barcodes/{sample}/pass_qc_barcodes_gex.bam"),
        #bam = _results("nf_gex_results/prune/{sample}.before-dedup.bam"),
        blacklist = config["blacklist"],
    output:
        bw =  _results("tracks_by_sample/gex.{sample}.TPM.{strand}.bw"),
    params:
        strand = "{strand}"
    conda:
        "deeptools"
    shell:
        """
        bamCoverage -b {input.bam} \
            -o {output.bw} \
            --normalizeUsing BPM \
            --filterRNAstrand {params.strand} \
            --numberOfProcessors 8 \
            --blackListFileName {input.blacklist}
        """


rule deeptools_matrix:
    input:
        bws = expand(_results("tracks_by_sample/gex.{{sample}}.TPM.{strand}.bw"), strand=['forward','reverse']),
        tss = _results("features/tss_nonoverlapping_genes.bed"),
        intergenic = _results("features/atac_summits_intergenic.bed"),
    output:
        matrix = _results("deeptools/{sample}/matrix.mat.gz"),
    conda:
        "deeptools",
    shell:
        """
        computeMatrix reference-point -S {input.bws} \
            -R {input.tss} {input.intergenic}\
            --referencePoint center \
            --upstream 1000 \
            --downstream 1000 \
            -o {output.matrix} \
            --binSize 50 \
            --numberOfProcessors 20
        """


rule deeptools_plot:
    input:
        matrix = _results("deeptools/{sample}/matrix.mat.gz"),
    output:
        png = _results("deeptools/{sample}/profile_plot_tss_vs_intergenic.png"),
    conda:
        "deeptools",
    shell:
        """
        plotProfile -m {input.matrix} \
            -out {output.png} \
            --numPlotsPerRow 1 \
            --perGroup \
            --plotTitle "Aggregate signal around TSS or intergenic ATAC-seq summits"
        """



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregate signal plots -- 5GEX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### USING BAM
rule ngs_plots:
    input:
        bam =  _results("cellranger/outs/possorted_genome_bam.bam"),
        tss = "results/multiome/features/tss_nonoverlapping_genes.bed",
        intergenic = "results/multiome/features/atac_summits_intergenic.bed",
    output:
        config = _results("ngsplot/{sample}/ngsplot_config.txt"),
        pdf = _results("ngsplot/{sample}.avgprof.pdf"),
    params:
        outdir = _results("ngsplot/{sample}"),
    shell:
        """
        ## make config
        echo -e "{input.bam}\t{input.tss}\tTSS" > {output.config}
        echo -e "{input.bam}\t{input.intergenic}\tIntergenic_Summit" >> {output.config}

        ## run ngs
        export NGSPLOT=/lab/work/ccrober/sw/ngsplot
        ngs.plot.r -G hg38 -R bed -C {output.config} -O {params.outdir} -L 1000
        """

### USING STRANDED BIGWIGS
rule bam_to_bigwig:
    input:
        bam = _results("cellranger/outs/possorted_genome_bam.bam"),
        blacklist = config['blacklist'],
    output:
        bw =  _results("tracks/{sample}.TPM.{strand}.bw"),
    params:
        strand = "{strand}",
    conda:
        "deeptools"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} \
            --normalizeUsing BPM \
            --filterRNAstrand {params.strand} \
            --numberOfProcessors 20 \
            --blackListFileName {input.blacklist}
        """

rule deeptools_matrix:
    input:
        bws = expand(_results("scafe/ctss_to_bigwig/{{sample}}/wig/{{sample}}.cpm.{strand}.bw"), strand=['fwd','rev']),
        tss = "results/multiome/features/tss_nonoverlapping_genes.bed",
        intergenic = "results/multiome/features/atac_summits_intergenic.bed",
    output:
        matrix = _results("deeptools/{sample}/matrix.mat.gz"),
    conda:
        "deeptools",
    shell:
        """
        computeMatrix reference-point -S {input.bws} \
            -R {input.tss} {input.intergenic}\
            --referencePoint center \
            --upstream 1000 \
            --downstream 1000 \
            -o {output.matrix} \
            --missingDataAsZero \
            --numberOfProcessors 20
        """


rule deeptools_plot:
    input:
        matrix = _results("deeptools/{sample}/matrix.mat.gz"),
    output:
        png = _results("deeptools/{sample}/profile_plot_tss_vs_intergenic.png"),
    conda:
        "deeptools",
    shell:
        """
        plotProfile -m {input.matrix} \
            -out {output.png} \
            --numPlotsPerRow 1 \
            --perGroup \
            --plotTitle "Aggregate signal around TSS or intergenic ATAC-seq summits"
        """
