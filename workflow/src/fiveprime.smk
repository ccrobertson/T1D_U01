#!/usr/bin/env python3

from os.path import join
import os
import pandas as pd
from functools import partial



#GOAL is to create a function to extract sample names from demuxlet file
#so that all we have to do is provide a config with "name","demuxlet_report",
#and "fastq_dir" to run the pipeline on a new set of 5GEX samples
def iterate_samples(demuxlet_report):
    d = pd.read_csv(demuxlet_report)
    samples = d["Sample_ID"].tolist()
    return samples

name = config["name"]
samples = iterate_samples(config["demuxlet_report"])

_results = partial(os.path.join, "results/fiveprime", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


rule all:
    input:
        _results("library_config.json"),
        expand(_results("tracks_by_sample/gex.{sample}.TPM.bw"), sample=samples)


rule symlinks:
    input:
        demuxlet_report = config["demuxlet_report"],
    output:
        json = _results("library_config.json"),
    params:
        indir = config["fastq_dir"],
        #want complete path for symlinks?
        outdir = os.path.join("/lab/work/ccrober/T1D_U01", _results("fastq")),
        samples = iterate_samples(config["demuxlet_report"]),
    shell:
        """
        createSymLinks () {{
          local target_dir=$1
          local link_dir=$2
          local sample=$3
          ln -s --force $target_dir/$sample/*R1_001.fastq.gz $link_dir/${{sample}}_R1_001.fastq.gz
          ln -s --force $target_dir/$sample/*R2_001.fastq.gz $link_dir/${{sample}}_R2_001.fastq.gz
        }}

        mkdir -p {params.outdir}
        for sample in {params.samples}
        do
          echo $sample
          createSymLinks {params.indir} {params.outdir} $sample
        done

        python workflow/scripts/build_fiveprime_json.py --demuxlet_report {input.demuxlet_report} --fastq_dir {params.outdir} --modality GEX > {output.json}
        """

### NOTE: Have not tried executing this step within the snakemake pipeline
# There might be issues with relative paths, and it may be necessary to provide
# complete paths to get it to work properly. Also not sure how the singularity
# pull from Sylabs within the nf pipeline will work. Might be ok!
rule nf_rna:
    input:
        json = _results("nf_gex_config.json"),
    output:
        _results("nf_gex_results/qc/{sample}-hg38.qc.txt"),
    params:
        outdir = _results("nf_gex_results")
    shell:
        """
        module load golang/1.13.5 singularity/3.5.1
        cd {params.outdir}
        nextflow run -resume -params-file {input.json} --chemistry V2 --results {params.outdir} pipelines/snRNAseq-NextFlow/main.nf
        """

rule cell_ranger:
    input:
    output:
    shell:
        """
        cellranger-arc

        """

rule bam_to_bigwig:
    input:
        bam = _results("nf_gex_results/prune/{sample}-hg38.before-dedup.bam"),
        blacklist = config['blacklist'],
    output:
        bw =  _results("tracks_by_sample/gex.{sample}.TPM.bw"),
    conda:
        "general"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} --normalizeUsing BPM --numberOfProcessors 20 --blackListFileName {input.blacklist}
        """


# rule qc_plots_by_sample:
#     input:
#         txt = _results("nf_gex_results/qc/{sample}-hg38.qc.txt"),
#     output:
#          _results("fiveprime_qc/{sample}/barcodes_empty.txt"),
#          _results("fiveprime_qc/{sample}/barcodes_nuclei.txt"),
#          _results("fiveprime_qc/{sample}/qc_grid_unfiltered.png"),
#          _results("fiveprime_qc/{sample}/qc_grid_prefiltered.png"),
#          _results("fiveprime_qc/{sample}/qc_grid_postfiltering.png"),
#     params:
#         outdir = _results("fiveprime_qc/{sample}"),
#         sample = "{sample}",
#     shell:
#         """
#         Rscript workflow/scripts/fiveprime_qc.R --input_txt {input.txt} --outdir {params.outdir} --sample {params.sample}
#         """
#
#
# rule prepare_rna_counts:
#     input:
#         _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx"),
#         _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv"),
#         barcodes_nuclei = _results("fiveprime_qc/{sample}/barcodes_nuclei.txt"),
#         barcodes_empty = _results("fiveprime_qc/{sample}/barcodes_empty.txt"),
#     output:
#         _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw/genes.tsv"),
#         counts_nuclei = _results("counts_by_sample_gex/{sample}/counts_nuclei.rds"),
#         counts_empty = _results("counts_by_sample_gex/{sample}/counts_empty.rds"),
#     params:
#         input_10x_dir = _results("nf_gex_results/starsolo/{sample}/{sample}.Solo.out/GeneFull_ExonOverIntron/raw"),
#     conda:
#         "Renv"
#     shell:
#         """
#         ln -s --relative --force {params.input_10x_dir}/features.tsv {params.input_10x_dir}/genes.tsv
#         Rscript workflow/scripts/prepare_rna_counts.R \
#             --input_10x_dir {params.input_10x_dir} \
#             --barcodes_nuclei {input.barcodes_nuclei} \
#             --barcodes_empty {input.barcodes_empty} \
#             --counts_nuclei {output.counts_nuclei} \
#             --counts_empty {output.counts_empty}
#         """


# #alignment
# rule alignment:
#     input:
#         insert_fastqs
#         barcode_fastqs
#     ouput:
#         bam
#         counts
#     params:
#         sample=
#         star_index=
#         gtf_file=
#         soloUMIlen=
#         barcode_whitelist=
#     shell:
#         """
#         STAR --soloBarcodeReadLength 0 \
#             --runThreadN 10 \
#             --outFileNamePrefix {params.sample}. \
#             --genomeLoad NoSharedMemory \
#             --runRNGseed 789727 \
#             --readFilesCommand gunzip -c \
#             --outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN \
#             --genomeDir {params.star_index} \
#             --outSAMtype BAM SortedByCoordinate \
#             --outSAMunmapped Within KeepPairs \
#             --sjdbGTFfile {params.gtf_file} \
#             --soloType Droplet \
#             --soloUMIlen {params.soloUMIlen} \
#             --soloFeatures Transcript3p Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS SJ Velocyto \
#             --soloMultiMappers Uniform PropUnique EM Rescue \
#             --soloUMIfiltering MultiGeneUMI \
#             --soloCBmatchWLtype 1MM_multi_pseudocounts \
#             --soloCellFilter None \
#             --soloCBwhitelist {params.barcode_whitelist} \
#             --readFilesIn {insert_fastq} {barcode_fastq}
#         """
#
#
#
#
# #Step 3 - SCAFE
# #rule scafe:
# #    input:
#
# ### What is the fragment length distribution?
# # cannot currently align R1
# # what proportion of R2 include TSO sequence -- meaning we sequenced through the insert?
# # oligo sequence is available here:  https://support.10xgenomics.com/permalink/6iTwTmSBVBPTYMXt3oVmD
# # The oligo sequences start on page 94.
