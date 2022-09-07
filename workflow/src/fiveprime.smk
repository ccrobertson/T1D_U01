#!/usr/bin/env python3

from os.path import join
import os
import pandas as pd
from functools import partial
import math


name = config["name"]
_results = partial(os.path.join, "results/fiveprime", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


def iterate_samples(demuxlet_report):
    d = pd.read_csv(demuxlet_report)
    samples = d["Sample_ID"].tolist()
    return samples

def calculate_max_missing_by_batch(samplename, genotyped_donors):
    n_subjects = len(config["batches"][samplename].keys())
    n_donors = pd.read_table(genotyped_donors, header=None).shape[0]
    max_missing_rate = math.ceil((n_subjects-1)/(n_donors+n_subjects)*100, )/100
    return max_missing_rate

samples = iterate_samples(config["demuxlet_report"])
configfile: "results/multiome/Sample_islet_list_for_multiomics_Batches_long_format_fixedHPAP093_with_libname.json"


rule all:
    input:
        _results("library_config.json"),
        #expand(_results("tracks_by_sample/gex.{sample}.TPM.bw"), sample=samples)
        expand(_results("tracks/{sample}.TPM.{strand}.bw"),strand = ['forward', 'reverse'], sample="Sample_5125-NM-1-5GEX"),
        #expand(_results("ngsplot/{sample}/ngsplot_multiplot.png"), sample="Sample_5125-NM-1-5GEX"),
        expand(_results("deeptools/{sample}/profile_plot_tss_vs_intergenic.png"), sample="Sample_5125-NM-1-5GEX"),
        _results("demultiplex/corrplot_r_by_donor.png"),




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


# rule cell_ranger:
#     input:
#     output:
#     conda:
#         "cellranger"
#     shell:
#         """
#         cellranger
#         """


rule souporcell_fix_vcf:
    input:
        cluster_genotypes = _results("souporcell/cluster_genotypes.vcf"),
    output:
        cluster_genotypes_reformatted = _results("souporcell/cluster_genotypes_reformatted.vcf.gz"),
    params:
        sample = "NM1",
    conda:
        "genetics"
    shell:
        """
        awk -v OFS='\t' -v sample={params.sample} '$1!~/^#CHROM/ {{print $0}} $1~/^#CHROM/ {{for(i=10; i<=NF; ++i) $i=sample"_souporcell_"$i; print $0 }}' {input.cluster_genotypes} | grep -v BACKGROUND | bgzip > {output.cluster_genotypes_reformatted}
        tabix -p vcf {output.cluster_genotypes_reformatted}
        """

rule demultiplex:
    input:
        souporcell_vcf = _results("souporcell/cluster_genotypes_reformatted.vcf.gz"),
        imputed_vcf = "results/imputation/2022_02_16_T1D_genotypes/imputation_results/chrALL.donors_only.maf_gt_0.01__rsq_gt_0.95.dose.vcf.gz",
    output:
        kin0 =  _results("demultiplex/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("demultiplex/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("demultiplex/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("demultiplex/souporcell_clusters_and_donors_filtered.raw"),
    params:
        prefix = _results("demultiplex/souporcell_clusters_and_donors"),
        prefix_filtered = _results("demultiplex/souporcell_clusters_and_donors_filtered"),
        max_missing = calculate_max_missing_by_batch("Sample_5124-NM-1-hg38", config["genotyped_donors"])
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


rule demultiplex_plots:
    input:
        kin0 = _results("demultiplex/souporcell_clusters_and_donors_filtered.kin0"),
        kinMat = _results("demultiplex/souporcell_clusters_and_donors_filtered.king"),
        kinIds = _results("demultiplex/souporcell_clusters_and_donors_filtered.king.id"),
        genoraw = _results("demultiplex/souporcell_clusters_and_donors_filtered.raw"),
    output:
        _results("demultiplex/corrplot_r_by_donor.png"),
    params:
        outdir = _results("demultiplex"),
        donors = "HPAP105,HPAP093,ICRH139,ICRH142,HPAP107",
    conda:
        "Renv"
    shell:
        """
        Rscript workflow/scripts/demultiplexing_plots.R \
            --kin0 {input.kin0} \
            --kinMat {input.kinMat} \
            --kinIds {input.kinIds} \
            --donors "HPAP105,HPAP093,ICRH139,ICRH142,HPAP107" \
            --genoraw {input.genoraw} \
            --outdir {params.outdir}
        """


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregate signal plots
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




# rule scafe:
#     input:
#     output:
#     container:
#         "workflow/envs/scafe_latest.sif"
#     shell:
#         """
#         singularity shell --writable  workflow/envs/scafe_latest.sif
#         PATH=/SCAFE/scripts:$PATH
#         LC_ALL=C
#
#         SINGULARITYENV_APPEND_PATH=/SCAFE/scripts singularity exec --writable --cleanenv workflow/envs/scafe_latest.sif scafe.workflow.sc.solo --help
#
#         #singularity exec --writable --env APPEND_PATH=/SCAFE/scripts --env LC_ALL=C  scafe.workflow.sc.solo --help
#         """

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
# rule starsolo:
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
