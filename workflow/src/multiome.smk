#!/usr/bin/env python3

#IONICE = 'ionice -c2 -n7'


#RESULTS = config["results"]
RESULTS = "results/freeze1"
OUTDIR = join(RESULTS, "multiome")


# rule all:
#     input:
#         expand(join(RIBO_QC_DIR, "{sample}.txt"), sample=samples),
#         expand(join(OUTDIR, "atac_processed.rds"), sample=samples),
#         expand(join(OUTDIR, "rna_processed.rds"), sample=samples),
#         expand(
#             join(FIRST_PASS_DIR,
#                  "liger__k{k}__l{l}__knn{knn_k}__res{res}.html"),
#             k=kvals, l=lambdas, knn_k=knn_ks, res=resolutions
#         ),


rule filter_barcodes:
    input:
        gex = join(RESULTS, "nf_atac_results/ataqv/single-nucleus/{sample}.txt"),
        atac = join(RESULTS,"nf_gex_results/qc/{sample}.qc.txt"),
    output:
        barcodes = join(OUTDIR, "barcode_filtering/barcodes_{sample}.txt"),
    params:
        umi_min = config["umi_min"],
        pct_max_ribo = config["pct_max_ribo"],
        min_genes = config["ngenes_min"],
        emptyDrops = config["dropletUtils_rna"],
        mito_knees = config["mitochondrial_knees_rna"],
        tss_min = config["tss_min"],
        tss_max = config["tss_max"],
        hqaa_min = config["hqaa_min"],
        pct_max = config["pct_max_hqaa"],
        autosome_max = config["autosome_max"],
        mito_knees = config["mitochondrial_knees"],
    shell:
        """
        Rscript scripts/cross_modality_qc.R -ig {input.gex} -ia {input.atac} \
            -o {output.barcodes} --sample {wildcards.sample} \
            --umi_min {params.umi_min}  --ngenes {params.min_genes} \
            --pct_max_ribo {params.pct_max_ribo} \
            --emptyDrops {params.emptyDrops} \
            --mito_knees={params.mito_knees} \
            --tss_min {params.tss_min} --tss_max {params.tss_max} \
            --hqaa_min {params.hqaa_min} --pct_max {params.pct_max} \
            --autosome_max {params.autosome_max}
        """
