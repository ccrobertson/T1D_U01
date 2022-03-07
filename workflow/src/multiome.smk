#!/usr/bin/env python3

#IONICE = 'ionice -c2 -n7'


#RESULTS = config["results"]
RESULTS = "results/freeze1"
OUTDIR = join(RESULTS, "{sample}")

rule filter_barcodes_rna:
    input:
        rna = join("data/rna/nucleus-qc", "{sample}-rna.txt"),
    output:
        rna = join(OUTDIR, "barcodes_rna.txt"),
    params:
        umi_min = config["umi_min"],
        pct_max_ribo = config["pct_max_ribo"],
        min_genes = config["ngenes_min"],
        emptyDrops = config["dropletUtils_rna"],
        mito_knees = config["mitochondrial_knees_rna"],
    shell:
        """
        Rscript scripts/select_barcodes_rna_dropletUtils.R -i {input.rna} \
            -o {output.rna} --sample {wildcards.sample} \
            --umi_min {params.umi_min}  --ngenes {params.min_genes} \
            --pct_max_ribo {params.pct_max_ribo} \
            --emptyDrops {params.emptyDrops} \
            --mito_knees={params.mito_knees}
        """


rule filter_barcodes_atac:
    input:
        atac = config["atac_qc_file"],
    output:
        atac = join(OUTDIR, "barcodes_atac.txt"),
    params:
        tss_min = config["tss_min"],
        tss_max = config["tss_max"],
        hqaa_min = config["hqaa_min"],
        pct_max = config["pct_max_hqaa"],
        autosome_max = config["autosome_max"],
        mito_knees = config["mitochondrial_knees"],
    shell:
        """
        Rscript scripts/select_barcodes_atac.R -i {input.atac} \
            -o {output.atac} --sample {wildcards.sample} \
            --tss_min {params.tss_min} --tss_max {params.tss_max} \
            --hqaa_min {params.hqaa_min} --pct_max {params.pct_max} \
            --autosome_max {params.autosome_max} \
            --mito_knees={params.mito_knees}
        """
