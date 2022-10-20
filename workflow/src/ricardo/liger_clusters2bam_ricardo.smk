#!/usr/bin/env python3
from os.path import join
from glob import glob


def get_bams_from_condition(condition, clus):
    samples_dic = {
        "healthy": samples_healthy,
        "healthymock": samples_healthymock,
        "pret1d": samples_pret1d,
        "t1d": samples_t1d,
        "aab": samples_aab,
        "mock": samples_mock,
        "cytokines": samples_cytokines,
        "cvb4": samples_cvb4,
        
    }
    outfiles = expand(
        join(OUTBAM_DIR, "{sample}__atac__{clus}.bam"),
        sample=samples_dic[condition], clus=clus
    )
    return outfiles


def list_other_clusters(cluster):
    return [i for i in clusters if i != cluster]


# Config
unique_id = config["liger_obj"]
# CLUSTERS_N = config["n_clusters"]
MIN_N = config["min_n"]
IONICE = "ionice -c2 -n7"
EXT = config["summit_ext"]


# Directories
RESULTS = join(config["results"], "clusters", unique_id)
# LIGEROBJ_DIR = config["liger_dir"]
ATACBAM_DIR = config["atac_bam_dir"]
BARCODES_DIR = join(RESULTS, "barcodes")
OUTBAM_DIR = join(RESULTS, "bam_files")
OUTBAM_COND_DIR = join(RESULTS, "bam_files_conditions")
MACS2_DIR = join(RESULTS, "macs2")
MACS2_IND_DIR = join(RESULTS, "macs2_per_sample")
MACS2_COND_DIR = join(RESULTS, "macs2_per_condition")
COUNTS_DIR = join(RESULTS, "counts")
VSIGNAL_DIR = join(RESULTS, "vplots")
ATAQV_DIR = join(RESULTS, "ataqv")
VIEWER_DIR = join(RESULTS, "viewer")
BROWSER_DIR = join(RESULTS, "genome_browser")

# Wildcards
samples = config["samples"]
samples_healthy = config["samples_healthy"]
samples_healthymock = config["samples_healthymock"]
samples_pret1d = config["samples_pret1d"]
samples_t1d = config["samples_t1d"]
samples_aab = config["samples_aab"]
samples_mock = config["samples_mock"]
samples_cytokines = config["samples_cytokines"]
samples_cvb4 = config["samples_cvb4"]

beta_cluster = config["beta_cluster"]

# clusters = [x for x in range(0, CLUSTERS_N)]
# clusters = list(
#     set(
#         [
#         re.sub(r".*_atac_([A-Z0-9]+)\.txt", r"\1", i) 
#         for i in glob(join(BARCODES_DIR, "*.txt"))
#         ]
#     )
# )

clusters = config["clusters"]

conditions = ["healthymock", "aab", "cytokines", "cvb4"]

# ligerex = join(LIGEROBJ_DIR, unique_id) + ".rds"



rule all:
    input:
        expand(
            join(OUTBAM_COND_DIR, "{cond}_atac__{clus}.bam.bai"),
            cond=conditions, clus=clusters
        ),
        # expand(
        #     join(OUTBAM_DIR, "{sample}__atac__{cluster}.bam"),
        #     sample=samples, cluster=clusters
        # ),
        # expand(
        #     join(OUTBAM_DIR, "{sample}__atac__{cluster}.bam.bai"),
        #     sample=samples, cluster=clusters
        # ),
        # # expand(
        # #     join(
        # #         MACS2_IND_DIR, 
        # #         "{sample}__atac__{cluster}_treat_pileup.normalized.bw"
        # #     ),
        # #     sample=samples, cluster=clusters
        # # ),
        # expand(
        #     join(
        #         MACS2_IND_DIR, 
        #         "{sample}__atac__{cluster}_peaks.broadPeak.noblacklist"
        #     ),
        #     sample=samples, cluster=clusters
        # ),
        # expand(
        #     join(
        #         MACS2_IND_DIR, 
        #         "{sample}__atac__{cluster}_summits.noblacklist.ext{ext}.bed"
        #     ),
        #     sample=samples, cluster=clusters, ext=EXT
        # ),
        # expand(
        #     join(
        #         MACS2_COND_DIR, 
        #         "{cond}_atac__{clus}_treat_pileup.normalized.bw"
        #     ),
        #     cond=conditions, clus=beta_cluster
        # ),
        # expand(
        #     join(
        #         MACS2_COND_DIR, 
        #         "{cond}_atac__{clus}_peaks.broadPeak.noblacklist"
        #     ),
        #     cond=conditions, clus=clusters
        # ),
        # expand(
        #     join(
        #         MACS2_COND_DIR, 
        #         "{cond}_atac__{clus}_summits.noblacklist.ext{ext}.bed"
        #     ),
        #     cond=conditions, clus=beta_cluster, ext=EXT
        # ),
        # expand(join(OUTBAM_DIR, "atac__{cluster}.bam"), cluster=clusters),
        # expand(
        #     join(MACS2_DIR, "{cluster}_treat_pileup.normalized.bw"), 
        #     cluster=clusters
        # ),
        # expand(
        #     join(MACS2_DIR, "{cluster}_summits.noblacklist.ext{ext}.bed"),
        #     cluster=clusters, ext=EXT
        # ),
        # # expand(join(ATAQV_DIR, "{cluster}.ataqv.json.gz"), cluster=clusters),
        # expand(join(VSIGNAL_DIR, "{cluster}.png"), cluster=clusters),
        # join(VIEWER_DIR, "index.html"),
        # expand(
        #     join(MACS2_DIR, "{cluster}_peaks" + ".min{}.bed".format(MIN_N)),
        #     cluster=clusters
        # ),
        # expand(join(MACS2_DIR, "{cluster}_summits.ext{ext}.min{MIN_N}.bed"),
        #        cluster=clusters, MIN_N=MIN_N, ext=EXT),
        # expand(join(MACS2_DIR, "{cluster}_peaks.min{MIN_N}.bed"),
        #        cluster=clusters, MIN_N=MIN_N),
        # expand(join(MACS2_DIR, "{cluster}_summits.ext{ext}.cluster_specific.bed"),
        #        cluster=clusters, ext=EXT),
        # expand(join(MACS2_DIR, "{cluster}_peaks.cluster_specific.bed"),
        #        cluster=clusters),
        # expand(join(COUNTS_DIR, "{cluster}_peaks.min{MIN_N}.txt"),
        #        cluster=clusters, MIN_N=MIN_N),
        # expand(join(COUNTS_DIR, "{cluster}_summits.ext{ext}.min{MIN_N}.txt"),
        #        cluster=clusters, MIN_N=MIN_N, ext=EXT),
        # join(BROWSER_DIR, "trackDb.peaks.txt"),
        # join(BROWSER_DIR, "trackDb.summits150.txt"),
        # join(BROWSER_DIR, "trackDb.bigwig.txt"),
        # join(BROWSER_DIR, "trackDb.bigwig_conditions.txt"),



rule fetch_atac_barcodes:
    """Fake rule - just touch files with 0 barcodes. Use Rmd file to run"""
    output:
        expand(
            join(BARCODES_DIR, "barcodes_{sample}_atac_{cluster}.txt"),
            sample=samples, cluster=clusters
        ),
    params:
        OUTBAM_DIR
    shell:
        """
        # Rscript get_cluster_barcodes.R {input} {params}
        touch {output}
        """

rule get_bam_from_barcodes:
    input:
        bam = join(ATACBAM_DIR, "{sample}-hg19.pruned.bam"),
        barcodes = join(BARCODES_DIR, "barcodes_{sample}_atac_{cluster}.txt"),
    output:
        join(OUTBAM_DIR, "{sample}__atac__{cluster}.bam"),
    shell:
        """
        {IONICE} python scripts/filter-bam-by-barcode.py {input.bam} \
            {output} {input.barcodes}
        """

rule bamtobed_individual:
    input:
        rules.get_bam_from_barcodes.output,
    output:
        temp(join(OUTBAM_DIR, "{sample}__atac__{cluster}.bam"))
    shell:
        """
        {IONICE} bedtools bamtobed -i {input} > {output}
        """

rule index_bams:
    input:
        rules.get_bam_from_barcodes.output,
    output:
        join(OUTBAM_DIR, "{sample}__atac__{cluster}.bam.bai"),
    shell:
        """
        {IONICE} samtools index {input}
        """

rule merge_bams:
    input:
        lambda wildcards: expand(
            join(OUTBAM_DIR, "{sample}__atac__{cluster}.bam"),
            sample=samples, cluster=wildcards.cluster
        ),
    output:
        join(OUTBAM_DIR, "atac__{cluster}.bam")
    shell:
        """
        {IONICE} samtools merge {output}.tmp {input} &&
        {IONICE} samtools sort -o {output} -T {output}xx {output}.tmp &&
        rm {output}.tmp
        """

rule index_bams2:
    input:
        rules.merge_bams.output,
    output:
        join(OUTBAM_DIR, "atac__{cluster}.bam.bai"),
    shell:
        """
        {IONICE} samtools index {input}
        """

rule bamtobed_cluster:
    input:
        rules.merge_bams.output,
    output:
        temp(join(OUTBAM_DIR, "atac__{cluster}.bed"))
    shell:
        """
        {IONICE} bedtools bamtobed -i {input} > {output}
        """

# Process BAMs per sample

rule call_broad_peaks_individual:
    input:
        rules.bamtobed_individual.output,
    output:
        bed = join(MACS2_IND_DIR,
                   "{sample}__atac__{cluster}_peaks.broadPeak"),
        bdg = join(MACS2_IND_DIR,
                   "{sample}__atac__{cluster}_treat_pileup.bdg"),
    params:
        name = "{sample}__atac__{cluster}",
        stderr_location = join(MACS2_IND_DIR, 
                               "{sample}__atac__{cluster}.macs2.out"),
        outdir = MACS2_IND_DIR,
        
    shell:
        """
        {IONICE} macs2 callpeak -t {input} --outdir {params.outdir} \
            -n {params.name} -g hs \
            --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad \
            --keep-dup all &> {params.stderr_location}
        """

rule call_summits_individual:
    input:
        rules.bamtobed_individual.output,
    output:
        join(MACS2_IND_DIR, "{sample}__atac__{cluster}_summits.bed"),
    params:
        name = "{sample}__atac__{cluster}",
        stderr_location = join(MACS2_IND_DIR, 
                               "{sample}__atac__{cluster}.narrow.macs2.out"),
        outdir = MACS2_IND_DIR,
    shell:
        """
        {IONICE} macs2 callpeak -t {input} --outdir {params.outdir} \
            -n {params.name} -g hs \
            --nomodel --shift -100 --seed 762873 --extsize 200 --call-summits \
            --keep-dup all &> {params.stderr_location}
        """

rule blacklist_filter_individual:
    input:
        bed = rules.call_broad_peaks_individual.output.bed,
        bl1 = config["blacklists"][0],
        bl2 = config["blacklists"][1],
    output:
        join(MACS2_IND_DIR,
             "{sample}__atac__{cluster}_peaks.broadPeak.noblacklist")
    shell:
        """
        {IONICE} intersectBed -v -a {input.bed} -b {input.bl1} | \
            intersectBed -v -a stdin -b {input.bl2} > {output}
        """

rule extend_summits_individual:
    input:
        bed = rules.call_summits_individual.output,
        bl1 = config["blacklists"][0],
        bl2 = config["blacklists"][1],
    params:
        genome = "/lab/data/reference/human/hg19/hg19.chrom_sizes"
    output:
        join(MACS2_IND_DIR,
             "{sample}__atac__{cluster}_summits.noblacklist.ext{ext}.bed")
    shell:
        """
        bedtools slop -i {input.bed} -g {params} -b {EXT} | sortBed | \
            mergeBed -c 4 -o distinct | \
            intersectBed -v -a stdin -b {input.bl1} | \
            intersectBed -v -a stdin -b {input.bl2} > {output}
        """

rule normalize_bedgraph_individual:
    input:
        rules.call_broad_peaks_individual.output.bdg,
    output:
        join(MACS2_IND_DIR, 
            "{sample}__atac__{cluster}_treat_pileup.normalized.bdg")
    params:
        out = join(MACS2_IND_DIR, "{sample}__atac__{cluster}.macs2.out")
    shell:
        """
        {IONICE} cat {input} | \
            awk -v NTAGS=$(grep 'total tags in treatment' {params.out} | \
            awk '{{print $NF}}') '{{$4=$4*(10000000/NTAGS); print}}' | \
            perl -pe 's/\\s+/\\t/g; s/$/\\n/' | grep -v "_" | \
            LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule bedgraph_to_bigwig_individual:
    input:
        rules.normalize_bedgraph_individual.output
    output:
        join(MACS2_IND_DIR, 
            "{sample}__atac__{cluster}_treat_pileup.normalized.bw")
    params:
        config["genome_size"]
    shell:
        """
        {IONICE} bedtools slop -i {input} -g {params} -b 0 | \
            awk '$2 < $3 {{print}}' > {output}.tmp &&\
        {IONICE} bedGraphToBigWig {output}.tmp {params} {output} &&\
        rm {output}.tmp
        """


# Process BAMs per condition

rule merge_bams_condition:
    input:
        lambda wildcards: get_bams_from_condition(
            wildcards.cond, wildcards.clus
        ),
    output:
        join(OUTBAM_COND_DIR, "{cond}_atac__{clus}.bam"),
    shell:
        """
        {IONICE} samtools merge -f {output}.tmp {input} &&
        {IONICE} samtools sort -o {output} \
            -T {output}xx {output}.tmp &&
        rm {output}.tmp
        """

rule index_bams_condition:
    input:
        rules.merge_bams_condition.output,
    output:
        join(OUTBAM_COND_DIR, "{cond}_atac__{clus}.bam.bai"),
    shell:
        """
        {IONICE} samtools index {input}
        """

rule bamtobed_condition:
    input:
        rules.merge_bams_condition.output,
    output:
        temp(join(OUTBAM_DIR, "{cond}__atac__{clus}.bed"))
    shell:
        """
        {IONICE} bedtools bamtobed -i {input} > {output}
        """

rule call_broad_peaks_condition:
    input:
        rules.bamtobed_condition.output,
    output:
        bed = join(MACS2_COND_DIR,
                   "{cond}_atac__{clus}_peaks.broadPeak"),
        bdg = join(MACS2_COND_DIR,
                   "{cond}_atac__{clus}_treat_pileup.bdg"),
    params:
        name = "{cond}_atac__{clus}",
        stderr_location = join(MACS2_COND_DIR, 
                               "{cond}_atac__{clus}.macs2.out"),
        outdir = MACS2_COND_DIR,
    shell:
        """
        {IONICE} macs2 callpeak -t {input} --outdir {params.outdir} \
            -n {params.name} -g hs \
            --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad \
            --keep-dup all &> {params.stderr_location}
        """

rule call_summits_condition:
    input:
        rules.bamtobed_condition.output,
    output:
        join(MACS2_COND_DIR, "{cond}_atac__{clus}_summits.bed"),
    params:
        name = "{cond}_atac__{clus}",
        stderr_location = join(MACS2_COND_DIR, 
                               "{cond}_atac__{clus}.narrow.macs2.out"),
        outdir = MACS2_COND_DIR,
    shell:
        """
        {IONICE} macs2 callpeak -t {input} --outdir {params.outdir} \
            -n {params.name} -g hs \
            --nomodel --shift -100 --seed 762873 --extsize 200 --call-summits \
            --keep-dup all &> {params.stderr_location}
        """

rule blacklist_filter_condition:
    input:
        bed = rules.call_broad_peaks_condition.output.bed,
        bl1 = config["blacklists"][0],
        bl2 = config["blacklists"][1],
    output:
        join(MACS2_COND_DIR,
             "{cond}_atac__{clus}_peaks.broadPeak.noblacklist")
    shell:
        """
        {IONICE} intersectBed -v -a {input.bed} -b {input.bl1} | \
            intersectBed -v -a stdin -b {input.bl2} > {output}
        """

rule extend_summits_conditions:
    input:
        bed = rules.call_summits_condition.output,
        bl1 = config["blacklists"][0],
        bl2 = config["blacklists"][1],
    params:
        genome = "/lab/data/reference/human/hg19/hg19.chrom_sizes"
    output:
        join(MACS2_COND_DIR,
             "{cond}_atac__{clus}_summits.noblacklist.ext{ext}.bed")
    shell:
        """
        bedtools slop -i {input.bed} -g {params} -b {EXT} | sortBed | \
            mergeBed -c 4 -o distinct | \
            intersectBed -v -a stdin -b {input.bl1} | \
            intersectBed -v -a stdin -b {input.bl2} > {output}
        """

rule normalize_bedgraph_condition:
    input:
        rules.call_broad_peaks_condition.output.bdg,
    output:
        join(MACS2_COND_DIR, 
            "{cond}_atac__{clus}_treat_pileup.normalized.bdg")
    params:
        out = join(MACS2_COND_DIR, "{cond}_atac__{clus}.macs2.out")
    shell:
        """
        {IONICE} cat {input} | \
            awk -v NTAGS=$(grep 'total tags in treatment' {params.out} | \
            awk '{{print $NF}}') '{{$4=$4*(10000000/NTAGS); print}}' | \
            perl -pe 's/\\s+/\\t/g; s/$/\\n/' | grep -v "_" | \
            LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule bedgraph_to_bigwig_condition:
    input:
        rules.normalize_bedgraph_condition.output
    output:
        join(MACS2_COND_DIR, 
            "{cond}_atac__{clus}_treat_pileup.normalized.bw")
    params:
        config["genome_size"]
    shell:
        """
        {IONICE} bedtools slop -i {input} -g {params} -b 0 | \
            awk '$2 < $3 {{print}}' > {output}.tmp &&\
        {IONICE} bedGraphToBigWig {output}.tmp {params} {output} &&\
        rm {output}.tmp
        """


# Process pseudo-bulk per cluster

rule call_broad_peaks:
    input:
        rules.bamtobed_cluster.output,
    output:
        bed = join(MACS2_DIR, "{cluster}_peaks.broadPeak"),
        bdg = join(MACS2_DIR, "{cluster}_treat_pileup.bdg"),
    params:
        name = "{cluster}",
        stderr_location = join(MACS2_DIR, "{cluster}.macs2.out"),
        outdir = MACS2_DIR
    shell:
        """
        {IONICE} macs2 callpeak -t {input} --outdir {params.outdir} \
            -n {params.name} -g hs \
            --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad \
            --keep-dup all &> {params.stderr_location}
        """

rule call_summits:
    input:
        rules.bamtobed_cluster.output,
    output:
        join(MACS2_DIR, "{cluster}_summits.bed")
    params:
        name = "{cluster}",
        stderr_location = join(MACS2_DIR, "{cluster}.narrow.macs2.out"),
        outdir = MACS2_DIR
    shell:
        """
        {IONICE} macs2 callpeak -t {input} --outdir {params.outdir} \
            -n {params.name} -g hs \
            --nomodel --shift -100 --seed 762873 --extsize 200 --call-summits \
            --keep-dup all &> {params.stderr_location}
        """

rule blacklist_filter:
    input:
        bed = rules.call_broad_peaks.output.bed,
        bl1 = config["blacklists"][0],
        bl2 = config["blacklists"][1],
    output:
        join(MACS2_DIR, "{cluster}_peaks.broadPeak.noblacklist")
    shell:
        """
        {IONICE} intersectBed -v -a {input.bed} -b {input.bl1} | \
            intersectBed -v -a stdin -b {input.bl2} > {output}
        """

rule extend_summits:
    input:
        bed = rules.call_summits.output,
        bl1 = config["blacklists"][0],
        bl2 = config["blacklists"][1],
    params:
        genome = "/lab/data/reference/human/hg19/hg19.chrom_sizes"
    output:
        join(MACS2_DIR, "{cluster}_summits.noblacklist.ext{ext}.bed")
    shell:
        """
        bedtools slop -i {input.bed} -g {params} -b {EXT} | sortBed | \
            mergeBed -c 4 -o distinct | \
            intersectBed -v -a stdin -b {input.bl1} | \
            intersectBed -v -a stdin -b {input.bl2} > {output}
        """


rule normalize_bedgraph:
    input:
        rules.call_broad_peaks.output.bdg,
    output:
        join(MACS2_DIR, "{cluster}_treat_pileup.normalized.bdg")
    params:
        out = join(MACS2_DIR, "{cluster}.macs2.out")
    shell:
        """
        {IONICE} cat {input} | \
            awk -v NTAGS=$(grep 'total tags in treatment' {params.out} | \
            awk '{{print $NF}}') '{{$4=$4*(10000000/NTAGS); print}}' | \
            perl -pe 's/\\s+/\\t/g; s/$/\\n/' | grep -v "_" | \
            LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

rule bedgraph_to_bigwig:
    input:
        rules.normalize_bedgraph.output
    output:
        join(MACS2_DIR, "{cluster}_treat_pileup.normalized.bw")
    params:
        config["genome_size"]
    shell:
        """
        {IONICE} bedtools slop -i {input} -g {params} -b 0 | \
            awk '$2 < $3 {{print}}' > {output}.tmp &&\
        {IONICE} bedGraphToBigWig {output}.tmp {params} {output} &&\
        rm {output}.tmp
        """

rule ataqv:
    input:
        md_bam = rules.merge_bams.output,
        bai = rules.index_bams2.output,
        peaks = rules.blacklist_filter.output,
    output:
        metrics = join(ATAQV_DIR, "{cluster}.ataqv.json.gz"),
        stdout_destination = join(ATAQV_DIR, "{cluster}.ataqv.out")
    params:
        name = "{cluster}",
        description = "{cluster}", 
        organism = "human",
        tss_file = config["tss_file"],
        bl = " ".join(
            ["--excluded-region-file {}".format(x) 
             for x in config["blacklists"]]
        ),
    shell:
        """
        {IONICE} ataqv --peak-file {input.peaks} \
            --name {params.description} --metrics-file {output.metrics} \
            {params.bl} --tss-file {params.tss_file} \
            --ignore-read-groups {params.organism} {input.md_bam} \
            > {output.stdout_destination}
        """

rule atqv_session:
    input:
        expand(rules.ataqv.output.metrics, cluster=clusters),
    output:
        join(VIEWER_DIR, "index.html"),
    params:
        outdir = VIEWER_DIR,
    shell:
        """
        mkarv --force -d 'LIGER intergration' {params} {input}
        """

rule vsignal:
    input:
        bam = rules.merge_bams.output,
        bai = rules.index_bams2.output,

    output:
        vsignal = join(VSIGNAL_DIR, "{cluster}.vsignal.gz"),
        vplot = join(VSIGNAL_DIR, "{cluster}.png")
    threads:
        4
    params:
        handle = join(VSIGNAL_DIR, "{cluster}"),
        vplot_name = "{cluster}" + "_CTCF_cohesin",
        flags = "-r 500 -f 3 -F 4 -F 8 -q 30",
        motif = config["ctcf_cohesin"]
    resources:
        io_limit = 1
    shell:
        """
        {IONICE} measure_signal -p {threads} {params.flags} {input.bam} \
            {params.motif} | gzip -c > {output.vsignal}
        Rscript ~albanus/scripts/makeVplots.R -f {output.vsignal} \
            -o {params.handle} -n {params.vplot_name} --ylim 1.5 --split
        """

# rule generate_report:
#     input:
#         ligerex,
#     output:
#         join(RESULTS, "report.html"),
#     shell:
#         """
#         tmp="report_$(date '+%Y_%m_%d_%H_%M_%S_%N').html"
#         Rscript -e "rmarkdown::render('scripts/liger_object_qc.Rmd', \
#             params=list(infile='{input}'), output_file='${{tmp}}')" &&
#         mv scripts/${{tmp}} {output}
#         """

rule peak_reproducibility:
    input:
        rules.blacklist_filter.output,
    output:
        join(MACS2_DIR, "{cluster}_peaks.broadPeak.reproducible"),
    params:
        glob_pattern = join(MACS2_IND_DIR, "*__{cluster}_peaks*.noblacklist"),
    shell:
        """
        intersectBed -wb -wa  -a {input} -b {params} | cut -f 1-4,14 | \
            sed 's/__.*//g' | uniq > {output}
        """

rule make_master_peaks:
    input:
        rules.peak_reproducibility.output
    output:
        join(MACS2_DIR, "{cluster}_peaks" + ".min{}.bed".format(MIN_N))
    params:
        min_n = MIN_N
    shell:
        """
        Rscript scripts/subset_reproducible_peaks.R {input} {output} {params}
        """

rule cluster_specific_peaks:
    input:
        a = rules.make_master_peaks.output,
        b = lambda wildcards: expand(
            rules.make_master_peaks.output,
            cluster = list_other_clusters(wildcards.cluster)
        )
    output:
        join(MACS2_DIR, "{cluster}_peaks.cluster_specific.bed")
    shell:
        """
        {IONICE} intersectBed -v -a {input.a} -b {input.b} > {output}
        """

rule count_reads_per_peak:
    input:
        bam = [join(OUTBAM_DIR, i + "__atac__{cluster}.bam") 
               for i in sorted(samples)],
        bai = [join(OUTBAM_DIR, i + "__atac__{cluster}.bam.bai") 
               for i in sorted(samples)],
        bed = rules.make_master_peaks.output,
    output:
        join(COUNTS_DIR, "{cluster}_peaks" + ".min{}.txt".format(MIN_N))
    params:
        bam_order = join(COUNTS_DIR, "{cluster}_peaks.bams.txt")
    shell:
        """
        {IONICE} sortBed -i {input.bed} | \
            multiBamCov -bams {input.bam} -bed stdin > {output} && \
        for f in {input.bam}; do base=`basename $f`; echo $base; \
            done > {params.bam_order}
        """

rule summit_reproducibility:
    input:
        rules.extend_summits.output,
    output:
        join(MACS2_DIR, "{cluster}_summits.ext{ext}.reproducible"),
    params:
        glob_pattern = join(
            MACS2_IND_DIR, "*__{cluster}_summits*.noblacklist.ext{ext}.bed"
        ),
    shell:
        """
        intersectBed -wb -wa  -a {input} -b {params} | cut -f 1-4,9 | \
            sed 's/__.*//g' | sortBed | uniq > {output}
        """

rule make_master_peaks_summits:
    input:
        rules.summit_reproducibility.output
    output:
        join(
            MACS2_DIR, "{cluster}_summits.ext{ext}" + ".min{}.bed".format(MIN_N)
        )
    params:
        min_n = MIN_N
    shell:
        """
        Rscript scripts/subset_reproducible_peaks.R {input} {output} {params}
        """

rule cluster_specific_summits:
    input:
        a = rules.make_master_peaks_summits.output,
        b = lambda wildcards: expand(
            rules.make_master_peaks_summits.output,
            cluster=list_other_clusters(wildcards.cluster), ext=EXT
        )
    output:
        join(MACS2_DIR, "{cluster}_summits.ext{ext}.cluster_specific.bed")
    shell:
        """
        {IONICE} intersectBed -v -a {input.a} -b {input.b} > {output}
        """

rule count_reads_per_summit:
    input:
        bam = [join(OUTBAM_DIR, i + "__atac__{cluster}.bam") 
               for i in sorted(samples)],
        bai = [join(OUTBAM_DIR, i + "__atac__{cluster}.bam.bai") 
               for i in sorted(samples)],
        bed = rules.make_master_peaks_summits.output,
    output:
        join(COUNTS_DIR, "{cluster}_summits.ext{ext}" + ".min{}.txt".format(MIN_N))
    params:
        bam_order = join(COUNTS_DIR, "{cluster}_summits.ext{ext}.bams.txt")
    shell:
        """
        {IONICE} sortBed -i {input.bed} | \
            multiBamCov -bams {input.bam} -bed stdin > {output} && \
        for f in {input.bam}; do base=`basename $f`; echo $base; \
            done > {params.bam_order}
        """

rule trackDb_bigwig:
    input:
        expand(rules.bedgraph_to_bigwig.output, cluster=clusters),
    output:
        join(BROWSER_DIR, "trackDb.bigwig.txt"),
    params:
        BROWSER_DIR,
    shell:
        """
        url_base="https://theparkerlab.med.umich.edu"
        url_base="${{url_base}}/data/albanus/2020_nih_islets_sn_t1d"
        for f in {input}
        do
            full_path=`readlink -e ${{f}}`
            base=`basename ${{f}}`
            name=`basename ${{f}} _treat_pileup.normalized.bw`
            outlink="{params}/${{base}}"
            ln -fs ${{full_path}} ${{outlink}}

            echo "track type=bigWig db=hg19 name='${{name}}' visibility=full \
                color=255,128,0 alwaysZero=on maxHeightPixels=50:50:50 \
                windowingFunction=mean smoothingWindow=3 autoScale=off \
                viewLimits=0:20 \
                bigDataUrl=${{url_base}}/data_freeze_2021_01_11/${{base}}"
        done > {output}
        """

rule trackDb_bigwig_conditions:
    input:
        expand(
            join(
                MACS2_COND_DIR, 
                "{cond}_atac__{clust}_treat_pileup.normalized.bw"
            ),
            cond=conditions, clust=beta_cluster
        ),
    output:
        join(BROWSER_DIR, "trackDb.bigwig_conditions.txt"),
    params:
        BROWSER_DIR,
    shell:
        """
        url_base="https://theparkerlab.med.umich.edu"
        url_base="${{url_base}}/data/albanus/2020_nih_islets_sn_t1d"
        for f in {input}
        do
            full_path=`readlink -e ${{f}}`
            base=`basename ${{f}}`
            name=`basename ${{f}} _treat_pileup.normalized.bw`
            outlink="{params}/${{base}}"
            ln -fs ${{full_path}} ${{outlink}}
            echo "track type=bigWig db=hg19 name='${{name}}' visibility=full \
                color=255,128,0 alwaysZero=on maxHeightPixels=50:50:50 \
                windowingFunction=mean smoothingWindow=3 autoScale=off \
                viewLimits=0:20 \
                bigDataUrl=${{url_base}}/data_freeze_2021_01_11/${{base}}"
        done > {output}
        """

rule trackDb_peaks:
    input:
        expand(
            join(MACS2_DIR, "{cluster}_peaks.min{MIN_N}.bed"),
            cluster=clusters, MIN_N=MIN_N
        )
    output:
        join(BROWSER_DIR, "trackDb.peaks.txt"),
    params:
        BROWSER_DIR,
    shell:
        """
        url_base="https://theparkerlab.med.umich.edu"
        url_base="${{url_base}}/data/albanus/2020_nih_islets_sn_t1d"
        for f in {input}
        do
            name=`basename ${{f}} _peaks.min2.bed`
            base=`basename ${{f}}`
            outfile="{params}/${{base}}"
            echo "track name='${{name}} peaks'" > ${{outfile}}
            cut -f 1-3 ${{f}} >> ${{outfile}}
            echo "${{url_base}}/data_freeze_2021_01_11/${{base}}" 
        done > {output}
        """

rule trackDb_summits:
    input:
        expand(
            join(MACS2_DIR, "{cluster}_summits.ext150.min{MIN_N}.bed"),
            cluster=clusters, MIN_N=MIN_N
        )
    output:
        join(BROWSER_DIR, "trackDb.summits150.txt"),
    params:
        BROWSER_DIR,
    shell:
        """
        url_base="https://theparkerlab.med.umich.edu"
        url_base="${{url_base}}/data/albanus/2020_nih_islets_sn_t1d"
        for f in {input}
        do
            name=`basename ${{f}} _summits.ext150.min2.bed`
            base=`basename ${{f}}`
            outfile="{params}/${{base}}"
            echo "track name='${{name}} summits'" > ${{outfile}}
            cut -f 1-3 ${{f}} >> ${{outfile}}
            echo "${{url_base}}/data_freeze_2021_01_11/${{base}}" 
        done > {output}
        """