#!/usr/bin/env python3

from os.path import join
import os
from functools import partial


_results = partial(os.path.join, "results/five_prime")
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")


#Step 1 - fastqc
#rule fastqc:


#Step 2 - alignment -- for now using cellranger output
rule alignment:
    input:
        insert_fastqs
        barcode_fastqs
    ouput:
        bam
        counts
    params:
        sample=
        star_index=
        gtf_file=
        soloUMIlen=
        barcode_whitelist=
    shell:
        """
        STAR --soloBarcodeReadLength 0 \
            --runThreadN 10 \
            --outFileNamePrefix {params.sample}. \
            --genomeLoad NoSharedMemory \
            --runRNGseed 789727 \
            --readFilesCommand gunzip -c \
            --outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN \
            --genomeDir {params.star_index} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within KeepPairs \
            --sjdbGTFfile {params.gtf_file} \
            --soloType Droplet \
            --soloUMIlen {params.soloUMIlen} \
            --soloFeatures Transcript3p Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS SJ Velocyto \
            --soloMultiMappers Uniform PropUnique EM Rescue \
            --soloUMIfiltering MultiGeneUMI \
            --soloCBmatchWLtype 1MM_multi_pseudocounts \
            --soloCellFilter None \
            --soloCBwhitelist {params.barcode_whitelist} \
            --readFilesIn {insert_fastq} {barcode_fastq}
        """




#Step 3 - SCAFE
#rule scafe:
#    input:

### What is the fragment length distribution?
# cannot currently align R1
# what proportion of R2 include TSO sequence -- meaning we sequenced through the insert?
# oligo sequence is available here:  https://support.10xgenomics.com/permalink/6iTwTmSBVBPTYMXt3oVmD
# The oligo sequences start on page 94.
