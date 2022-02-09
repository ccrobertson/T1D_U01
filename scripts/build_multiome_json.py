#!/usr/bin/env python

import argparse
import json
import os


parser = argparse.ArgumentParser()
parser.add_argument('--sample-file', dest="samplefile", default = None, required = True, help = 'List of sample IDs determined by demultiplexing.')
parser.add_argument('--fastq-dir', dest="dir", default = None, required = True, help = 'Location of fastq files.')

args = parser.parse_args()

f = open(args.samplefile)
#print(f.read())

#for each sample create dictionary
libraries = {}
for x in f:
    print(x)
    libraries[x] = {}
    libraries[x]["genomes"] = ["hg38"]
    libraries[x]["readgroups"] = {}
    for l in [ 'a', 'b', 'c', 'd']:
        libraries[x]["readgroups"][l] = {}
        index=args.dir + "/" + str.strip(x) + "-ATAC_" + l + "_R2_001.fastq.gz"
        pe1=args.dir + "/" + str.strip(x) + "-ATAC_" + l + "_R1_001.fastq.gz"
        pe2=args.dir + "/" + str.strip(x) + "-ATAC_" + l + "_R3_001.fastq.gz"
        if os.path.isfile(index) and os.path.isfile(pe1) and os.path.isfile(pe2):
            libraries[x]["readgroups"][l]["index"] = index
            libraries[x]["readgroups"][l]["1"] = pe1
            libraries[x]["readgroups"][l]["2"] = pe2
        else:
            print("ERROR: Missing some fastq files.")
            exit()


#save to json file
y = json.dumps(libraries)
print(y)
