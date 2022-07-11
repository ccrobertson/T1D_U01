#!/usr/bin/env python

# Write a script to read simple list of sample names from a file,
# create a nested python dictionary (library -> readgroups -> fastq files)
# and then write to stdout in json format to be used as the nextflow -params-file


import argparse
import json
import pandas as pd
import os


parser = argparse.ArgumentParser()
parser.add_argument('--demuxlet_report', dest="demuxlet_report", default = None, required = True, help = 'Demuxlet report provided by seqcore.')
parser.add_argument('--fastq_dir', dest="dir", default = None, required = True, help = 'Location of fastq files.')
parser.add_argument('--modality', dest="modality", default = None, required = True, help = 'Create json for ATAC or GEX data?')

args = parser.parse_args()

d = pd.read_csv(args.demuxlet_report)
samples = d["Sample_ID"]

#for each sample create dictionary
libraries = {}
for sample_id in samples:
    #sample_id = str.strip(x)
    libraries[sample_id] = {}
    libraries[sample_id]["genome"] = ["hg38"]
    libraries[sample_id]["readgroups"] = {}
    if args.modality == "GEX":
        #create GEX dictionary
        for l in [ 'a' ]:  ### Note this is just a dummy value to accommodate expected json format expected by nf pipeline
            group = sample_id + "_" + l
            libraries[sample_id]["readgroups"][group] = {}
            r1=args.dir + "/" + sample_id + "_R1_001.fastq.gz"
            r2=args.dir + "/" + sample_id + "_R2_001.fastq.gz"
            if os.path.isfile(r1) and os.path.isfile(r2):
                libraries[sample_id]["readgroups"][group]["1"] = r1
                libraries[sample_id]["readgroups"][group]["2"] = r2
            else:
                print("ERROR: Missing expected GEX fastq files.")
                exit()
    else:
        print("ERROR: invalid modality specified.")
        exit()

#print to stdout in json format
print(json.dumps({"libraries": libraries}, indent=4))



# #for each sample create dictionary
# atac_libraries = {}
# gex_libraries = {}
# for x in f:
#     #print(x)
#     sample_id = str.strip(x)
#     #create ATAC dictionary
#     atac_libraries[sample_id] = {}
#     atac_libraries[sample_id]["genome"] = ["hg38"]
#     atac_libraries[sample_id]["readgroups"] = {}
#     for l in [ 'a', 'b', 'c', 'd']:
#         group = sample_id + "_" + l
#         atac_libraries[sample_id]["readgroups"][group] = {}
#         index=args.atacdir + "/" + sample_id + "-ATAC_" + l + "_R2_001.fastq.gz"
#         pe1=args.atacdir + "/" + sample_id + "-ATAC_" + l + "_R1_001.fastq.gz"
#         pe2=args.atacdir + "/" + sample_id + "-ATAC_" + l + "_R3_001.fastq.gz"
#         if os.path.isfile(index) and os.path.isfile(pe1) and os.path.isfile(pe2):
#             atac_libraries[sample_id]["readgroups"][group]["index"] = index
#             atac_libraries[sample_id]["readgroups"][group]["1"] = pe1
#             atac_libraries[sample_id]["readgroups"][group]["2"] = pe2
#         else:
#             print("ERROR: Missing expected ATAC fastq files.")
#             exit()
#     #create GEX dictionary
#     gex_libraries[sample_id] = {}
#     gex_libraries[sample_id]["genomes"] = ["hg38"]
#     gex_libraries[sample_id]["readgroups"] = {}
#     for l in [ 'a' ]:  ### Note this is just a dummy value to accommodate expected json format expected by nf pipeline
#         group = sample_id + "_" + l
#         gex_libraries[sample_id]["readgroups"][group] = {}
#         r1=args.gexdir + "/" + sample_id + "-3GEX_R1_001.fastq.gz"
#         r2=args.gexdir + "/" + sample_id + "-3GEX_R2_001.fastq.gz"
#         print(r1)
#         print(r2)
#         if os.path.isfile(r1) and os.path.isfile(r2):
#             gex_libraries[sample_id]["readgroups"][group]["1"] = r1
#             gex_libraries[sample_id]["readgroups"][group]["2"] = r2
#         else:
#             print("ERROR: Missing expected GEX fastq files.")
#             exit()
#
#
#
# #save to json file
# print(json.dumps({"libraries": atac_libraries}, indent=4))
# print(json.dumps({"libraries": gex_libraries}, indent=4))
