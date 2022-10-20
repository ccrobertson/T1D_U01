#!/usr/bin/env python

import argparse
import json
import pandas as pd
import glob
import re
import os


parser = argparse.ArgumentParser()
parser.add_argument('--fastq_dir', dest="fastq_dir", default = None, required = True, help = 'Location of fastq files.')
parser.add_argument('--library_key', dest="library_key", default = None, required = True, help = 'TSV file mapping Seqcore_SampleID in bcl2fastq report to Master_SampleID for analysis pipeline.')
parser.add_argument('--modality', dest="modality", default = None, required = True, help = 'Which modality pipeline is this for?')
args = parser.parse_args()


### Map SeqcoreIDs to Master_SampleIDs
library_key = pd.read_table(args.library_key)

def map_SeqcoreID_to_SampleID(seqcore_id):
    master_id = library_key.loc[library_key['Seqcore_SampleID']==seqcore_id]['Master_SampleID'].tolist()[0]
    return master_id

def map_SampleID_to_SeqcoreID(master_id):
    seqcore_id = library_key.loc[library_key['Master_SampleID']==master_id]['Seqcore_SampleID'].tolist()[0]
    return seqcore_id


### Get unique values in a list
def unique(l):
    ul = []
    for x in l:
        if x not in ul:
            ul.append(x)
    return ul



### Which samples to include in json?
seqcore_readgroup_ids = bcl2fq_report['Sample_ID']
if args.modality == 'multiome_ATAC':
    seqcore_sample_ids = unique([re.sub("-ATAC_[abcd]","", x) for x in seqcore_readgroup_ids])
elif args.modality == 'multiome_3GEX':
    seqcore_sample_ids = unique([re.sub("-3GEX","", x) for x in seqcore_readgroup_ids])
else:
    seqcore_sample_ids = unique(seqcore_readgroup_ids)



### For each seqcore sample create dictionary
libraries = {}
for seqcore_id in seqcore_sample_ids:
    master_id = map_SeqcoreID_to_SampleID(seqcore_id)
    libraries[master_id] = {}
    libraries[master_id]["genome"] = ["hg38"]
    if args.modality == 'ATAC':
        libraries[master_id]["readgroups"] = {}
        for l in [ 'a', 'b', 'c', 'd']:
            group = seqcore_id + "-ATAC_" + l
            libraries[master_id]["readgroups"][group] = {}
            dir = args.fastq_dir + "/" + group
            index=glob.glob(dir + "/*_R2_001.fastq.gz")[0]
            pe1=glob.glob(dir + "/*_R1_001.fastq.gz")[0]
            pe2=glob.glob(dir + "/*_R3_001.fastq.gz")[0]
            if os.path.isfile(index) and os.path.isfile(pe1) and os.path.isfile(pe2):
                libraries[master_id]["readgroups"][group]["index"] = index
                libraries[master_id]["readgroups"][group]["1"] = pe1
                libraries[master_id]["readgroups"][group]["2"] = pe2
            else:
                print("ERROR: Missing expected ATAC fastq files.")
                exit()
    elif args.modality == "3GEX":
        libraries[master_id]["readgroups"] = {}
        for l in [ 'a' ]:  ### Note this is just a dummy value to accommodate expected json format expected by nf pipeline
            group = master_id + "_" + l
            libraries[master_id]["readgroups"][group] = {}
            R1 = glob.glob(args.fastq_dir + "/*" + seqcore_id + "-3GEX/*" + "R1_001.fastq.gz")[0]
            R2 = glob.glob(args.fastq_dir + "/*" + seqcore_id + "-3GEX/*" + "R2_001.fastq.gz")[0]
            if os.path.isfile(R1) and os.path.isfile(R2):
                libraries[master_id]["readgroups"][group]["1"] = R1
                libraries[master_id]["readgroups"][group]["2"] = R2
            else:
                print("ERROR: Missing expected GEX fastq files.")
                exit()
    else:
        print("ERROR: invalid modality specified.")
        exit()

#print to stdout in json format
print(json.dumps({"libraries": libraries}, indent=4))
