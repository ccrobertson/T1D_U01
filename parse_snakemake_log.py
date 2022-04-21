#!/usr/bin/env python

import sys
import os
import re
import pandas as pd
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_colwidth', 100)
from shutil import get_terminal_size
pd.set_option('display.width', get_terminal_size()[0])

log = sys.argv[1]
#log = "logs/snakemake.log"


f = open(log)
submissions = {}
count = 0
for x in f:
    #Is this a new instance?
    if re.search("^\[(.+)\]$", x):
        count += 1
        submissions[count] = {}
        m = re.search("^\[(.+)\]$", x)
        submissions[count]["date"] = m.group(1)

    #Get job info
    if re.search("jobid: ([0-9]+)", x):
        m = re.search("jobid: ([0-9]+)", x)
        submissions[count]["jobid"] = m.group(1)
    if re.search("rule (\w+):", x):
        m = re.search("rule (\w+):", x)
        submissions[count]["rule"] = m.group(1)
    if re.search("wildcards: (\w+)", x):
        m = re.search("wildcards: (.+)", x)
        submissions[count]["wildcards"] = m.group(1)

    #Get submission info
    if re.search("Submitted job ([0-9]+) with external jobid 'Submitted batch job ([0-9]+)'", x):
        m = re.search("Submitted job ([0-9]+) with external jobid 'Submitted batch job ([0-9]+)'", x)
        submissions[count]["slurmid"] = m.group(2)
    if re.search("Error", x):
        submissions[count]["slurmid"] = "Error"

    #Get finished job info
    if re.search("Finished job [0-9]+.", x):
        m = re.search("Finished job ([0-9]+).", x)
        submissions[count]["jobid"] = m.group(1)
        submissions[count]["slurmid"] = "Finished"




out = pd.DataFrame.from_dict(submissions, orient='index')
if "slurmid" in list(out.columns):
    print(out[["date","rule","wildcards","jobid","slurmid"]])
else:
    print(out[["date","rule","wildcards","jobid"]])
