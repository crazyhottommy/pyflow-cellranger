#!/usr/bin/env python3


import json
import os
import csv
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("meta", help="Required. the FULL path to the tab delimited meta file containing run info")
args = parser.parse_args()

assert args.meta is not None, "please provide the path to the meta file"

FILES = defaultdict(lambda: defaultdict(list))

with open(args.meta, "r") as f:
    reader = csv.reader(f, delimiter = "\t")
    # skip the header
    header = next(reader)
    for row in reader:
        run_id = row[0].strip()
        flow_cell = row[1].strip()
        sample = row[2].strip()
        # This is name for the fastq folder
        batch_name = row[3].strip()
        csv_file = row[4].strip()
        ## now just assume the file name in the metafile contained in the fastq file path
        FILES[sample][batch_name].append(run_id)
        FILES[sample][batch_name].append(flow_cell)
        FILES[sample][batch_name].append(csv_file)

sample_num = len(FILES.keys())
print ("total {} unique samples will be processed".format(sample_num))
print ("------------------------------------------")
runs = []
for sample_name in sorted(FILES.keys()):
    batches = sorted(FILES[sample_name].keys())
    print ("{sample_name} has {n} runs: {batches}".format(sample_name = sample_name, n = len(FILES[sample_name]), batches= " ".join(batches)))
    for batch in batches:
    	run_id = FILES[sample_name][batch][0]
    	flow_cell = FILES[sample_name][batch][1]
    	csv_file = FILES[sample_name][batch][2]
    	run_info = [batch, run_id, flow_cell, csv_file]
    	if run_info not in runs:
    		runs.append([batch, run_id, flow_cell, csv_file])
print()

for run in runs:
	batch = run[0]
	run_id = run[1]
	flow_cell = run[2]
	csv_file = run[3]
	print("{batch}'s run id is {run_id} with flow cell {flow_cell} and csv file {csv_file}". \
		format(batch = batch, run_id = run_id, flow_cell = flow_cell, csv_file = csv_file))

print ("------------------------------------------")
print("check the samples.json file for the meta data info")
print()

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)
        
