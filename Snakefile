shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

import subprocess
import csv
import os
configfile: "config.yaml"

CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())

## how many runs
if config['download']:	
	RUNS = []
	for sample in SAMPLES:
		for batch in FILES[sample].keys():
			run = FILES[sample][batch][0]
			if run not in RUNS:
				RUNS.append(run + "/RunInfo.xml")

# for the same batch, there is only one run id and flow cell name. need to check
RUN_DICT = {}
for sample in SAMPLES:
	RUN_DICT.update(FILES[sample])

FASTQ = []
for batch in RUN_DICT:
	FASTQ.append("mkfastq_stamps/" + batch + ".stamp")

COUNT = []
for sample in SAMPLES:
	COUNT.append("count_stamps/" + sample + ".stamp")


TARGETS = []
if config['download']:
	TARGETS.extend(RUNS)

TARGETS.extend(FASTQ)
TARGETS.extend(COUNT)

if config['make_seurat']:
	SEURAT_OBJ = expand("seurat_objs/{sample}_seurat.rds", sample = SAMPLES)
	TARGETS.extend(SEURAT_OBJ)

localrules: all
rule all:
    input: TARGETS

## if the bcl files are in BaseSpace, need to download
## the downloaded folder is named with the run_id


rule download_raw:
	output: "{run_id}/RunInfo.xml"
	log: "00log/downloading_runid_{run_id}.log"
	threads: 1
	params: jobname = "{run_id}"
	message: "downloading/symlink run id {wildcards.run_id} "
	run:
		if config['download']:
			shell("source activate py27")
			shell("python scripts/BaseSpaceRunDownloader_v2.py -r {wildcards.run_id}  -a {config[token]} 2> {log}")
		## if not downloaded from baseSpace, soft link the directory 
		else:
			with open(config['run_dirs']) as f:
				reader = csv.reader(f, delimiter = "\t")
				# skipt the header
				header = next(reader)
				run_dir_dict = {}
				for row in reader:
					run_id = row[0].strip()
					run_dir = row[1].strip()
					run_dir_dict[run_id] = run_dir
		# snakemake autocreate the output folder first, I need to remove it first, otherwise ln -s will link to the subfolder
		# of run_id
		cmd1 = "rm -rf {run_id}".format(run_id = wildcards.run_id)			
		cmd2 = "ln -rs {run_dir} {run_id}".format(run_dir = run_dir_dict[wildcards.run_id], run_id = wildcards.run_id)
		shell(cmd1)
		shell(cmd2)


def get_run_id(wildcards):
	return RUN_DICT[wildcards.batch_name][0] + "/RunInfo.xml"

def get_run_info(wildcards):
	return RUN_DICT[wildcards.batch_name]


# snakemake creates the folder if it does not exisit. but cellranger fails if the folder is 
# not created by itself. I asked on twitter, and Shaun Jackman answered:
# Instead of making the target a file in a subdirectory like `sample/output` use `touch sample.stamp`
# https://twitter.com/sjackman/status/1056285115393798144
# I also found an example by Devon Ryan https://github.com/maxplanck-ie/10X_snakepipe/blob/master/Snakefile
rule make_fastq:
	input: get_run_id 
	output: "mkfastq_stamps/{batch_name}.stamp"
	log: "00log/{batch_name}_cellranger_mkfastq.log"
	threads: 8
	params: run_info = get_run_info
	message: "making fastqs for {input} using {threads} threads"
	shell:
		"""
        module load bcl2fastq2
		cellranger mkfastq --id={wildcards.batch_name} \
                     --run={params.run_info[0]} \
                     --samplesheet={params.run_info[2]} \
                     --jobmode=local \
                     --localmem=30 \
                     --localcores=8 2> {log} 
        touch {output}
		"""

def get_count_input(wildcards):
	count_input = []
	for batch in FILES[wildcards.sample].keys():
		count_input.append("mkfastq_stamps/" + batch + ".stamp")
	return count_input

def get_fastq_per_sample(wildcards):
	fastqs = []
	for batch in FILES[wildcards.sample].keys():
		fastqs.append(batch + "/outs/fastq_path/" + FILES[wildcards.sample][batch][1] + "/" + wildcards.sample)
	return ",".join(fastqs)

rule count:
	input: get_count_input
	output: "count_stamps/{sample}.stamp"
	log: "00log/{sample}_cellranger_count.log"
	threads: 12
	params: 
		fastqs = get_fastq_per_sample,
		custom = config.get("cellranger_count_args", "")
	message: "counting for {wildcards.sample} from {params} using {threads} threads"
	shell:
		"""
		cellranger count --id={wildcards.sample} \
                 --transcriptome={config[transcriptome]} \
                 --fastqs={params.fastqs} \
                 --sample={wildcards.sample} \
                 --expect-cells={config[expect_cells]} \
                 {params.custom} \
                 --localcores=12 \
                 --localmem=40 2> {log}
  		touch {output}
		"""

# transcriptome: /n/regal/informatics_public/reference_genome_by_tommy/cellranger_ref/mm10-2.1.0_premrna
# the last bit is the folder containing counts
count_dir = config['transcriptome'].split('/')[-1]

# params directive takes input functions https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html

if config['make_seurat']:
	rule make_seurat_obj:
		input: "count_stamps/{sample}.stamp"
		output: "seurat_objs/{sample}_seurat.rds"
		log: "00log/{sample}_make_seurat_obj.log"
		threads: 4
		params:
			make_seurat_obj_args = config.get("make_seurat_obj_args", ""),
			cellranger_count_dir = lambda wildcards: "{sample}/outs/filtered_gene_bc_matrices/{count_dir}".format(sample = wildcards.sample, count_dir = count_dir)
		message: "making seurat object for {sample}"
		script: "scripts/make_seurat.R"


onsuccess:
	print("Success: Snakemake completed!")
	shell("mail -s 'Snakemake workflow completed: Have a beer!' {config[email]} < {log}")

onerror:
	print("Error: Snakemake aborted!")
	shell("mail -s 'Snakemake workflow aborted: Have a coffee and see log inside!' {config[email]} < {log}")
	shell("exit 1")
