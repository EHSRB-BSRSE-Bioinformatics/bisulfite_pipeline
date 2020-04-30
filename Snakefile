import re
import os
import pandas as pd
from glob import glob
import pathlib

config_file_name = "config.yaml"

configfile: config_file_name

print("config file: " + config_file_name)

# set up conditions and samples
sample_info = config["samples"]

SAMPLES = sum([i.split(" ") for i in list(sample_info.values())], [])
CONDITIONS = list(sample_info.keys())

print("conditions: " + str(CONDITIONS))
print("samples: " + str(SAMPLES))

# set up directories
root = Path("/")

genome_dir = root / config["genome_root_dir"] / config["genome"]

input_dir = Path(config["input_dir"])
output_dir = Path(config["output_dir"])

trim_dir = output_dir / "trimmed"
align_dir = output_dir / "aligned"
log_dir = output_dir / "logs"

print("using genome: " + str())

print("input dir: " + str(input_dir))
print("output dir: " + str(output_dir))

print("trimmed sequences dir: " + str(trim_dir))
print("alignments dir: " + str(align_dir))
print("snakemake log dir: " + str(log_dir))

shell("cp -f Snakefile {log_dir}")
shell("cp -f {config_file_name} {log_dir}")

# run whole pipeline

rule all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),
        output_dir / "multiqc_report.html",

# sequence trimming

rule trim_galore_all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES)

rule trim_galore:
    input:
        input_dir / "{sample}.fastq.gz"
    output: 
        trim_dir / "{sample}_trimmed.fq.gz"
    run:
        shell("trim_galore --illumina --rrbs -q 20 --output_dir {trim_dir} {input}")


# alignment with bismark

# note: if this fails with a cryptic bowtie error, try running the bowtie2-build jobs directly:
# cd config["genome_root_dir"] / config["genome"] / Bisulfite_Genome/CT_conversion
# bowtie2-build -f genome_mfa.CT_conversion.fa BS_CT --threads 2
# cd ../GA_conversion
# bowtie2-build -f genome_mfa.GA_conversion.fa BS_GA --threads 2
rule bismark_index:
    message: "running bismark_genome_preparation"
    input:
        ancient(genome_dir)
    output:
        genome_dir / "Bisulfite_Genome"
    shell:
        "bismark_genome_preparation {genome_dir}"

rule bismark_all:
    message: "running bismark on all samples"
    input:
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),

rule bismark:
    input:
        genome_dir / "Bisulfite_Genome",
        pair1 = str(trim_dir / "{sample}_trimmed.fq.gz")
    output:
        align_dir / "{sample}_trimmed_bismark_bt2.bam"
    run:
        shell("bismark --output_dir {align_dir} {genome_dir} {input.pair1}")

# multiqc

rule multiqc:
    message: "running multiqc"
    input:
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES)
    output:
        output_dir / "multiqc_report.html"
    run:
        shell("multiqc -fz {input_dir} {trim_dir} {align_dir}")
        shell("mv multiqc_report.html multiqc_data.zip {output_dir}")

# stats via methylkit

rule methylkit:
    script:
        "./run_methylkit.R"