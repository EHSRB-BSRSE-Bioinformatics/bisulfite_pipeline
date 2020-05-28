import re
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
fastqc_dir = output_dir / "fastqc"
align_dir = output_dir / "aligned"
procesed_dir = output_dir / "processed"
log_dir = output_dir / "logs"

output_dir.mkdir(parents=True, exist_ok=True)
trim_dir.mkdir(parents=True, exist_ok=True)
fastqc_dir.mkdir(parents=True, exist_ok=True)
align_dir.mkdir(parents=True, exist_ok=True)
procesed_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

print("using genome: " + str(genome_dir))

print("input dir: " + str(input_dir))
print("output dir: " + str(output_dir))

print("trimmed sequences dir: " + str(trim_dir))
print("fastqc dir: " + str(fastqc_dir))
print("alignments dir: " + str(align_dir))
print("snakemake log dir: " + str(log_dir))

shell("cp -rf Snakefile {log_dir}")
shell("cp -rf {config_file_name} {log_dir}")

# gene names file, needed for annotation

rule gene_names:
    input: genome_dir / config["annotation_filename"]
    output: genome_dir / "geneid_to_name.txt"
    run:
        shell("zcat {input} | awk 'BEGIN{FS=\"\t\"}{split($9,a,\";\"); if($3~\"gene\") print) a[1]\"\t\"a[3]}' | sed 's/gene_source \"ensembl\"//' | tr ' ' '\t' | cut -f 2,5 | tr -d '\"' > {output}")

# run whole pipeline

rule all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),
        output_dir / "multiqc_report.html",
    log: log_dir / "all.log"

# sequence trimming

rule trim_galore_all:
    input:
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES)
    log: log_dir / "trim_galore_all.log"

rule trim_galore:
    input:
        input_dir / "{sample}.fastq.gz"
    output: 
        trim_dir / "{sample}_trimmed.fq.gz"
    log: log_dir / "trim_galore.{sample}.log"
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
    log: log_dir / "bismark_index.log"
    shell:
        "bismark_genome_preparation {genome_dir}"

rule bismark_all:
    message: "running bismark on all samples"
    input:
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES),
    log: log_dir / "bismark_all.log"

rule bismark:
    input:
        ancient(genome_dir / "Bisulfite_Genome"),
        pair1 = str(trim_dir / "{sample}_trimmed.fq.gz")
    output:
        align_dir / "{sample}_trimmed_bismark_bt2.bam"
    log: log_dir / "bismark.{sample}.log"
    run:
        shell("bismark --output_dir {align_dir} {genome_dir} {input.pair1}")

rule bismark_extract_methylation:
    input:
        align_dir / "{sample}_trimmed_bismark_bt2.bam"
    output:
        procesed_dir / "{sample}_trimmed_bismark_bt2.CpG_report.txt.gz"
    log: log_dir / "bismark_methylextract.{sample}.log"
    run:
        shell("bismark_methylation_extractor -s --gazillion --bedGraph --buffer_size 10G --cytosine_report --genome_folder {genome_dir} --gzip --output {procesed_dir} {input}")

# fastqc

rule fastqc_all:
    input:
        expand(str(fastqc_dir / "{sample}_fastqc.zip"), sample=SAMPLES)

rule fastqc:
    input:
        input_dir / "{sample}.fastq.gz"
    output: 
        fastqc_dir / "{sample}_fastqc.zip"
    run:
        shell("fastqc -o {fastqc_dir} {input}")


# multiqc

rule multiqc:
    message: "running multiqc"
    input:
        expand(str(fastqc_dir / "{sample}_fastqc.zip"), sample=SAMPLES),
        expand(str(trim_dir / "{sample}_trimmed.fq.gz"), sample=SAMPLES),
        expand(str(align_dir / "{sample}_trimmed_bismark_bt2.bam"), sample=SAMPLES)
    output:
        output_dir / "multiqc_report.html"
    log: log_dir / "multiqc.log"
    run:
        shell("multiqc -fz {input_dir} {trim_dir} {align_dir} {fastqc_dir}")
        shell("mv multiqc_report.html multiqc_data.zip {output_dir}")

# stats via methylkit

rule methylkit:
    input:
        genome_dir / "geneid_to_name.txt",
        expand(str(procesed_dir / "{sample}_trimmed_bismark_bt2.CpG_report.txt.gz"), sample=SAMPLES)
    log: log_dir / "methylkit.log"
    script:
        "./run_methylkit.R"