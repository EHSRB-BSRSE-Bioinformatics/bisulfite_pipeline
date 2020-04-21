import re
import os
import pandas as pd
from glob import glob
import pathlib

configfile: "config.yaml"

SAMPLES_WT = pd.read_csv("metadata.csv").query('desc == "WT"')["SRR"].values
SAMPLES_TKO = pd.read_csv("metadata.csv").query('desc == "TKO"')["SRR"].values
SAMPLES = SAMPLES_WT.tolist() + SAMPLES_TKO.tolist()

print("WT samples: ")
print(SAMPLES_WT)
print("TKO samples: ")
print(SAMPLES_TKO)

root = Path("/")

project_root = Path(config["project_dir"])

rule trim_galore_all:
    message: "running trim galore on all samples"
    input:
        expand(str(project_root / "trimmed/{sample}_1_val_1.fq.gz"), sample=SAMPLES),
        expand(str(project_root / "trimmed/{sample}_2_val_2.fq.gz"), sample=SAMPLES)

rule trim_galore:
    message: "running trim galore for sample"
    input:
        pair1 = str(project_root / "fastq" / "{sample}_1.fastq.gz"),
        pair2 = str(project_root / "fastq" / "{sample}_2.fastq.gz")
    output: 
        project_root / "trimmed" / "{sample}_1_val_1.fq.gz",
        project_root / "trimmed" / "{sample}_2_val_2.fq.gz"
    run:
        output_dir = str(project_root / "trimmed")
        shell("trim_galore --illumina --rrbs --paired -q 20 --output_dir {output_dir} {input.pair1} {input.pair2}")

# note: if this fails with a cryptic bowtie error, try running the bowtie2-build jobs directly:
# cd config["genome_root_dir"] / config["genome"] / Bisulfite_Genome/CT_conversion
# bowtie2-build -f genome_mfa.CT_conversion.fa BS_CT --threads 2
# cd ../GA_conversion
# bowtie2-build -f genome_mfa.GA_conversion.fa BS_GA --threads 2
rule bismark_index:
    message: "running bismark_genome_preparation"
    input:
        ancient(root / config["genome_root_dir"] / config["genome"] )
    shell:
        "bismark_genome_preparation {input}"

rule bismark_all:
    message: "running bismark on all samples"
    input:
        expand(str(project_root / "aligned/{sample}_1_val_1_bismark_bt2_pe.bam"), sample=SAMPLES),


rule bismark:
    message: "running bismark"
    input:
        pair1 = str(project_root / "trimmed" / "{sample}_1_val_1.fq.gz"),
        pair2 = str(project_root / "trimmed" / "{sample}_2_val_2.fq.gz"),
    output:
        project_root / "aligned" / "{sample}_1_val_1_bismark_bt2_pe.bam"
    run:
        genome_dir = str(root / config["genome_root_dir"] / config["genome"])
        output_dir = str(project_root / "aligned")
        shell("bismark --output_dir {output_dir} {genome_dir} -1 {input.pair1} -2 {input.pair2}")

# rule bowtie_se:
#     input:
#         dir=ancient("ReadData/{sample}/")
#     output:
#         "ReadData/{sample}/{sample}_Bowtie_{DB,[A-Za-z0-9]+}.sam"
#     params:
#         DB_data = config["DB_data"], DB = config["DB"], N_MISMATCH = config["n_mismatch"]
#     run:
#         infiles = glob("%s/*_1.fastq.gz" %(input)) # Only the forward read
#         shell("bowtie2 -a -N {params.N_MISMATCH} -x {params.DB_data} -U <(zcat %s) | awk -F \"\\t\" '$3 != \"*\"' > {output}" %(infiles[0]))


# rule bowtie_se_all:
#     input:
#         expand("ReadData/{sample}/{sample}_Bowtie_Sub_se_{DB}.sam", sample = SAMPLES, DB = config["DB"]),
