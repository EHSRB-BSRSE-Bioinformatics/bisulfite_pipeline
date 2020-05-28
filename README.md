# Setting up the development environment

The following assumes a conda installation. As well, FastQC requires java (`java --version` to see if it's installed)

Add bioconda and conda-forge channels:

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

1. Install SRA tools

`conda install sra-tools`

2. Install Trim Galore & requirements

`conda install cutadapt`

Get FastQC here: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Unpack it and put it in (e.g.) `~/bin/FastQC`, then make it accessible:

`chmod 755 ~/bin/FastQC/fastqc`

`ln -s ~/bin/FastQC/fastqc ./fastqc`

`conda install trim-galore`

3. Install Bismark & requirements (including Bowtie2)

`conda install bowtie2`

`conda install bismark`

4. Install MultiQC

`conda install multiqc`

5. Install methylKit & genomation from Bioconductor

(in R)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylKit")
BiocManager::install("genomation")
BiocManager::install("GenomicFeatures")
```

Note that methylKit version >=1.14.1 is required. 1.14.0 in particular had a bug.

6. Install some other R packages

`install.packages('tidyverse')`
`install.packages('ggrepel')`

# Downloading data

1. Download genome in fasta format, e.g. from http://useast.ensembl.org/info/data/ftp/index.html

(The specific file is ftp://ftp.ensembl.org/pub/release-100/fasta/mesocricetus_auratus/dna/Mesocricetus_auratus.MesAur1.0.dna_sm.toplevel.fa.gz)
 
2. Download annotation in gtf format, e.g. from ftp://ftp.ensembl.org/pub/release-100/gtf/mesocricetus_auratus/Mesocricetus_auratus.MesAur1.0.100.gtf.gz

Important note: the genome and annotation should use the same chromosome/contig naming convention.


# notes

Probably needed?

`conda install samtools`

`conda install pysam`
