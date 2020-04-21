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


Probably needed?

`conda install samtools`

`conda install pysam`