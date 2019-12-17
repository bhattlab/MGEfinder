FROM continuumio/miniconda3:4.6.14

RUN conda install -c bioconda python=3.6 htslib pysam snakemake emboss bowtie2 samtools cd-hit

RUN pip install mgefinder

WORKDIR /data/
