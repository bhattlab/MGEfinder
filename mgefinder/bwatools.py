import warnings
warnings.filterwarnings("ignore")
import sys
from glob import glob
from snakemake import shell
import click
from os.path import isfile
import pysam
from mgefinder import misc


def index_genome(genome_path, silence=True):
    if silence:
        shell('bwa index {genome_path} &> /dev/null;'.format(genome_path=genome_path))
    else:
        shell('bwa index {genome_path}'.format(genome_path=genome_path))

    return genome_is_indexed(genome_path)

def genome_is_indexed(genome_path):

    amb = genome_path + '.amb'
    ann = genome_path + '.ann'
    bwt = genome_path + '.bwt'
    pac = genome_path + '.pac'
    sa = genome_path + '.sa'

    index_files = [amb, ann, bwt, pac, sa]

    indexed = True
    files = glob(genome_path+'*')
    for f in index_files:
        if f not in files:
            indexed = False

    return indexed

def align_to_genome_pe(fastq1, fastq2, genome_path, out_sam, threads=1, verbose=False):
    if verbose:
        command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} > {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, fastq2=fastq2, out_sam=out_sam, threads=threads)
    else:
        command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} 2> /dev/null 1> {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, fastq2=fastq2, out_sam=out_sam, threads=threads)

    logger.debug("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False

def align_to_genome_se(fastq1, genome_path, out_sam, threads=1, verbose=False):
    if verbose:
        command = "bwa mem -t {threads} {genome_path} {fastq1} > {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, out_sam=out_sam, threads=threads)
    else:
        command = "bwa mem -t {threads} {genome_path} {fastq1} 2> /dev/null 1> {out_sam}; ".format(
            genome_path=genome_path, fastq1=fastq1, out_sam=out_sam, threads=threads)

    logger.debug("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False

def align_to_genome_fasta_pe(fasta1, fasta2, genome_path, out_sam, threads=1, verbose=False, additional_flags=''):
    if verbose:
        command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} {additional_flags} > {out_sam}; ".format(
            genome_path=genome_path, fastq1=fasta1, fastq2=fasta2, out_sam=out_sam, threads=threads,
            additional_flags=additional_flags)
    else:
        command = "bwa mem -t {threads} {genome_path} {fastq1} {fastq2} {additional_flags} 2> /dev/null 1> {out_sam}; ".format(
            genome_path=genome_path, fastq1=fasta1, fastq2=fasta2, out_sam=out_sam, threads=threads,
            additional_flags=additional_flags)

    logger.debug("Executing command: %s" % command)
    shell(command)

    if isfile(out_sam):
        return True
    else:
        return False


def align_fasta_to_genome(fasta, genome_path, out_bam, threads=1, silence=True, additional_flags=''):

    if silence:
        command1 = "rm -f {out_bam}.sam {out_bam}.sam.tmp; " \
                  "bwa mem {additional_flags} -t {threads} {genome_path} {fasta} 2> /dev/null 1> {out_bam}.sam;"
        command2 = "samtools view -h -F 4 {out_bam}.sam.secseq 2> /dev/null 1> {out_bam}.sam.tmp; " \
                  "samtools sort {out_bam}.sam.tmp 2> /dev/null 1> {out_bam}; " \
                  "samtools index {out_bam} 2> /dev/null;" \
                  "rm -f {out_bam}.sam {out_bam}.sam.secseq {out_bam}.sam.tmp".format(
            genome_path=genome_path, fasta=fasta, out_bam=out_bam, threads=threads, additional_flags=additional_flags)

    else:
        command1 = "rm -f {out_bam}.sam {out_bam}.sam.tmp; " \
                  "bwa mem {additional_flags} -t {threads} {genome_path} {fasta} 1> {out_sam}.sam;"
        command2 = "samtools view -h -F 4 {out_bam}.sam.secseq 1> {out_bam}.sam.tmp; " \
                  "samtools sort {out_bam}.sam.tmp 1> {out_bam}; " \
                  "samtools index {out_bam};" \
                  "rm -f {out_bam}.sam {out_bam}.sam.secseq {out_bam}.sam.tmp".format(
            genome_path=genome_path, fasta=fasta, out_bam=out_bam, threads=threads, additional_flags=additional_flags)

    #print(command)
    shell(command1)
    add_sequence_to_secondary_alignment('{out_bam}.sam'.format(out_bam=out_bam),
                                        '{out_bam}.sam.secseq'.format(out_bam=out_bam))
    shell(command2)

    if isfile(out_bam):
        return True
    else:
        return False

def add_sequence_to_secondary_alignment(sam_file_in, sam_file_out):
    outfile = open(sam_file_out, 'w')
    infile = pysam.AlignmentFile(sam_file_in, 'r')
    outfile = pysam.AlignmentFile(sam_file_out, "w", template=infile)

    current_seq = None
    current_seq_reverse = None
    for read in infile:
        if read.query_sequence is not None:
            current_seq = read.query_sequence
            current_seq_reverse = read.is_reverse
        else:
            if current_seq_reverse == read.is_reverse:
                read.query_sequence = current_seq
            else:
                read.query_sequence = misc.revcomp(current_seq)

        outfile.write(read)
    outfile.close()