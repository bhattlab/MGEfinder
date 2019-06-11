import warnings
warnings.filterwarnings("ignore")
import sys
from glob import glob
from snakemake import shell
import pygogo as gogo
from os.path import isfile

def remove_secondary_alignments(in_bam, out_bam, delete_in_bam=False):
    shell('samtools view -b -h -F 0x900 {in_bam} > {out_bam}'.format(in_bam=in_bam, out_bam=out_bam))

    if delete_in_bam:
        shell('rm {in_bam}'.format(in_bam=in_bam))

    if isfile(out_bam):
        return True
    else:
        return False


def sort_coordinate(in_bam, out_bam, delete_in_bam=False):
    shell('samtools sort {in_bam} > {out_bam}'.format(in_bam=in_bam, out_bam=out_bam))

    if delete_in_bam:
        shell('rm {in_bam}'.format(in_bam=in_bam))

    if isfile(out_bam):
        return True
    else:
        return False

def index(in_bam):
    shell('samtools index {in_bam}'.format(in_bam=in_bam))

    if isfile(in_bam+'.bai'):
        return True
    else:
        return False