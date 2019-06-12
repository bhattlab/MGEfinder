import warnings
warnings.filterwarnings("ignore")
import sys
import click
from mgefinder import bwatools, samtools
from os.path import dirname, basename, join, isfile
import pysam
from snakemake import shell
import pygogo as gogo

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


def read_sam_pairs(samfile):
    while True:
        yield(next(samfile), next(samfile))

def format_for_mgefinder(in_sam, out_bam, read_format='rb', delete_in_sam=False):

    samfile = pysam.AlignmentFile(in_sam, read_format)
    outbam= pysam.AlignmentFile(out_bam, "wb", template=samfile)
    sam_pairs = read_sam_pairs(samfile)

    for p1, p2 in sam_pairs:

        p1_qual = p1.tostring(samfile).split('\t')[10]
        p2_qual = p2.tostring(samfile).split('\t')[10]

        if p1.query_name != p2.query_name:

            click.echo("Fatal Error: Mismatched reads in bwa_tools.format()")
            sys.exit()

        if not p1.is_unmapped and p2.is_unmapped:
            p1.set_tag('MT', p2.query_sequence)
            p1.set_tag('MQ', p2_qual)
            outbam.write(p1)
        elif p1.is_unmapped and not p2.is_unmapped:
            p2.set_tag('MT', p1.query_sequence)
            p2.set_tag('MQ', p1_qual)
            outbam.write(p2)
        elif not p1.is_unmapped and not p2.is_unmapped:
            p1.set_tag('MT', p2.query_sequence)
            p1.set_tag('MQ', p2_qual)
            p2.set_tag('MT', p1.query_sequence)
            p2.set_tag('MQ', p1_qual)
            outbam.write(p1)
            outbam.write(p2)

    samfile.close()
    outbam.close()
    if delete_in_sam:
        shell('rm {in_sam}'.format(in_sam=in_sam))

    if isfile(out_bam):
        return out_bam
    else:
        return None

def _formatbam(in_sam, out_bam, single_end, keep_tmp_files):

    logger.info("Removing secondary alignments...")
    tmp_cleaned_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.cleaned.bam.tmp')
    sorted_query_name = samtools.remove_secondary_alignments(in_sam, tmp_cleaned_bam, delete_in_bam=False)
    if not sorted_query_name:
        logger.info("Fatal error: Failed to remove secondary alignments...")
        sys.exit()
    logger.info("Successfully removed secondary alignments...\n")

    if not single_end:
        logger.info("Formatting bam file for use by mgefinder...")
        tmp_formatted_bam = join(dirname(out_bam), '.'.join(basename(out_bam).split('.')[:-1]) + '.formatted.bam.tmp')
        reformatted = format_for_mgefinder(tmp_cleaned_bam, tmp_formatted_bam, delete_in_sam=not keep_tmp_files)
        if not reformatted:
            logger.error("Fatal error: SAM file reformatting failed.")
            sys.exit()
        logger.info("SAM file successfully reformatted...\n")
    else:
        tmp_formatted_bam = tmp_cleaned_bam

    logger.info("Sorting the BAM file by chromosomal location...")
    sorted_bam = samtools.sort_coordinate(tmp_formatted_bam, out_bam, delete_in_bam=not keep_tmp_files)
    if not sorted_bam:
        logger.error("Fatal error: Failed to sort the BAM file.")
        sys.exit()
    logger.info("BAM file successfully sorted...\n")

    logger.info("Index the sorted BAM file...")
    indexed = samtools.index(out_bam)
    if not indexed:
        logger.info("Fatal error: Failed to index sorted BAM file")
        sys.exit()
    logger.info("BAM file successfully indexed...\n")


if __name__ == '__main__':
    _formatbam()