import warnings
warnings.filterwarnings("ignore")
from mgefinder.sctools import *

def get_left_softclipped_reads_at_site(bam_file, contig, left_site, get_quals=False, softclip_only=False):
    left_softclipped_reads = []
    left_softclipped_quals = []

    start = left_site + 1
    end = left_site + 2

    contig_len = contig_length(bam_file, contig)
    if end > contig_len:
        end = contig_len

    for pu in bam_file.pileup(contig, start, end, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_left_softclipped_lenient_at_site(read, contig, left_site):

                if softclip_only:
                    keep_seq = left_softclipped_sequence(read)
                else:
                    keep_seq = read.query_sequence
                    keep_quals = get_query_qualities_ascii(read, bam_file)

                left_softclipped_reads.append(keep_seq)
                left_softclipped_quals.append(keep_quals)

    if get_quals:
        return left_softclipped_reads, left_softclipped_quals
    else:
        return left_softclipped_reads



def get_right_softclipped_reads_at_site(bam_file, contig, right_site, get_quals=False, softclip_only=False):
    right_softclipped_reads = []
    right_softclipped_quals = []

    start = right_site - 1
    end = right_site

    if start < 0:
        start = 0

    for pu in bam_file.pileup(contig, start, end, truncate=True):

        for pr in pu.pileups:
            read = pr.alignment
            if is_right_softclipped_lenient_at_site(read, contig, right_site):
                if softclip_only:
                    keep_seq = right_softclipped_sequence(read)
                else:
                    keep_seq = read.query_sequence
                    keep_quals = get_query_qualities_ascii(read, bam_file)

                right_softclipped_reads.append(keep_seq)
                right_softclipped_quals.append(keep_quals)

    if get_quals:
        return right_softclipped_reads, right_softclipped_quals
    else:
        right_softclipped_reads


def get_right_unmapped_reads(bam_file, contig, right_site, get_quals=False, search_region_length=500):
    right_unmapped_reads = []
    right_unmapped_quals = []

    start = right_site - search_region_length
    end = right_site

    if start < 0:
        start = 0

    for read in bam_file.fetch(contig, start, end):

        if (not read.is_reverse) and read.mate_is_unmapped:
            right_unmapped_reads.append(read.get_tag('MT'))
            right_unmapped_quals.append(read.get_tag('MQ'))

    if get_quals:
        return right_unmapped_reads, right_unmapped_quals
    else:
        return right_unmapped_reads


def get_left_unmapped_reads(bam_file, contig, left_site, get_quals=False, search_region_length=500):
    left_unmapped_reads = []
    left_unmapped_quals = []

    start = left_site + 1
    end = left_site + 1 + search_region_length

    contig_len = contig_length(bam_file, contig)
    if end > contig_len:
        end = contig_len

    for read in bam_file.fetch(contig, start, end):

        if read.is_reverse and read.mate_is_unmapped:
            left_unmapped_reads.append(read.get_tag('MT'))
            left_unmapped_quals.append(read.get_tag('MQ'))

    if get_quals:
        return left_unmapped_reads, left_unmapped_quals
    else:
        return left_unmapped_reads


def contig_length(bam, contig):
    return dict(zip(bam.references, bam.lengths))[contig]

def count_runthrough_reads(bam, contig, site, min_qual=20, min_alignment_inner_length=21):

    if site  < 0 or site >= contig_length(bam, contig):
        return 0
    else:
        count = 0
        for read in bam.fetch(contig, site, site+1):

            if read.mapping_quality < min_qual:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            count += 1

        return count

def count_softclipped_reads(bam, contig, site, min_qual=20, min_alignment_inner_length=21):

    if site  < 0 or site >= contig_length(bam, contig):
        return 0
    else:
        count = 0
        for read in bam.fetch(contig, site+1, site+2):

            if read.mapping_quality < min_qual:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            if is_left_softclipped_lenient_at_site(read, contig, site):
                count += 1

        for read in bam.fetch(contig, site - 1, site):

            if read.mapping_quality < min_qual:
                continue

            if not read_meets_min_alignment_inner_length(read, min_alignment_inner_length):
                continue

            if is_right_softclipped_lenient_at_site(read, contig, site):
                count += 1


        return count


def get_query_qualities_ascii(read, bam):
    return read.tostring(bam).split('\t')[10]


def query_qualities_to_phred(quals):
    return [ord(q)-33 for q in quals]


def get_perc_identity(read):
    total = read.query_alignment_length
    identity = (total - read.get_tag('NM')) / total
    return identity


def get_bam_contig_dict(bam_file):
    contig_dict = {}
    for contig in bam_file.header['SQ']:
        contig_dict[contig['SN']] = contig['LN']
    return contig_dict

def get_insertion_length(position, read, reverse=False):
    refpos = read.get_reference_positions(full_length=True)

    insert_length = 0

    if not reverse:

        for pos in refpos[refpos.index(position)+1:]:

            if pos is not None:
                break
            else:
                insert_length += 1

    else:
        refpos = refpos[::-1]
        for pos in refpos[refpos.index(position)+1:]:

            if pos is not None:
                break
            else:
                insert_length += 1

    return insert_length
