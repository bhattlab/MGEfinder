import sys
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import pandas as pd
from mgefinder import fastatools
from mgefinder import bowtie2tools
from mgefinder import sctools
from mgefinder import misc
from mgefinder import pysamtools
import click
import pysam
from Bio import SeqIO
from collections import OrderedDict
from os.path import dirname, join
from random import randint
from snakemake import shell
from collections import defaultdict


def _inferseq_database(pairsfile, inferseq_database, min_perc_identity, max_internal_softclip_prop,
                       max_edge_distance, output_file, keep_intermediate):

    index_database(inferseq_database)

    database_dict = {rec.id: rec.seq for rec in SeqIO.parse(inferseq_database, 'fasta')}

    tmp_dir = dirname(output_file)

    pairs = pd.read_csv(pairsfile, sep='\t', keep_default_na=False, na_values=[
        '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A','N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan'])

    handle_empty_pairsfile(pairs, output_file)

    click.echo("Aligning pairs to database...")
    assembly_termini_fasta_prefix = write_termini_to_align_to_database(pairs, tmp_dir)
    assembly_outbam = join(tmp_dir, 'mgefinder.inferseq_database.' + str(randint(0, 1e20)) + '.bam')
    bowtie2tools.align_fasta_to_genome(
        assembly_termini_fasta_prefix + '.fasta',
        inferseq_database, assembly_outbam, silence=True,
        additional_flags='--all --score-min G,1,5'
    )

    click.echo("Inferring sequences from pairs aligned to database...")
    sequences_inferred_database = infer_sequences_database(assembly_outbam, database_dict,
                                                           min_perc_identity, max_internal_softclip_prop, max_edge_distance)

    if not keep_intermediate:
        shell('rm {fasta_prefix}* {outbam}*'.format(fasta_prefix=assembly_termini_fasta_prefix, outbam=assembly_outbam))

    method1 = make_dataframe(sequences_inferred_database, method='inferred_database')

    all_inferred_results = method1.sort_values(
        by=['pair_id', 'method']
    )

    all_inferred_results.loc[:, 'pair_id'] = list(map(str, map(int, list(all_inferred_results['pair_id']))))
    all_inferred_results = all_inferred_results.query("inferred_seq_length > 0")

    click.echo("Writing results to file %s..." % output_file)

    if not output_file:
        output_file = 'mgefinder.inferseq_database.tsv'

    if pairs.shape[0] > 0:
        sample_id = list(pairs['sample'])[0]
        all_inferred_results.insert(0, 'sample', sample_id)
    else:
        all_inferred_results.insert(0, 'sample', None)

    all_inferred_results.to_csv(output_file, sep='\t', index=False)


def get_inferred_sequences(pairs, genome_dict, add_softclipped_bases=False):

    inferred_sequences = []
    for read1, read2 in pairs:

        name = read1.reference_name + ':' + str(read1.reference_start) + '-' + str(read2.reference_end)
        inferred_sequence = genome_dict[read1.reference_name][read1.reference_start:read2.reference_end]

        if add_softclipped_bases:
            inferred_sequence = sctools.left_softclipped_sequence_strict(read1) + inferred_sequence + sctools.right_softclipped_sequence_strict(read2)

        if read1.is_read2:
            inferred_sequence = misc.revcomp(inferred_sequence)

        inferred_sequences.append((name, len(inferred_sequence), inferred_sequence))

    return inferred_sequences


def make_dataframe(inferred_sequences, method=None):
    outdict = OrderedDict([("pair_id", []), ("method", []), ("loc", []),
                           ("inferred_seq_length", []), ("inferred_seq", [])])
    for pair_id in inferred_sequences:
        for result in inferred_sequences[pair_id]:
            outdict['pair_id'].append(int(pair_id.split('_')[0]))
            outdict['method'].append(method)
            outdict['loc'].append(str(result[0]))
            outdict['inferred_seq_length'].append(result[1])
            outdict['inferred_seq'].append(''.join(result[2]))

    outdf = pd.DataFrame.from_dict(outdict)

    return outdf


def infer_sequences_database(bam_file, database_dict, min_perc_identity, max_internal_softclip_prop, max_edge_distance):
    bam = pysam.AlignmentFile(bam_file, 'rb')

    keep_reads = prefilter_reads(bam, database_dict, min_perc_identity, max_internal_softclip_prop, max_edge_distance)
    keep_pairs = get_pairs(keep_reads, database_dict, max_edge_distance)

    print("TOTAL PAIRS AFTER FILTER 1: %d" % count_total_pairs(keep_pairs))

    for pair_id in keep_pairs:
        keep_pairs[pair_id] = keep_best_alignment_score(keep_pairs[pair_id])
    print("TOTAL PAIRS AFTER FILTER 2: %d" % count_total_pairs(keep_pairs))

    inferred_sequences = defaultdict(list)
    for pair in keep_pairs:
        inferred_sequences[pair] = get_inferred_sequences(keep_pairs[pair], database_dict, add_softclipped_bases=True)

    return inferred_sequences


def count_total_pairs(pairs):
    count = 0
    for pair in pairs:
        count += len(pairs[pair])
    return count


def prefilter_reads(bam, database_dict, min_perc_identity, max_internal_softclip_prop, max_edge_distance):
    keep_reads = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for read in bam:


        if pysamtools.get_perc_identity(read) < min_perc_identity:
            continue

        if not read.is_reverse:

            if sctools.is_right_softclipped_strict(read) and \
                sctools.right_softclipped_position(read) < len(database_dict[read.reference_name]) and \
                sctools.right_softclip_proportion(read) > max_internal_softclip_prop:
                continue

            elif read.reference_start > max_edge_distance:
                continue

            elif sctools.is_left_softclipped_strict(read) and \
                abs(0 - sctools.left_softclip_reference_start(read)) > max_edge_distance:
                continue

        if read.is_reverse:

            if sctools.is_left_softclipped_strict(read) and \
                sctools.left_softclipped_position(read) >= 0 and \
                sctools.left_softclip_proportion(read) > max_internal_softclip_prop:
                continue

            elif (len(database_dict[read.reference_name]) - read.reference_end) > max_edge_distance:
                continue

            elif sctools.is_right_softclipped_strict(read) and \
                abs(0 - (len(database_dict[read.reference_name]) - sctools.right_softclip_reference_end(read))) > max_edge_distance:
                continue

        pair_id, terminus_id = read.query_name.split('_')

        keep_reads[pair_id][read.reference_name][terminus_id].append(read)

    return keep_reads


def get_pairs(reads, database_dict, length_difference_max):

    keep_reads = defaultdict(list)
    for pair in reads:

        for ref in reads[pair]:

            if len(reads[pair][ref]['1']) > 0 and \
                len(reads[pair][ref]['2']) > 0 and \
                reads_mapped_both_ends(reads[pair][ref]['1'], reads[pair][ref]['2'], database_dict, length_difference_max):

                pairs = match_pairs(reads[pair][ref]['1'], reads[pair][ref]['2'])


                for p in pairs:

                    if not pair_mapped_both_ends(p, database_dict, length_difference_max):
                        continue

                    if p[0].is_reverse and (not p[1].is_reverse):
                        keep_reads[pair].append((p[1], p[0]))

                    elif (not p[0].is_reverse) and p[1].is_reverse:
                        keep_reads[pair].append((p[0], p[1]))

    return keep_reads


def reads_mapped_both_ends(reads1, reads2, database_dict, length_difference_max):
    fiveprime_read1_count = 0
    fiveprime_read2_count = 0
    threeprime_read1_count = 0
    threeprime_read2_count = 0


    for read in reads1 + reads2:
        if read.is_reverse and (len(database_dict[read.reference_name]) - read.reference_end) <= length_difference_max:

            if read.query_name.split('_')[-1] == '1':
                threeprime_read1_count += 1
            else:
                threeprime_read2_count += 1
        elif not read.is_reverse and read.reference_start <= length_difference_max:
            if read.query_name.split('_')[-1] == '1':
                fiveprime_read1_count += 1
            else:
                fiveprime_read2_count += 1

    if (fiveprime_read1_count > 0 and threeprime_read2_count > 0) or \
        (threeprime_read1_count > 0 and fiveprime_read2_count > 0):
        return True
    else:
        return False

def pair_mapped_both_ends(pair, database_dict, length_difference_max):
    has_fiveprime_read1 = False
    has_fiveprime_read2 = False
    has_threeprime_read1 = False
    has_threeprime_read2 = False

    read1, read2 = pair


    if read1.is_reverse and (len(database_dict[read1.reference_name]) - read1.reference_end) <= length_difference_max:

        has_threeprime_read1 = True

    elif not read1.is_reverse and read1.reference_start <= length_difference_max:

        has_fiveprime_read1 = True


    if read2.is_reverse and (len(database_dict[read2.reference_name]) - read2.reference_end) <= length_difference_max:

        has_threeprime_read2 = True

    elif not read2.is_reverse and read2.reference_start <= length_difference_max:

        has_fiveprime_read2 = True


    if (has_fiveprime_read1 and has_threeprime_read2) or \
            (has_threeprime_read1 and has_fiveprime_read2):
        return True
    else:
        return False


def match_pairs(reads1, reads2):
    pairwise_dist = []
    for i in range(len(reads1)):
        read1 = reads1[i]
        if read1.is_reverse:
            read1_pos = read1.reference_end
        else:
            read1_pos = read1.reference_start
        for j in range(len(reads2)):
            read2 = reads2[j]

            if read2.is_reverse:
                read2_pos = read2.reference_end
            else:
                read2_pos = read1.reference_start

            pairwise_dist.append((i, j, abs(read1_pos - read2_pos)))

    pairwise_dist_sorted = sorted(pairwise_dist, key=lambda x: x[2], reverse=True)
    pairs = []

    paired1 = set()
    paired2 = set()

    for s in pairwise_dist_sorted:
        if s[0] not in paired1 and s[1] not in paired2:
            pairs.append((reads1[s[0]], reads2[s[1]]))
            paired1.add(s[0])
            paired2.add(s[1])

    return pairs


def keep_best_alignment_score(reads):
    keep_reads = []

    best_score = 0
    for read1, read2 in reads:
        read_combined_alignment_score = read1.get_tag('AS') + read2.get_tag('AS')
        if read_combined_alignment_score > best_score:
            best_score = read_combined_alignment_score

    for read1, read2 in reads:
        read_combined_alignment_score = read1.get_tag('AS') + read2.get_tag('AS')
        if read_combined_alignment_score == best_score:
            keep_reads.append((read1, read2))

    return keep_reads


def index_database(inferseq_database):

    if not bowtie2tools.genome_is_indexed(inferseq_database):
        click.echo("Indexing inferseq database...")
        bowtie2tools.index_genome(inferseq_database)
    click.echo("Database has been indexed...")


def write_termini_to_align_to_database(pairs, tmp_dir):
    fasta_prefix = join(tmp_dir, 'mgefinder.inferseq_database.' + str(randint(0, 1e20)))

    writetermini = get_termini(pairs)

    fastatools.write_termini_to_unpaired_fasta(writetermini, fasta_prefix)

    return fasta_prefix


def get_termini(pairs):
    termini = []
    for index, p in pairs.iterrows():
        pair_id, seq_5p, seq_3p = p['pair_id'], p['seq_5p'], p['seq_3p']
        seq_5p, seq_3p = seq_5p.upper().rstrip('N'), seq_3p.upper().lstrip('N')
        termini.append({'pair_id':str(pair_id),
                       'seq_5p': seq_5p,
                       'seq_3p': seq_3p})
    return termini


def handle_empty_pairsfile(pairs, output_file):
    if pairs.shape[0] == 0:
        outfile = pd.DataFrame(columns=['pair_id', 'method', 'loc', 'inferred_seq_length', 'inferred_seq'])

        if not output_file:
            output_file = 'mgefinder.inferseq_database.tsv'

        outfile.to_csv(output_file, sep='\t', index=False)
        click.echo("Empty pairs file, exiting...")
        sys.exit()