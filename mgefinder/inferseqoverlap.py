import sys
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import pandas as pd
import click
from collections import OrderedDict
from collections import defaultdict
from Bio.pairwise2 import format_alignment
from Bio import pairwise2

verbose=True


def _inferseq_overlap(pairsfile, min_overlap_score, min_overlap_perc_identity, min_inferseq_size, output_file):

    pairs = pd.read_csv(pairsfile, sep='\t', keep_default_na=False, na_values=[
        '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A','N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan'])

    handle_empty_pairsfile(pairs, output_file)

    click.echo("Inferring sequences from overlap...")
    sequences_inferred_from_overlap = infer_sequences_overlap(pairs, min_overlap_score, min_overlap_perc_identity)
    method1 = make_dataframe(sequences_inferred_from_overlap, method='inferred_overlap')

    method1.loc[:, 'pair_id'] = list(map(str, map(int, list(method1['pair_id']))))
    method1 = method1.query("inferred_seq_length >= @min_inferseq_size")

    click.echo("Writing results to file %s..." % output_file)

    if pairs.shape[0] > 0:
        sample_id = list(pairs['sample'])[0]
        method1.insert(0, 'sample', sample_id)
    else:
        method1.insert(0, 'sample', None)
    method1.to_csv(output_file, sep='\t', index=False)


def make_dataframe(inferred_sequences, method=None):
    outdict = OrderedDict([("pair_id", []), ("method", []), ("loc", []),
                           ("inferred_seq_length", []), ("inferred_seq", [])])
    for pair_id in inferred_sequences:
        for result in inferred_sequences[pair_id]:
            outdict['pair_id'].append(pair_id)
            outdict['method'].append(method)
            outdict['loc'].append(str(result[0]))
            outdict['inferred_seq_length'].append(result[1])
            outdict['inferred_seq'].append(''.join(result[2]))

    outdf = pd.DataFrame.from_dict(outdict)

    return outdf


def infer_sequences_overlap(pairs, min_overlap_score, min_overlap_perc_identity):
    sequences = defaultdict(list)
    for index, row in pairs.iterrows():
        pair_id, seq_5p, seq_3p = row['pair_id'], row['seq_5p'], row['seq_3p']
        merged_assembly, seq_5p_ostart, seq_5p_oend, seq_3p_ostart, seq_3p_oend  = find_overlap(
            seq_5p, seq_3p, min_overlap_score, min_overlap_perc_identity)
        if merged_assembly is not None:
            loc = 'seq_5p:' + str(seq_5p_ostart) + '-' + str(seq_5p_oend) + ';' + 'seq_3p:' + str(seq_3p_ostart) + '-' + str(seq_3p_oend)
            sequences[pair_id].append((loc, len(merged_assembly), merged_assembly))

    return sequences

def find_overlap(seq_5p, seq_3p, min_overlap_score, min_overlap_perc_identity, verbose=False):
    align_info = get_best_sliding_alignment(seq_5p, seq_3p)
    best_score, best_score_mismatches, best_r_start, best_r_end, best_q_start, best_q_end = align_info
    align_length = best_r_end - best_r_start

    merged_assembly = None

    if best_score >= min_overlap_score and (1-best_score_mismatches / float(align_length)) > min_overlap_perc_identity:
        if verbose:
            click.echo("They merged!")
            alignments = pairwise2.align.localms(seq_5p, seq_3p, 1, -1, -5, -5)
            print(format_alignment(*alignments[0]))

        else:
            merged_assembly = merge_overlapping_sequences(seq_5p[:best_q_start], seq_5p[best_q_start:best_q_end],
                                                          seq_3p[best_r_start:best_r_end], seq_3p[best_r_end:])


    return merged_assembly, best_q_start, best_q_end, best_r_start, best_r_end


def get_best_sliding_alignment(query_read, ref_read):

    #print("Query Read:", query_read)
    #print("Reference Read:", ref_read)

    max_len = min(len(query_read), len(ref_read))
    max_len_query, max_len_ref = False, False

    if max_len == len(query_read):
        max_len_query = True
    else:
        max_len_ref = True

    best_score = -sys.maxsize
    best_score_mismatches = 0
    best_r_start = 0
    best_r_end = 0
    best_q_start = 0
    best_q_end = 0

    q_start = len(query_read) - 1
    q_end = len(query_read)
    r_start = 0
    r_end = 1

    reached_max = False

    while q_start != q_end:
        score = 0
        mismatches = 0
        for i in range(q_end - q_start):

            if query_read[q_start:q_end][i] == ref_read[r_start:r_end][i]:
                score += 1

            else:
                score -= 1
                mismatches += 1

        if score > best_score:
            best_score = score
            best_score_mismatches = mismatches
            best_r_start = r_start
            best_r_end = r_end
            best_q_start = q_start
            best_q_end = q_end

        if q_end - q_start == max_len:
            reached_max = True

        if not reached_max:
            q_start -= 1
            r_end += 1

        elif reached_max and max_len_query:
            if r_end == len(ref_read):
                r_start += 1
                q_end -= 1
            else:
                r_start += 1
                r_end += 1

        elif reached_max and max_len_ref:
            if q_start == 0:
                r_start += 1
                q_end -= 1
            else:
                q_start -= 1
                q_end -= 1

    return best_score, best_score_mismatches, best_r_start, best_r_end, best_q_start, best_q_end

def handle_empty_pairsfile(pairs, output_file):
    if pairs.shape[0] == 0:
        outfile = pd.DataFrame(columns=['pair_id', 'method', 'loc', 'inferred_seq_length', 'inferred_seq'])

        if not output_file:
            output_file = 'mgefinder.inferseq_overlap.tsv'

        outfile.to_csv(output_file, sep='\t', index=False)
        click.echo("Empty pairs file, exiting...")
        sys.exit()

def merge_overlapping_sequences(start_seq, overlap_seq1, overlap_seq2, end_seq):

    merged_length = len(start_seq) + len(overlap_seq1) + len(end_seq)
    merged_seq = start_seq

    for i in range(len(overlap_seq1)):
        if len(merged_seq) > merged_length/2:
            merged_seq += overlap_seq2[i]
        else:
            merged_seq += overlap_seq1[i]
    merged_seq += end_seq

    return(merged_seq)


if __name__ == '__main__':
    inferseq_overlap()