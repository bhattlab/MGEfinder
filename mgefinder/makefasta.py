import sys
import pandas as pd
import pygogo as gogo
from collections import defaultdict
from mgefinder import fastatools
from mgefinder import misc

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def _makefasta(clusterseq, summarize_clusters, output_prefix):

    fasta_all_path = output_prefix + '.all_seqs.fna'
    fasta_repr_path = output_prefix + '.repr_seqs.fna'

    logger.info('Reading clusterseq...')
    clusterseq = pd.read_table(clusterseq)

    logger.info('Reading summarized clusters...')
    summarize_clusters = pd.read_table(summarize_clusters)

    logger.info('Creating fasta for all unique sequences...')
    make_all_unique_fasta(clusterseq, summarize_clusters, fasta_all_path)

    logger.info('Creating fasta for representative sequences...')
    make_repr_cluster_fasta(summarize_clusters, fasta_repr_path)


def make_repr_cluster_fasta(summarize_clusters, fasta_repr_path):

    sequences = []
    names = []
    for cluster, seqid, group, repr_seq in zip(
            summarize_clusters.cluster, summarize_clusters.repr_seqid, summarize_clusters.group,
            summarize_clusters.repr_seq
    ):
        name = '_'.join([cluster, group, seqid])
        names.append(name)
        sequences.append(repr_seq)

    fastatools.write_sequences_to_fasta(sequences, fasta_repr_path, names)


def make_all_unique_fasta(clusterseq, summarize_clusters, fasta_all_path):

    keep_clusters = set(list(summarize_clusters.cluster))
    unique_seqs = defaultdict(set)
    for cluster, group, seqid, seq in zip(
            clusterseq.cluster, clusterseq.group, clusterseq.seqid, clusterseq.inferred_seq):

        if cluster in keep_clusters:
            unique_seqs['_'.join([cluster, group, seqid])].add(seq)

    names = list(unique_seqs.keys())
    sequences = [list(unique_seqs[name])[0] for name in names]

    fastatools.write_sequences_to_fasta(sequences, fasta_all_path, names)
