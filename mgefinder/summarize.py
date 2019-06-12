import sys
import pandas as pd
from collections import defaultdict
import pygogo as gogo
import numpy as np
from mgefinder import misc

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

def _summarize(clusterseq, genotypes, output_prefix):

    summary_clusters_path = output_prefix + '.clusters.tsv'
    summary_groups_path = output_prefix + '.groups.tsv'

    logger.info('Reading clusterseq...')
    clusterseq = pd.read_table(clusterseq)

    logger.info('Reading genotypes...')
    genotypes = pd.read_table(genotypes)


    # GET RID OF THIS WHEN YOU FIX THE CLUSTERING
    clusterseq = clusterseq[-clusterseq.cluster.str.startswith('UClust')]
    genotypes = genotypes[-genotypes.cluster.str.startswith('UClust')]


    logger.info("Summarizing clusters...")
    cluster_summary = summarize_clusters(clusterseq, genotypes)

    logger.info("Summarizing groups...")
    group_summary = summarize_groups(clusterseq, genotypes)

    logger.info("Writing results to files %s and %s" % (summary_clusters_path, summary_groups_path))
    cluster_summary.to_csv(summary_clusters_path, sep='\t', index=False)
    group_summary.to_csv(summary_groups_path, sep='\t', index=False)


def summarize_clusters(clusterseq, genotypes):

    clusters2groups = get_clusters2groups(clusterseq.cluster, clusterseq.group)
    num_unique_sites = get_num_unique_sites(genotypes.cluster,
                                            zip(genotypes.contig, genotypes.pos_5p, genotypes.pos_3p),
                                            genotypes.conf)
    conf_counts = get_conf_counts(genotypes.cluster, genotypes.conf)
    seq_lengths = get_seq_lengths(clusterseq.cluster, clusterseq.inferred_seq)
    repr_seqs = get_repr_seq(clusterseq.cluster, clusterseq.seqid, clusterseq.inferred_seq)
    sample_stats = get_sample_stats(genotypes.cluster, genotypes['sample'])

    df1 = pd.merge(clusters2groups, num_unique_sites, how='inner', on=['cluster']).drop_duplicates()
    df2 = pd.merge(df1, seq_lengths, how='inner', on=['cluster']).drop_duplicates()
    df3 = pd.merge(df2, sample_stats, how='inner', on=['cluster']).drop_duplicates()
    df4 = pd.merge(df3, conf_counts, how='inner', on=['cluster']).drop_duplicates()
    final = pd.merge(df4, repr_seqs, how='inner', on=['cluster']).drop_duplicates().sort_values(
        ['num_unique_sites_all'], ascending=False
    )

    return final


def summarize_groups(clusterseq, genotypes):

    num_unique_sites = get_num_unique_sites(genotypes.group, zip(genotypes.contig, genotypes.pos_5p, genotypes.pos_3p),
                                            genotypes.conf)
    conf_counts = get_conf_counts(genotypes.group, genotypes.conf)
    seq_lengths = get_seq_lengths(clusterseq.group, clusterseq.inferred_seq)
    repr_seqs = get_repr_seq(clusterseq.group, clusterseq.seqid, clusterseq.inferred_seq)
    sample_stats = get_sample_stats(genotypes.group, genotypes['sample'])
    repr_cluster2seq = get_repr_cluster2seq(clusterseq.cluster, clusterseq.inferred_seq)
    repr_seqs = (pd.merge(repr_cluster2seq, repr_seqs, how='inner', on='repr_seq')[[
        'cluster', 'repr_cluster', 'repr_seqid', 'repr_seq_length', 'repr_seq'
    ]])

    df1 = pd.merge(num_unique_sites, seq_lengths, how='inner', on=['cluster']).drop_duplicates()
    df2 = pd.merge(df1, sample_stats, how='inner', on=['cluster']).drop_duplicates()
    df3 = pd.merge(df2, conf_counts, how='inner', on=['cluster']).drop_duplicates()
    final = pd.merge(df3, repr_seqs, how='inner', on=['cluster']).drop_duplicates().sort_values(
        ['num_unique_sites_all'], ascending=False
    )

    final.columns = ['group'] + list(final.columns[1:])

    return final


def get_clusters2groups(clusters, groups, header=['cluster', 'group']):

    clusters2groups = list(zip(clusters, groups))
    clusters2groups = pd.DataFrame(clusters2groups, columns = header)

    return clusters2groups


def get_num_unique_sites(clusters, sites, conf, header=['cluster', 'num_unique_sites_all', 'num_unique_sites_unambig']):

    unique_sites_all = defaultdict(set)
    unique_sites_unambig = defaultdict(set)

    for cluster, site, conf in zip(clusters, sites, conf):

        ambig = False
        if conf == 'A':
            ambig = True

        for clust in cluster.split(';'):

            unique_sites_all[clust].add(site)

            if not ambig:
                unique_sites_unambig[clust].add(site)


    unique_sites_all = [(cluster, len(unique_sites_all[cluster])) for cluster in unique_sites_all]
    unique_sites_unambig = [(cluster, len(unique_sites_unambig[cluster])) for cluster in unique_sites_unambig]

    unique_sites_all = pd.DataFrame(unique_sites_all, columns=header[:-1])
    unique_sites_unambig = pd.DataFrame(unique_sites_unambig, columns=[header[0]] + [header[-1]])

    out = pd.merge(unique_sites_all, unique_sites_unambig, how='left', on='cluster').fillna(0)

    return out


def get_conf_counts(clusters, confs, header=['cluster', 'IAwFC', 'IAwHC', 'IO', 'IAwoC', 'IDB',
                                             'ArSC', 'ArMS', 'ArML', 'A']):

    conf_counts = defaultdict(lambda: {lab: 0 for lab in header[1:]})

    for cluster, conf in zip(clusters, confs):

        for clust in cluster.split(';'):
            conf_counts[clust][conf] += 1

    conf_counts = [tuple([cluster] + [conf_counts[cluster][lab] for lab in header[1:]]) for cluster in conf_counts]

    conf_counts = pd.DataFrame(conf_counts, columns=header)

    return conf_counts


def get_seq_lengths(clusters, seqs, header=['cluster', 'num_unique_seqs', 'mean_length', 'min_length', 'max_length']):

    seq_lengths1 = defaultdict(set)

    for cluster, seq in zip(clusters, seqs):
        seq_lengths1[cluster].add(seq)

    seq_lengths2 = defaultdict(lambda: {'unique_seqs': set(), 'num_unique_seqs': 0,
                                        'mean_length': 0, 'min_length': 0, 'max_length': 0})
    for cluster in seq_lengths1:
        for seq in seq_lengths1[cluster]:
            seq_lengths2[cluster]['unique_seqs'].add(tuple(sorted([seq, misc.revcomp(seq)])))

    for cluster in seq_lengths2:

        all_seq_lengths = list(map(lambda x: len(x[0]), seq_lengths2[cluster]['unique_seqs']))
        seq_lengths2[cluster]['num_unique_seqs'] = len(all_seq_lengths)
        seq_lengths2[cluster]['mean_length'] = np.mean(all_seq_lengths)
        seq_lengths2[cluster]['min_length'] = np.min(all_seq_lengths)
        seq_lengths2[cluster]['max_length'] = np.max(all_seq_lengths)

    seq_lengths2 = [tuple([cluster] + [seq_lengths2[cluster][lab] for lab in header[1:]]) for cluster in seq_lengths2]
    seq_lengths2 = pd.DataFrame(seq_lengths2, columns=header)

    return seq_lengths2


def get_repr_seq(clusters, seqids, seqs, header=['cluster', 'repr_seqid', 'repr_seq_length', 'repr_seq']):

    seq_counts = defaultdict(lambda: defaultdict(int))

    for cluster, seqid, seq in zip(clusters, seqids, seqs):
        seq_counts[cluster][(seqid, seq)] += 1


    repr_seq = dict()
    for clust in seq_counts:

        max_seq = max(seq_counts[clust].items(), key=lambda x: x[1])

        repr_seq[clust] = max_seq[0]

    repr_seq = [(clust, repr_seq[clust][0], len(repr_seq[clust][1]), repr_seq[clust][1]) for clust in repr_seq]

    repr_seq = pd.DataFrame(repr_seq, columns=header)

    return repr_seq


def get_sample_stats(clusters, samples, header=['cluster', 'num_samples', 'mean_copy', 'max_copy']):

    sample_counts = defaultdict(lambda: defaultdict(int))

    for cluster, samp in zip(clusters, samples):

        for clust in cluster.split(';'):
            sample_counts[clust][samp] += 1


    sample_stats = defaultdict(dict)
    for cluster in sample_counts:

        counts = [x[1] for x in sample_counts[cluster].items()]
        sample_stats[cluster]['num_samples'] = len(sample_counts[cluster])
        sample_stats[cluster]['mean_copy'] = np.mean(counts)
        sample_stats[cluster]['max_copy'] = np.max(counts)

    sample_stats = [(cluster, sample_stats[cluster]['num_samples'],
                     sample_stats[cluster]['mean_copy'], sample_stats[cluster]['max_copy'])
                    for cluster in sample_stats]

    sample_stats = pd.DataFrame(sample_stats, columns=header)

    return sample_stats

def get_repr_cluster2seq(clusters, seqs, header=['repr_cluster', 'repr_seq']):

    repr_cluster2seq = set()

    for cluster, seq in zip(clusters, seqs):

        repr_cluster2seq.add((cluster, seq))

    repr_cluster2seq = pd.DataFrame(list(repr_cluster2seq), columns=header)

    return repr_cluster2seq