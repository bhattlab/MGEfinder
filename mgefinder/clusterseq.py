import sys
import warnings
warnings.filterwarnings("ignore")
import pygogo as gogo
import pandas as pd
from mgefinder import fastatools, cdhittools
from Bio import SeqIO
from snakemake import shell
import networkx as nx
from collections import defaultdict
import click
from os.path import isfile

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger
from tqdm import tqdm
from os.path import basename


def _clusterseq(inferseqfiles, minimum_size, maximum_size, threads, memory, output_file):

    logger.info("Parsing inferseq files")
    if len(inferseqfiles) == 1 and is_path_list(inferseqfiles[0]):
        inferseqfiles = [l.strip() for l in open(inferseqfiles[0], 'r')]

    logger.info("Combining the inferseq files...")
    inferseq = combine_inferseq_files(inferseqfiles, minimum_size, maximum_size)

    seq_cluster = SequenceClusterer(inferseq, threads, memory, output_file)

    if inferseq.shape[0] == 0:
        logger.info("No flanks found in the input file...")

        clustered_seqs = seq_cluster.get_header_dataframe()

        if output_file:
            clustered_seqs.to_csv(output_file, sep='\t', index=False)


    else:

        seq_cluster.run_cluster_seqs()
        clustered_seqs = seq_cluster.make_final_dataframe()
        seq_cluster.remove_int_files()

        if output_file:
            logger.info("Saving results to file %s" % output_file)
            clustered_seqs.to_csv(output_file, sep='\t', index=False)

        return clustered_seqs

    return clustered_seqs


def combine_inferseq_files(inferseq_files, minimum_size, maximum_size):
    click.echo('Loading file {num1}/{num2}: {f}'.format(num1=1, num2=len(inferseq_files), f=inferseq_files[0]))

    f0 = inferseq_files[0]
    inferseq = pd.read_table(f0)

    for i, f in enumerate(inferseq_files[1:]):
        click.echo('Loading file {num1}/{num2}: {f}'.format(num1=i + 2, num2=len(inferseq_files), f=f))

        df = pd.read_table(f)

        inferseq = inferseq.append(df)

    inferseq = inferseq[(inferseq.inferred_seq_length >= minimum_size) & (inferseq.inferred_seq_length <= maximum_size)]
    
    return inferseq


def is_path_list(f):
    with open(f) as infile:
        for l in infile:

            if len(l.strip().split()) > 1:
                return False

            if not isfile(l.strip()):
                return False

    return True

class SequenceClusterer:

    inferseq = None
    outfile = None

    threads = None
    memory = None

    seq_fasta_path = None
    cluster_100p_path = None
    cluster_mappings = None
    final_cluster_mappings = None

    def __init__(self, inferseq, threads, memory, output_file):

        self.inferseq = inferseq
        self.pair_id_list = list(zip(self.inferseq['sample'], self.inferseq['pair_id']))

        self.inferseq.fillna('N', inplace=True)
        self.outfile = output_file

        self.threads = threads
        self.memory = memory

        self.seq_fasta_path = output_file + '.clusterseq.fna'
        self.cluster_100p_path = output_file + '.clusterseq.100p.fna'
        self.cluster_perc_ident_path = output_file + '.clusterseq.perc_ident.fna'


    def run_cluster_seqs(self):

        self.create_inferseq_fasta()
        self.cluster_100p()
        self.cluster_perc_identity()
        self.cluster_mappings = self.map_seqs_to_clusters()
        self.cluster_shared_pairs()


    def make_cluster_mappings(self):
        self.cluster_mappings = self.map_seqs_to_clusters()
        self.cluster_shared_pairs()


    def create_inferseq_fasta(self):

        seqs = list(self.inferseq['inferred_seq'])

        fastatools.write_sequences_to_fasta(seqs, self.seq_fasta_path)


    def cluster_100p(self):

        click.echo("Clustering at 100% identity...")

        cdhittools.cluster_100p(self.seq_fasta_path, self.cluster_100p_path)


    def cluster_perc_identity(self):

        cdhittools.cluster_reciprocal_identity(self.cluster_100p_path, self.cluster_perc_ident_path,
                                               threads=self.threads, memory=self.memory,
                                               alignment_coverage=0.85, perc_identity=0.90)


    def map_seqs_to_clusters(self):

        starting_ids = [rec.id for rec in SeqIO.parse(self.seq_fasta_path, 'fasta')]
        cdhit_100p_clusters = cdhittools.parse_cdhit_output_100p(self.seq_fasta_path, self.cluster_100p_path)
        cdhit_perc_identity_clusters = cdhittools.parse_cdhit_output_perc_identity(self.cluster_perc_ident_path)


        final_cluster_mappings = []

        unclustered_count = 0
        for ident in starting_ids:

            if ident in cdhit_100p_clusters.keys() and cdhit_100p_clusters[ident] in cdhit_perc_identity_clusters:
                final_cluster_mappings.append((ident, 'c' + cdhit_perc_identity_clusters[cdhit_100p_clusters[ident]],
                                               's' + str(cdhit_100p_clusters[ident])))
            else:
                final_cluster_mappings.append((ident, 'uc' + str(unclustered_count),
                                               's' + str(cdhit_100p_clusters[ident])))
                unclustered_count += 1

        return final_cluster_mappings


    def cluster_shared_pairs(self):

        G = nx.Graph()

        pair_clusters = defaultdict(set)
        for pair_index in self.cluster_mappings:

            pair = self.pair_id_list[int(pair_index[0])]
            pair_clusters[pair].add(pair_index[1])



        for pair in pair_clusters.keys():

            nodes = list(pair_clusters[pair])

            for node in nodes:
                G.add_node(node)

            if len(nodes) > 1:

                for i in range(len(nodes)-1):
                    n1, n2 = nodes[i], nodes[i+1]
                    G.add_edge(n1, n2)


        final_clusters = dict()
        comps = list(nx.connected_components(G))
        for i in range(len(comps)):
            comp = comps[i]

            for clust in comp:
                final_clusters[clust] = "g" + str(i+1)

        self.final_cluster_mappings = final_clusters


    def make_final_dataframe(self):

        out_cluster = pd.DataFrame(self.inferseq)
        out_cluster['seqid'] = [clust[2] for clust in self.cluster_mappings]
        out_cluster['cluster'] = [clust[1] for clust in self.cluster_mappings]
        out_cluster['group'] = [self.final_cluster_mappings[clust[1]] for clust in self.cluster_mappings]
        out_cluster = out_cluster[self.get_header_list()]

        return out_cluster


    def get_header_list(self):
        header = ['sample', 'pair_id', 'method', 'loc', 'inferred_seq_length', 'seqid', 'cluster', 'group',
                  'inferred_seq']

        return header


    def get_header_dataframe(self):
        terminus_pairs = pd.DataFrame(columns=self.get_header_list())

        return terminus_pairs

    def remove_int_files(self):

        shell("rm %s" % ' '.join([
            self.seq_fasta_path+'*', self.cluster_100p_path + '*', self.cluster_perc_ident_path + '*'
        ]))