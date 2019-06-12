import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from mgefinder import fastatools, cdhittools
from Bio import SeqIO
from snakemake import shell
from collections import defaultdict
import click
from os.path import join, isfile
from os import makedirs, rename
from shutil import rmtree
from mgefinder.bowtie2tools import index_genome


def _makedatabase(inferseqfiles, minimum_size, maximum_size, threads, memory, force, output_dir, prefix):

    click.echo("Parsing inferseq files")
    if len(inferseqfiles) == 1 and is_path_list(inferseqfiles[0]):
        inferseqfiles = [l.strip() for l in open(inferseqfiles[0], 'r')]

    click.echo("Combining the inferseq files...")
    inferseq = combine_inferseq_files(inferseqfiles, minimum_size, maximum_size)

    database_maker = DatabaseMaker(inferseq, threads, memory, output_dir)

    try:
        makedirs(output_dir)

    except FileExistsError:

        if force:
            click.echo('Deleting old database directory...')
            rmtree(output_dir)
            makedirs(output_dir)
        else:
            click.echo('Output directory already exists. Use --force to overwrite directory.')
            sys.exit()

    if inferseq.shape[0] == 0:
        click.echo("No termini found in the input file...")

        clustered_seqs = database_maker.get_header_dataframe()

    else:

        cluster_fna = database_maker.run_cluster_seqs()
        outfile = join(output_dir, prefix+'.fna')
        rename(cluster_fna, outfile)
        database_maker.remove_int_files()
        click.echo("Indexing database for use with bowtie2")
        index_genome(outfile)

def is_path_list(f):
    with open(f) as infile:
        for l in infile:

            if len(l.strip().split()) > 1:
                return False

            if not isfile(l.strip()):
                return False

    return True


class DatabaseMaker:

    inferseq = None
    outfile = None

    threads = None
    memory = None

    seq_fasta_path = None
    cluster_100p_path = None
    cluster_mappings = None
    final_cluster_mappings = None

    def __init__(self, inferseq, threads, memory, output_dir):

        self.inferseq = inferseq
        self.pair_id_list = list(self.inferseq['pair_id'])

        self.inferseq.fillna('N', inplace=True)
        self.outdir = output_dir

        self.threads = threads
        self.memory = memory

        self.seq_fasta_path = output_dir + '/mgefinder.inferseq.fna'
        self.cluster_100p_path = output_dir + '/mgefinder.inferseq.100p.fna'
        self.cluster_perc_ident_path = output_dir + '/mgefinder.inferseq.perc_ident.fna'


    def run_cluster_seqs(self):

        self.create_inferseq_fasta()
        self.cluster_100p()
        self.cluster_perc_identity()

        return self.cluster_perc_ident_path


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
                                               alignment_coverage=0.99, perc_identity=0.99)


    def map_seqs_to_clusters(self):

        starting_ids = [rec.id for rec in SeqIO.parse(self.seq_fasta_path, 'fasta')]
        cdhit_100p_clusters = cdhittools.parse_cdhit_output_100p(self.seq_fasta_path, self.cluster_100p_path)
        cdhit_perc_identity_clusters = cdhittools.parse_cdhit_output_perc_identity(self.cluster_perc_ident_path)


        final_cluster_mappings = []

        unclustered_count = 0
        for ident in starting_ids:
            if ident in cdhit_100p_clusters.keys() and cdhit_100p_clusters[ident] in cdhit_perc_identity_clusters:
                final_cluster_mappings.append((ident, 'PCluster' + cdhit_perc_identity_clusters[cdhit_100p_clusters[ident]]))
            else:
                final_cluster_mappings.append((ident, 'UCluster' + str(unclustered_count)))
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
                final_clusters[clust] = "FCluster" + str(i+1)

        self.final_cluster_mappings = final_clusters


    def make_final_dataframe(self):

        out_cluster = pd.DataFrame(self.inferseq)
        out_cluster['cluster_perc_ident'] = [clust[1] for clust in self.cluster_mappings]
        out_cluster['cluster_final'] = [self.final_cluster_mappings[clust[1]] for clust in self.cluster_mappings]
        out_cluster = out_cluster[self.get_header_list()]

        return out_cluster

    def get_header_list(self):
        header = ['pair_id', 'method', 'loc', 'inferred_seq_length', 'cluster_perc_ident', 'cluster_final',
                  'inferred_seq']

        return header


    def get_header_dataframe(self):
        terminus_pairs = pd.DataFrame(columns=self.get_header_list())

        return terminus_pairs

    def remove_int_files(self):

        shell("rm %s" % ' '.join([
            self.seq_fasta_path+'*', self.cluster_100p_path + '*', self.cluster_perc_ident_path + '*'
        ]))




def combine_inferseq_files(inferseq_files, minimum_size, maximum_size):
    click.echo('Loading file {num1}/{num2}: {f}'.format(num1=1, num2=len(inferseq_files), f=inferseq_files[0]))
    inferseq = pd.read_csv(inferseq_files[0], sep='\t')
    for i, f in enumerate(inferseq_files[1:]):
        click.echo('Loading file {num1}/{num2}: {f}'.format(num1=i+2, num2=len(inferseq_files), f=f))
        inferseq = inferseq.append(pd.read_csv(f, sep='\t').
                                   query('inferred_seq_length >= @minimum_size').
                                   query('inferred_seq_length <= @maximum_size'))
    return inferseq
