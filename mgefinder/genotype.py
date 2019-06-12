import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from collections import defaultdict
import click
from os.path import isfile


def _genotype(clusterseq, pairfiles, filter_clusters_inferred_assembly, output_file):


    click.echo("Loading clusterseq file...")
    clusterseq = pd.read_table(clusterseq)[['sample', 'pair_id', 'method', 'seqid', 'cluster', 'group']]

    click.echo("Parsing pair files")
    if len(pairfiles) == 1 and is_path_list(pairfiles[0]):
        pairfiles = [l.strip() for l in open(pairfiles[0], 'r')]

    click.echo("Loading pair files...")
    pairs = combine_pair_files(pairfiles)


    genotyper = Genotyper(clusterseq, pairs, filter_clusters_inferred_assembly)


    if pairs.shape[0] == 0:
        click.echo("No termini found in the input file...")

        genotypes = genotyper.get_header_dataframe()

        if output_file:
            genotypes.to_csv(output_file, sep='\t', index=False)


    else:

        genotypes = genotyper.genotype()

        if output_file:
            click.echo("Saving results to file %s" % output_file)
            genotypes.to_csv(output_file, sep='\t', index=False)

        return genotypes


def combine_pair_files(pairfiles):

    keep_header = ['sample', 'pair_id', 'contig', 'pos_5p', 'pos_3p']

    click.echo('Loading file {num1}/{num2}: {f}'.format(num1=1, num2=len(pairfiles), f=pairfiles[0]))

    f0 = pairfiles[0]
    pairs = pd.read_table(f0)[keep_header]


    for i, f in enumerate(pairfiles[1:]):
        click.echo('Loading file {num1}/{num2}: {f}'.format(num1=i + 2, num2=len(pairfiles), f=f))

        df = pd.read_table(f)[keep_header]
        pairs = pairs.append(df)

    return pairs


def is_path_list(f):
    with open(f) as infile:
        for l in infile:

            if len(l.strip().split()) > 1:
                return False

            if not isfile(l.strip()):
                return False

    return True

class Genotyper:
    clusterseq = None
    pairs = None
    filter_clusters_inferred_assembly = None

    def __init__(self, clusterseq, pairs, filter_clusters_inferred_assembly=True):

        self.clusterseq = clusterseq
        self.pairs = pairs
        self.filter_clusters_inferred_assembly = filter_clusters_inferred_assembly

        self.pairs['pair_id'] = list(map(int, list(self.pairs['pair_id'])))
        self.pairs['pos_5p'] = list(map(int, list(self.pairs['pos_5p'])))
        self.pairs['pos_3p'] = list(map(int, list(self.pairs['pos_3p'])))
        self.clusterseq['pair_id'] = list(map(int, list(self.clusterseq['pair_id'])))

    def genotype(self):

        if self.filter_clusters_inferred_assembly:
            click.echo("Filtering out clusters that are never inferred from an assembly...")
            self.clusterseq = self.apply_filter_clusters_inferred_assembly(self.clusterseq)

        combined = pd.merge(self.pairs, self.clusterseq, how='left', on=['sample', 'pair_id'])

        no_inference = combined[combined.method.isnull()]
        has_inference = combined[combined.method.notnull()]

        self.report_null_genotypes(combined, no_inference, has_inference)
        genotypes = self.assign_genotypes_heuristic(has_inference)
        genotypes = self.resolve_ambiguous_genotypes(genotypes)

        return genotypes

    def apply_filter_clusters_inferred_assembly(self, clusterseq):

        method_counts = defaultdict(lambda: defaultdict(int))
        for cluster, method in zip(clusterseq.cluster, clusterseq.method):
            method_counts[cluster][method] += 1

        keep_clusters = set([clust for clust in method_counts
                             if method_counts[clust]['inferred_assembly_with_full_context'] > 0 or
                             method_counts[clust]['inferred_assembly_with_half_context'] > 0 or
                             method_counts[clust]['inferred_assembly_without_context'] > 0])

        exclude_clusters = set([clust for clust in set(clusterseq.cluster) if clust not in keep_clusters])

        click.echo(
            "Excluding %d clusters that were only inferred from the reference genome..." % len(exclude_clusters))

        clusterseq = clusterseq[[clust in keep_clusters for clust in clusterseq.cluster]]

        return clusterseq

    def report_null_genotypes(self, combined, no_inference, has_inference):

        total_count = combined[['sample', 'pair_id']].drop_duplicates().shape[0]
        total_no_inference = no_inference[['sample', 'pair_id']].drop_duplicates().shape[0]
        total_has_inference = has_inference[['sample', 'pair_id']].drop_duplicates().shape[0]

        click.echo(
            "Out of %d candidate insertions, %d had some inferred identity, while %d had no inferred identity." %
            (total_count, total_has_inference, total_no_inference)
        )

    def assign_genotypes_heuristic(self, data):

        click.echo("Assigning initial genotypes...")

        # IAwFC
        IAwFC = data[data.method == 'inferred_assembly_with_full_context']
        IAwFC['conf'] = 'IAwFC'
        IAwFC = self.combine_unresolved_seqs(IAwFC.drop(columns='method'))

        genotyped = IAwFC[['sample', 'pair_id']].drop_duplicates()
        genotyped = genotyped.convert_objects(convert_numeric=True)

        # IAwHC
        IAwHC = pd.merge(data, genotyped, on=['sample', 'pair_id'], how='outer', indicator=True)
        IAwHC = IAwHC[IAwHC._merge == 'left_only']
        IAwHC = IAwHC.drop(columns='_merge')
        IAwHC = IAwHC[IAwHC.method == 'inferred_assembly_with_half_context']
        IAwHC['conf'] = 'IAwHC'
        IAwHC = self.combine_unresolved_seqs(IAwHC.drop(columns='method'))

        genotyped = pd.concat([genotyped, IAwHC[['sample', 'pair_id']]]).drop_duplicates()
        genotyped = genotyped.convert_objects(convert_numeric=True)

        # IO
        IO = pd.merge(data, genotyped, on=['sample', 'pair_id'], how='outer', indicator=True)
        IO = IO[IO._merge == 'left_only']
        IO = IO.drop(columns='_merge')
        IO = IO[IO.method == 'inferred_overlap']
        IO['conf'] = 'IO'
        IO = self.combine_unresolved_seqs(IO.drop(columns='method'))

        genotyped = pd.concat([genotyped, IO[['sample', 'pair_id']]]).drop_duplicates()
        genotyped = genotyped.convert_objects(convert_numeric=True)

        # IAwoC
        IAwoC = pd.merge(data, genotyped, on=['sample', 'pair_id'], how='outer', indicator=True)
        IAwoC = IAwoC[IAwoC._merge == 'left_only']
        IAwoC = IAwoC.drop(columns='_merge')
        IAwoC = IAwoC[IAwoC.method == 'inferred_assembly_without_context']
        IAwoC['conf'] = 'IAwoC'
        IAwoC = self.combine_unresolved_seqs(IAwoC.drop(columns='method'))

        genotyped = pd.concat([genotyped, IAwoC[['sample', 'pair_id']]]).drop_duplicates()
        genotyped = genotyped.convert_objects(convert_numeric=True)

        # IDB
        IDB = pd.merge(data, genotyped, on=['sample', 'pair_id'], how='outer', indicator=True)
        IDB = IDB[IDB._merge == 'left_only']
        IDB = IDB.drop(columns='_merge')
        IDB = IDB[(IDB.method == 'inferred_database') | (IDB.method == 'inferred_reference')]
        IDB['conf'] = 'IDB'
        IDB = self.combine_unresolved_seqs(IDB.drop(columns='method'))

        return IAwFC.append(IAwHC).append(IO).append(IAwoC).append(IDB).drop_duplicates().convert_objects(convert_numeric=True)

    def resolve_ambiguous_genotypes(self, data):

        click.echo("Identifying ambiguous genotypes...")

        counts = data.groupby(['sample', 'pair_id']).size().rename('n').reset_index()

        data = pd.merge(data, counts, on=['sample', 'pair_id'], how='inner')
        data['ambiguous'] = data.n > 1
        data = data.drop(columns='n')
        

        ambiguous = data[data.ambiguous == True]
        not_ambiguous = data[data.ambiguous == False]

        if ambiguous.shape[0] == 0:
            return not_ambiguous

        click.echo("Resolving ambiguous genotypes where possible...")

        mobile_clusters = self.calculate_mobile_clusters(data)
        cluster_counts_per_site = self.count_clusters_per_site(data)

        unresolved = ambiguous.copy()

        resolved_all_sample_comparison, unresolved = self.resolve_all_sample_comparison(
            unresolved, cluster_counts_per_site
        )


        resolved_mobile, unresolved = self.resolve_mobile(unresolved, mobile_clusters)

        combined_unresolved = self.combine_unresolved_clusters(unresolved)

        out = pd.concat([not_ambiguous, resolved_all_sample_comparison,
                         resolved_mobile, combined_unresolved], sort=False)[self.get_header_list()]

        return out

    def combine_unresolved_seqs(self, unresolved):

        combined = defaultdict(set)

        for sample, pair_id, contig, pos_5p, pos_3p, group, conf, cluster, seqid in zip(
                unresolved['sample'], unresolved.pair_id, unresolved.contig, unresolved.pos_5p, unresolved.pos_3p,
                unresolved.group, unresolved.conf, unresolved.cluster, unresolved.seqid):

            combined[(sample, pair_id, contig, pos_5p, pos_3p, cluster, group, conf)].add(seqid)

        combined = [list(key) + [';'.join(sorted(list(combined[key])))] for key in combined.keys()]

        combined = pd.DataFrame(
            combined, columns=['sample', 'pair_id', 'contig', 'pos_5p', 'pos_3p', 'cluster', 'group', 'conf', 'seqid']
        )

        combined = combined[self.get_header_list()]

        return combined

    def combine_unresolved_clusters(self, unresolved):

        combined = defaultdict(lambda: [set(), set()])

        for sample, pair_id, contig, pos_5p, pos_3p, group, conf, cluster, seqid in zip(
                unresolved['sample'], unresolved.pair_id, unresolved.contig, unresolved.pos_5p, unresolved.pos_3p,
                unresolved.group, unresolved.conf, unresolved.cluster, unresolved.seqid):

            combined[(sample, pair_id, contig, pos_5p, pos_3p, group, conf)][0].add(seqid)
            combined[(sample, pair_id, contig, pos_5p, pos_3p, group, conf)][1].add(cluster)

        combined = [list(key) + [';'.join(sorted(list(combined[key][0])))] + [';'.join(sorted(list(combined[key][1])))]
                    for key in combined.keys()]

        combined = pd.DataFrame(
            combined, columns=['sample', 'pair_id', 'contig', 'pos_5p', 'pos_3p', 'group', 'conf', 'seqid', 'cluster']
        )

        combined = combined[self.get_header_list()]
        combined['conf'] = 'A'

        return combined


    def calculate_mobile_clusters(self, data):

        mobile_lenient = defaultdict(set)
        for contig, pos_3p, pos_5p, cluster in zip(data.contig, data.pos_3p, data.pos_5p, data.cluster):
            mobile_lenient[cluster].add((contig, pos_3p, pos_5p))

        high_conf = set(['IAwFC', 'IAwHC', 'IO', 'IAwoC'])
        mobile_strict = defaultdict(set)
        for contig, pos_3p, pos_5p, cluster, conf, in zip(data.contig, data.pos_3p, data.pos_5p, data.cluster, data.conf):
            if conf in high_conf:
                mobile_strict[cluster].add((contig, pos_3p, pos_5p))

        mobile_lenient = [(cluster, len(mobile_lenient[cluster])) for cluster in mobile_lenient]
        mobile_lenient = pd.DataFrame(mobile_lenient, columns=['cluster', 'mobile_lenient'])
        mobile_lenient['mobile_lenient'] = [True if i > 2 else False for i in list(mobile_lenient['mobile_lenient'])]

        mobile_strict = [(cluster, len(mobile_strict[cluster])) for cluster in mobile_strict]
        mobile_strict = pd.DataFrame(mobile_strict, columns=['cluster', 'mobile_strict'])
        mobile_strict['mobile_strict'] = [True if i > 2 else False for i in list(mobile_strict['mobile_strict'])]

        mobile = pd.merge(mobile_lenient, mobile_strict, how='left', on='cluster').fillna(0)

        return mobile

    def count_clusters_per_site(self, data):

        cluster_counts = defaultdict(int)
        for contig, pos_5p, pos_3p, cluster in zip(data.contig, data.pos_5p, data.pos_3p, data.cluster):
            cluster_counts[(contig, pos_5p, pos_3p, cluster)] += 1

        cluster_counts = [list(group) + [cluster_counts[group]] for group in cluster_counts]

        cluster_counts = pd.DataFrame(cluster_counts, columns=['contig', 'pos_5p', 'pos_3p', 'cluster', 'n'])

        return cluster_counts

    def resolve_all_sample_comparison(self, unresolved, cluster_counts):

        if unresolved.shape[0] == 0:
            return pd.DataFrame(columns=list(unresolved.columns)), unresolved

        resolved = (pd.merge(unresolved, cluster_counts, how='inner', on=['contig', 'pos_5p', 'pos_3p', 'cluster']).
                    groupby(['sample', 'pair_id']).
                    apply(lambda g: g[g['n'] == g['n'].max()]).
                    reset_index(drop=True).
                    groupby(['sample', 'pair_id']).
                    filter(lambda g: len(g) == 1).
                    assign(ambiguous=lambda x: False).
                    assign(conf=lambda x: 'ArSC').
                    drop(columns='n'))

        unresolved = (pd.merge(unresolved, cluster_counts, how='inner', on=['contig', 'pos_5p', 'pos_3p', 'cluster']).
                      groupby(['sample', 'pair_id']).
                      apply(lambda g: g[g['n'] == g['n'].max()]).
                      reset_index(drop=True).
                      groupby(['sample', 'pair_id']).
                      filter(lambda g: len(g) > 1).
                      drop(columns='n').reset_index(drop=True))

        return resolved, unresolved


    def resolve_mobile(self, unresolved, mobile_clusters):

        if unresolved.shape[0] == 0:
            return pd.DataFrame(columns=list(unresolved.columns)), unresolved

        # Resolve using mobile strict definition
        filtered = pd.merge(unresolved, mobile_clusters, how='inner', on='cluster')
        filtered = filtered[filtered.mobile_strict == True]

        if filtered.shape[0] == 0:
            return pd.DataFrame(columns=list(unresolved.columns)), unresolved

        resolved_mobile_strict = (filtered.
                    groupby(['sample', 'pair_id']).
                    filter(lambda g: len(g) == 1).
                    assign(ambiguous=lambda x: False).
                    assign(conf=lambda x: 'ArMS').
                    drop(columns=['mobile_strict', 'mobile_lenient']))


        # Resolve using mobile lenient definition
        filtered = pd.merge(unresolved, mobile_clusters, how='inner', on='cluster')
        filtered = filtered[filtered.mobile_lenient == False]

        filtered = pd.merge(filtered, resolved_mobile_strict[['sample', 'pair_id']], how='outer',
                            on=['sample', 'pair_id'], indicator=True)
        filtered = filtered[filtered._merge == 'left_only'].drop(columns='_merge')


        if filtered.shape[0] == 0:
            return pd.DataFrame(columns=list(unresolved.columns)), unresolved

        resolved_mobile_lenient = (filtered.
                                   groupby(['sample', 'pair_id']).
                                   filter(lambda g: len(g) == 1).
                                   assign(ambiguous=lambda x: False).
                                   assign(conf=lambda x: 'ArML').
                                   drop(columns=['mobile_strict', 'mobile_lenient']))

        # Now get the unresolved
        all_resolved = pd.concat([resolved_mobile_strict, resolved_mobile_lenient]).drop_duplicates()
        merged = pd.merge(unresolved, all_resolved[['sample', 'pair_id']], on=['sample', 'pair_id'], how='outer', indicator=True)
        merged = merged[merged._merge == 'left_only']
        unresolved = merged.drop(columns='_merge').assign(conf=lambda x: 'A')

        
        return all_resolved, unresolved

    def get_header_list(self):
        header = ['sample', 'pair_id', 'contig', 'pos_3p', 'pos_5p', 'seqid', 'cluster', 'group', 'conf']

        return header

    def get_header_dataframe(self):
        terminus_pairs = pd.DataFrame(columns=self.get_header_list())

        return terminus_pairs