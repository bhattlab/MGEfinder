import warnings
warnings.filterwarnings("ignore")
from mgefinder import sctools
from mgefinder.find import SoftclipParser, SoftclipSite
import pandas as pd
import pysam
import pygogo as gogo

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


def _recall(pairsfile, bamfile, min_alignment_quality, min_alignment_inner_length, large_insertion_cutoff, output_file):
    pairs = pd.read_csv(pairsfile, sep='\t')
    bam = pysam.AlignmentFile(bamfile)

    recaller = Recaller(bam, pairs, min_alignment_quality, min_alignment_inner_length, large_insertion_cutoff)

    recaller.parse_clipped_and_unclipped_read_info()
    recall_out = recaller.make_dataframe()

    if output_file:
        logger.info("Saving results to file %s" % output_file)
        recall_out.to_csv(output_file, sep='\t', index=False)

    return recall_out


class Recaller(SoftclipParser):

    pairs_dataframe = None

    def __init__(self, bam, pairs_dataframe, min_alignment_quality, min_alignment_inner_length, large_insertion_cutoff):
        SoftclipParser.__init__(self, bam, min_alignment_quality=min_alignment_quality,
                                min_alignment_inner_length=min_alignment_inner_length,
                                large_insertion_cutoff=large_insertion_cutoff)

        self.pairs_dataframe = pairs_dataframe
        self.load_pairs()


    def load_pairs(self):
        for index, row in self.pairs_dataframe.iterrows():
            self.softclipped_sites[row['contig']][row['pos_5p']] = SoftclipSite()
            self.softclipped_sites[row['contig']][row['pos_3p']] = SoftclipSite()


    def parse_clipped_and_unclipped_read_info(self):
        if verbose:
            logger.info("Getting clipped and unclipped read information near softclipped sites...")
            pass

        for contig in self.softclipped_sites:

            for pos in self.softclipped_sites[contig]:

                reads_at_site = self.get_reads_at_site(contig, pos)

                softclipped_5p_reads, softclipped_3p_reads = self.get_clipped_read_info_at_site(contig, pos, reads_at_site)

                runthrough_reads, small_insertion_5p_reads, small_insertion_3p_reads, \
                large_insertion_5p_reads, large_insertion_3p_reads, deletion_reads = \
                    self.get_unclipped_read_info_at_site(contig, pos, reads_at_site)

                self.softclipped_sites[contig][pos].add_softclip_5p_reads(softclipped_5p_reads)
                self.softclipped_sites[contig][pos].add_softclip_3p_reads(softclipped_3p_reads)

                self.softclipped_sites[contig][pos].add_runthrough_reads(runthrough_reads)
                self.softclipped_sites[contig][pos].add_small_insertion_5p_reads(small_insertion_5p_reads)
                self.softclipped_sites[contig][pos].add_small_insertion_3p_reads(small_insertion_3p_reads)
                self.softclipped_sites[contig][pos].add_large_insertion_5p_reads(large_insertion_5p_reads)
                self.softclipped_sites[contig][pos].add_large_insertion_3p_reads(large_insertion_3p_reads)
                self.softclipped_sites[contig][pos].add_deletion_reads(deletion_reads)

                upstream_deletion_reads, downstream_deletion_reads = None, None

                if pos - 1 >= 0:
                    upstream_deletion_reads = self.get_unclipped_read_info_at_site(
                        contig, pos - 1, reads_at_site, deletions_only=True
                    )

                if pos + 1 < self.contig_lengths[contig]:
                    downstream_deletion_reads = self.get_unclipped_read_info_at_site(
                        contig, pos + 1, reads_at_site, deletions_only=True
                    )

                if upstream_deletion_reads:
                    self.softclipped_sites[contig][pos].add_upstream_deletion_reads(upstream_deletion_reads)

                if downstream_deletion_reads:
                    self.softclipped_sites[contig][pos].add_downstream_deletion_reads(downstream_deletion_reads)


    def get_clipped_read_info_at_site(self, contig, pos, reads):

        softclipped_5p_reads = set()
        softclipped_3p_reads = set()

        for read in reads:

            if not self.passes_read_filters(read):
                continue

            if sctools.is_left_softclipped_lenient_at_site(read, contig, pos):
                softclipped_3p_reads.add(read)

            if sctools.is_right_softclipped_lenient_at_site(read, contig, pos):
                softclipped_5p_reads.add(read)

        return softclipped_5p_reads, softclipped_3p_reads


    def make_dataframe(self):

        column_names = ['contig', 'pos', 'orient', 'softclip_count_5p', 'softclip_count_3p', 'runthrough_count',
                        'small_insertion_count_5p', 'small_insertion_count_3p',
                        'large_insertion_count_5p', 'large_insertion_count_3p', 'deletion_count',
                        'upstream_deletion_count', 'downstream_deletion_count', 'total_count']

        outdata= dict()
        for contig in self.softclipped_sites:
            sorted_positions = sorted(list(self.softclipped_sites[contig].keys()))

            for pos in sorted_positions:
                site = self.softclipped_sites[contig][pos]

                outdata[len(outdata)] = [
                    contig, pos, '5p', site.get_softclip_5p_count(), site.get_softclip_3p_count(),
                    site.get_runthrough_count(), site.get_small_insertion_5p_count(), site.get_small_insertion_3p_count(),
                    site.get_large_insertion_5p_count(), site.get_large_insertion_3p_count(),
                    site.get_deletion_count(), site.get_upstream_deletion_count(), site.get_downstream_deletion_count(),
                    site.get_total_count()
                ]

                outdata[len(outdata)] = [
                    contig, pos, '3p', site.get_softclip_5p_count(), site.get_softclip_3p_count(),
                    site.get_runthrough_count(), site.get_small_insertion_5p_count(), site.get_small_insertion_3p_count(),
                    site.get_large_insertion_5p_count(), site.get_large_insertion_3p_count(),
                    site.get_deletion_count(), site.get_upstream_deletion_count(), site.get_downstream_deletion_count(),
                    site.get_total_count()
                ]

        out_df = pd.DataFrame.from_dict(outdata, orient='index', columns=column_names)

        return out_df
