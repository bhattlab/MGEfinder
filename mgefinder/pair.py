import warnings
warnings.filterwarnings("ignore")
import click
import pandas as pd
import numpy as np
from snakemake import shell
from random import randint
from mgefinder import fastatools, embosstools, pysamtools, sctools
from os.path import basename, join, dirname
from Bio import SeqIO
import pysam
from collections import defaultdict


def _pair(findfile, bamfile, genome, max_direct_repeat_length, min_alignment_quality,
                min_alignment_inner_length, max_junction_spanning_prop=0.15, large_insertion_cutoff=30,
                output_file=None):
    tmp_output_prefix = '.'.join(basename(output_file).split('.')[:-1])

    flanks = pd.read_csv(findfile, sep='\t')


    bam = pysam.AlignmentFile(bamfile, 'rb')

    flank_pairer = FlankPairer(flanks, bam, genome, max_direct_repeat_length, min_alignment_quality,
                               min_alignment_inner_length, max_junction_spanning_prop, large_insertion_cutoff,
                               tmp_dir=dirname(output_file), tmp_output_prefix=tmp_output_prefix)

    if flanks.shape[0] == 0:
        click.echo("No flanks found in the input file...")

        flank_pairs = flank_pairer.get_header_dataframe()
        flank_pairs.insert(0, 'sample', None)

        if output_file:
            flank_pairs.to_csv(output_file, sep='\t', index=False)
        return flank_pairs

    else:
        flank_pairs = flank_pairer.run_pair_flanks()
        sample_id = list(flanks['sample'])[0]
        flank_pairs.insert(0, 'sample', sample_id)

        if output_file:
            click.echo("Saving results to file %s" % output_file)
            flank_pairs.to_csv(output_file, sep='\t', index=False)

        return flank_pairs


class FlankPairer:

    flanks = None
    bam = None
    genome = None
    max_direct_repeat_length = None
    min_alignment_quality = None
    min_alignment_inner_length = None
    max_junction_spanning_prop = None
    truncated_flank_length = None
    ir_distance_from_end = None

    insertion_spanning_length = None
    large_insertion_cutoff = None

    tmp_dir = None
    tmp_output_prefix = None

    flank_pairs = None


    def __init__(self, flanks, bam, genome, max_direct_repeat_length,
                 min_alignment_quality, min_alignment_inner_length, max_junction_spanning_prop,
                 large_insertion_cutoff,
                 truncated_flank_length=40, ir_distance_from_end=15, insertion_spanning_length=10,
                 tmp_dir='/tmp', tmp_output_prefix='mgefinder'):
        self.flanks = flanks
        self.bam = bam
        self.genome = genome
        self.max_direct_repeat_length = max_direct_repeat_length
        self.min_alignment_quality = min_alignment_quality
        self.min_alignment_inner_length = min_alignment_inner_length
        self.max_junction_spanning_prop = max_junction_spanning_prop
        self.large_insertion_cutoff = large_insertion_cutoff
        self.truncated_flank_length = truncated_flank_length
        self.ir_distance_from_end = ir_distance_from_end
        self.insertion_spanning_length = insertion_spanning_length

        self.tmp_dir = tmp_dir
        self.tmp_output_prefix = tmp_output_prefix


    def run_pair_flanks(self):
        click.echo("Finding all flank pairs within %d bases of each other ..." % self.max_direct_repeat_length)
        pairs = self.pair_all_nearby_flanks(self.flanks)
        click.echo("Finding all inverted repeats at termini in %d candidate pairs..." % pairs.shape[0])
        pairs = self.check_pairs_for_ir(pairs)
        click.echo("Assigning pairs according to existence of inverted repeats, read count difference, and flank length difference...")
        assigned_pairs = self.assign_pairs(pairs)
        click.echo("Filtering out pairs with evidence of reads spanning both clipped junctions...")
        analyzed_pairs = self.count_insertion_spanning_reads(assigned_pairs)

        click.echo("Identified %d flank pairs in total..." % analyzed_pairs.shape[0])
        click.echo("Identified %d flank pairs with inverted repeats..." % analyzed_pairs.query('has_IR==True').shape[0])
        click.echo("Identified %d flank pairs with reads that span the insertion..." % analyzed_pairs.query('spanning_count > 0').shape[0])

        click.echo("Filtering sites with junction-spanning reads...")
        filtered_pairs = self.filter_junction_spanning(analyzed_pairs)
        click.echo("%d flank pairs remain after filtering..." % filtered_pairs.shape[0])


        click.echo("Getting direct repeats and surrounding genomic region...")
        final_pairs = self.get_direct_repeats(filtered_pairs)



        final_pairs['pair_id'] = list(range(1, final_pairs.shape[0]+1))

        return final_pairs


    def get_header_list(self):
        header = ['pair_id', 'contig', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p',
                  'total_count_5p', 'total_count_3p', 'spanning_count', 'has_IR', 'IR_length', 'IR_5p', 'IR_3p',
                  'seq_5p', 'seq_3p',  'direct_repeat_reference', 'direct_repeat_reads_consensus']

        return header

    def get_header_dataframe(self):
        flank_pairs = pd.DataFrame(columns=self.get_header_list())

        return flank_pairs


    def pair_all_nearby_flanks(self, flanks):

        column_names = ['contig', 'index_5p', 'index_3p', 'pos_5p', 'pos_3p', 'softclip_count_5p', 'softclip_count_3p',
                        'total_count_5p', 'total_count_3p']

        column_names += ['seq_5p', 'seq_3p']

        outpairs = dict()

        for index_5p, row in flanks.iterrows():

            if row.orient != '5p':
                continue

            contig, pos_5p = row['contig'], row['pos']
            softclip_count_5p, total_count_5p, seq_5p = row['consensus_softclip_count'], row['total_count'], row['consensus_seq']

            min_pos = pos_5p - self.max_direct_repeat_length - 1

            candidate_pairs = flanks.query('contig == @contig & pos >= @min_pos & pos < @pos_5p & orient == "3p"')

            for index_3p, row2 in candidate_pairs.iterrows():
                pos_3p, softclip_count_3p, total_count_3p, seq_3p = row2['pos'], row2['consensus_softclip_count'], \
                                                                         row2['total_count'], row2['consensus_seq']

                out = [contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, total_count_5p, total_count_3p]

                out += [seq_5p, seq_3p]
                outpairs[len(outpairs)] =  out

        outpairs = pd.DataFrame.from_dict(outpairs, orient='index', columns=column_names)
        return outpairs


    def check_pairs_for_ir(self, pairs):

        has_ir_all = []
        ir_5p_all = []
        ir_3p_all = []

        for index, row in pairs.iterrows():

            contig, index_5p, index_3p, pos_5p, pos_3p, softclip_count_5p, softclip_count_3p, \
            runthrough_count_5p, runthrough_count_3p, seq_5p, seq_3p = row

            trunc_seq_5p = self.truncate_sequence(seq_5p, self.truncated_flank_length, orient='5p')
            trunc_seq_3p = self.truncate_sequence(seq_3p, self.truncated_flank_length, orient='3p')

            combined_seq = str('N'*20).join([trunc_seq_5p, trunc_seq_3p])


            tmp_fasta_path = join(self.tmp_dir, self.tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.fasta')
            fastatools.write_sequences_to_fasta([combined_seq], tmp_fasta_path)

            tmp_einverted_outfile = join(self.tmp_dir, self.tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.out')
            tmp_einverted_outseq = join(self.tmp_dir, self.tmp_output_prefix + '.' + str(randint(0, 1e20)) + '.fa')
            embosstools.run_einverted(tmp_fasta_path, outfile=tmp_einverted_outfile, outseq=tmp_einverted_outseq)

            has_ir = False
            ir_length = 0
            keep_ir1 = None
            keep_ir2 = None
            for ir1, ir2 in embosstools.read_emboss_seq_results(tmp_einverted_outseq):
                if self.pair_has_ir(ir1, ir2, self.ir_distance_from_end, len(combined_seq)):
                    has_ir = True
                    if len(ir1.seq) > ir_length:
                        keep_ir1 = ir1.seq
                        keep_ir2 = ir2.seq

            shell('rm -f %s' % tmp_fasta_path)
            shell('rm -f %s' % tmp_einverted_outfile)
            shell('rm -f %s' % tmp_einverted_outseq)

            has_ir_all.append(has_ir)
            ir_3p_all.append(keep_ir1)
            ir_5p_all.append(keep_ir2)

        pairs['has_IR'] = has_ir_all
        pairs['IR_5p'] = ir_5p_all
        pairs['IR_3p'] = ir_3p_all

        return pairs


    def assign_pairs(self, pairs):
        pairs.loc[:, 'direct_repeat_length'] = pairs.loc[:, 'pos_5p'] - pairs.loc[:, 'pos_3p'] -1

        pairs['IR_length'] = np.array([len(seq) if seq is not None else 0 for seq in list(pairs.loc[:,'IR_5p'])])
        pairs['diffcount'] = abs(pairs['softclip_count_5p'] - pairs['softclip_count_3p'])
        pairs['difflength'] = abs(np.array(list(map(len, pairs['seq_5p']))) - np.array(list(map(len, pairs['seq_3p']))))
        pairs['ignore_pair'] = False
        sorted_pairs = pairs.sort_values(['IR_length', 'diffcount', 'difflength'], ascending=[False, True, True])

        keep_row = []

        for index in sorted_pairs.index:
            row = sorted_pairs.loc[index, :]

            if not row['ignore_pair']:
                keep_row.append(True)
                index_5p, index_3p = row['index_5p'], row['index_3p']
                has_pair_member = sorted_pairs.query("index_5p == @index_5p | index_3p == @index_3p")
                has_pair_member_indices = has_pair_member.query("index_5p != @index_5p | index_3p != @index_3p").index.values
                sorted_pairs.loc[has_pair_member_indices, 'ignore_pair'] = True
            else:
                keep_row.append(False)

        sorted_pairs.loc[:, 'keep_pair'] = np.array(keep_row)

        assigned_pairs = sorted_pairs.query("keep_pair == True").loc[:,
                         self.get_header_list()].sort_values(['contig', 'pos_5p', 'pos_3p'])

        return assigned_pairs


    def pair_has_ir(self, ir1, ir2, ir_distance_from_end, seqlen):
        if self.ir_near_5prime_end(ir1, ir_distance_from_end) and self.ir_near_3prime_end(ir2, ir_distance_from_end, seqlen):
            return True
        return False

    def ir_near_5prime_end(self, ir1, ir_distance_from_end):
        if ir1.ir_pos_5p <= ir_distance_from_end:
            return True
        return False

    def ir_near_3prime_end(self, ir2, ir_distance_from_end, seqlen):
        if ir2.ir_pos_3p >= (seqlen - ir_distance_from_end):
            return True
        return False

    def truncate_sequence(self, seq, truncated_seq_length, orient='5p'):
        truncated_seq = seq
        if len(truncated_seq) > truncated_seq_length:
            if orient == '5p':
                truncated_seq = truncated_seq[:truncated_seq_length]
            elif orient == '3p':
                truncated_seq = truncated_seq[-truncated_seq_length:]
        return truncated_seq

    def get_direct_repeats(self, flank_pairs):

        genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(self.genome, 'fasta')}
        positions = self.get_reference_direct_repeats(flank_pairs, genome_dict)
        positions = self.get_read_direct_repeats(positions, genome_dict)

        flank_pairs = flank_pairs.drop(['direct_repeat_reference', 'direct_repeat_reads_consensus'], axis=1).merge(positions, how='left')

        return flank_pairs


    def get_read_direct_repeats(self, positions, genome_dict, target_region_size=50):

        direct_repeats = []
        for index, row in positions.iterrows():
            contig, start, end = row['contig'], row['pos_3p'], row['pos_5p']

            direct_repeat_center = round((end + start) / 2)
            expanded_start = int(direct_repeat_center - (target_region_size / 2))
            expanded_end = int(direct_repeat_center + (target_region_size / 2))

            add_start_n = 0
            add_end_n = 0
            if expanded_start < 0:
                add_start_n = abs(expanded_start)
                expanded_start = 0
            if expanded_end > len(genome_dict[contig]):
                add_end_n = expanded_end - len(genome_dict[contig])
                expanded_end = len(genome_dict[contig])

            target_region = 'N' * add_start_n + genome_dict[contig][expanded_start:expanded_end] + 'N' * add_end_n
            target_region_positions = range(expanded_start, expanded_end)

            target_region_reads = self.initialize_target_region_reads(target_region, expanded_start, expanded_end)
            for read in self.bam.fetch(contig, expanded_start, expanded_end):
                ref_positions = read.get_reference_positions(full_length=True)
                read_query = read.query_sequence
                read_qualities = read.query_qualities

                i = 0
                for pos in ref_positions:
                    if pos is not None and pos in target_region_reads:
                        target_region_reads[pos][read_query[i]] += read_qualities[i]
                    i += 1

            consensus_target_region = self.get_consensus_target_region(target_region_reads)
            consensus_direct_repeat = consensus_target_region[target_region_positions.index((start+1)):target_region_positions.index(end)]

            direct_repeats.append(consensus_direct_repeat)

        positions['direct_repeat_reads_consensus'] = direct_repeats

        return positions


    def get_consensus_target_region(self, target_region_reads):
        start, end = min(target_region_reads.keys()), max(target_region_reads.keys())+1
        consensus = ''
        for pos in range(start, end):
            best_qual = 0
            best_nuc = ''
            for nuc in target_region_reads[pos]:
                qual = target_region_reads[pos][nuc]
                if qual > best_qual:
                    best_qual = qual
                    best_nuc = nuc
            consensus += best_nuc
        return consensus


    def initialize_target_region_reads(self, target_region, expanded_start, expanded_end):
        target_region_reads = defaultdict(lambda: defaultdict(int))
        i = 0
        for pos in range(expanded_start, expanded_end):
            target_region_reads[pos][target_region[i]] += 1
            i += 1
        return target_region_reads


    def get_reference_direct_repeats(self, flank_pairs, genome_dict, target_region_size=50):
        positions = flank_pairs.loc[:, ['contig', 'pos_5p', 'pos_3p']].drop_duplicates().reset_index(drop=True)

        direct_repeats = []
        for index, row in positions.iterrows():
            contig, start, end = row['contig'], row['pos_3p'], row['pos_5p']

            direct_repeat = genome_dict[contig][(start+1):end]
            direct_repeats.append(''.join(direct_repeat))

        positions['direct_repeat_reference'] = direct_repeats
        return positions


    def count_insertion_spanning_reads(self, assigned_pairs):

        spanning_reads = []
        contig_lengths = pysamtools.get_bam_contig_dict(self.bam)
        for index, row in assigned_pairs.iterrows():
            spanning_read_names = set()
            contig, pos_3p, pos_5p = row['contig'], row['pos_3p'], row['pos_5p']
            reads = self.get_reads_at_site(contig, pos_3p, pos_5p, self.bam, contig_lengths)

            for read in reads:

                if (read.reference_start < pos_3p - self.insertion_spanning_length and
                    read.reference_end > pos_5p + self.insertion_spanning_length + 1 and
                    not self.contains_large_insertion(read, pos_3p, pos_5p)):

                    spanning_read_names.add(read.query_name)

            spanning_reads.append(len(spanning_read_names))

        assigned_pairs['spanning_count'] = spanning_reads
        assigned_pairs = assigned_pairs.loc[:, self.get_header_list()]
        return assigned_pairs

    def get_reads_at_site(self, contig, pos_3p, pos_5p, bam, contig_lengths):

        reads = []

        start = pos_3p - 1
        end = pos_5p + 2

        if start < 0:
            start = 0
        if end > contig_lengths[contig]:
            end = contig_lengths[contig]

        for read in bam.fetch(contig, start, end):
            if self.passes_read_filters(read):
                reads.append(read)

        return reads


    def passes_read_filters(self, read):
        if read.mapping_quality < self.min_alignment_quality:
            return False
        elif not sctools.read_meets_min_alignment_inner_length(read, self.min_alignment_inner_length):
            return False
        else:
            return True


    def filter_junction_spanning(self, pairs):

        pairs['total'] = (pairs['total_count_5p'] + pairs['total_count_3p'])/2
        pairs['spanning_prop'] = pairs['spanning_count'] / pairs['total']

        prop_cutoff = self.max_junction_spanning_prop
        pairs = pairs.query('spanning_prop < @prop_cutoff')

        pairs = pairs[self.get_header_list()]

        return pairs


    def contains_large_insertion(self, read, pos_3p, pos_5p):


        large_insert_3p = self.identify_large_insertion_at_site(pos_3p, read)
        large_insert_5p = self.identify_large_insertion_at_site(pos_5p, read)

        if large_insert_3p or large_insert_5p:
            return True
        else:
            return False


    def identify_large_insertion_at_site(self, position, read):

        blocks = read.get_blocks()

        if len(blocks) > 1:
            for i in range(len(blocks)-1):
                block1, block2 = blocks[i], blocks[i+1]
                if block1[1] == block2[0]:

                    if position == block1[1]-1:
                        insertion_length = pysamtools.get_insertion_length(position, read)
                        if insertion_length >= self.large_insertion_cutoff:
                            return True
                        else:
                            return False
                    elif position == block2[0]:
                        insertion_length = pysamtools.get_insertion_length(position, read, reverse=True)
                        if insertion_length >= self.large_insertion_cutoff:
                            return True
                        else:
                            return False
                else:
                    return False

        return False


    def block_overlaps_site(self, block, position):

        if block[0] <= position and block[1] > position:
            return True
        else:
            return False