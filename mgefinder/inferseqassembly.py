import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import sys
import pandas as pd
from mgefinder import bowtie2tools
from mgefinder import sctools
from mgefinder import misc
from mgefinder.inferseq import InferSequence, AlignedPairs
import click
import pysam
from Bio import SeqIO
from collections import OrderedDict
from os.path import dirname
from collections import defaultdict


def _inferseq_assembly(pairsfile, bamfile, inferseq_assembly, inferseq_reference, min_perc_identity,
                       max_internal_softclip_prop, max_inferseq_size, min_inferseq_size, keep_intermediate, output_file):

    pairs = pd.read_csv(pairsfile, sep='\t', keep_default_na=False, na_values=[
        '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A', 'N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan'])

    handle_empty_pairsfile(pairs, output_file)

    index_genome(inferseq_assembly)
    tmp_dir = dirname(output_file)

    context_width = 25
    context_inferer = InferSequenceContext(
        pairs, inferseq_assembly, bamfile, inferseq_reference, min_perc_identity,
        max_internal_softclip_prop, max_inferseq_size,
        min_inferseq_size, keep_intermediate,
        context_width, 'inferred_assembly_with_context', tmp_dir
    )

    inferred_sequences_with_context = context_inferer.infer_sequences()

    nocontext_inferer = InferSequence(
        pairs, inferseq_assembly, min_perc_identity, max_internal_softclip_prop, max_inferseq_size,
        min_inferseq_size, keep_intermediate, 'inferred_assembly_without_context', tmp_dir
    )

    inferred_sequences_without_context = nocontext_inferer.infer_sequences()

    inferred_sequences = pd.concat([inferred_sequences_with_context, inferred_sequences_without_context]).sort_values(
        ['pair_id', 'loc']
    )

    click.echo("Inferred sequences for %d pairs..." % len(set(list(inferred_sequences['pair_id']))))
    click.echo("Writing results to file %s..." % output_file)
    if pairs.shape[0] > 0:
        sample_id = list(pairs['sample'])[0]
        inferred_sequences.insert(0, 'sample', sample_id)
    else:
        inferred_sequences.insert(0, 'sample', None)

    inferred_sequences.to_csv(output_file, sep='\t', index=False)


class InferSequenceContext(InferSequence):

    context_width = None
    ref_bam = None
    ref_genome_dict = None

    def __init__(self, pairs, genome_fasta, ref_bam, ref_genome_fasta, min_perc_identity, max_internal_softclip_prop,
                 max_inferseq_size, min_inferseq_size, keep_intermediate, context_width, method_name='inferred_sequence',
                 tmp_dir='/tmp'):

        InferSequence.__init__(self, pairs, genome_fasta, min_perc_identity, max_internal_softclip_prop,
                               max_inferseq_size, min_inferseq_size, keep_intermediate, method_name, tmp_dir)

        self.ref_bam = pysam.AlignmentFile(ref_bam, 'rb')
        self.context_width = context_width
        self.ref_genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(ref_genome_fasta, 'fasta')}

        self.all_aligned_pairs = defaultdict(lambda: AlignedPairsContext(self.context_width))


    def get_termini(self):

        context_termini = []
        for index, p in self.pairs.iterrows():
            pair_id, contig_name, pos_5p, pos_3p, seq_5p, seq_3p = p['pair_id'], p['contig'], p['pos_5p'], p['pos_3p'], \
                                                                   p['seq_5p'], p['seq_3p']
            seq_5p, seq_3p = seq_5p.rstrip('N'), seq_3p.lstrip('N')
            context_5p = self.get_sequence_context(self.ref_genome_dict[contig_name], contig_name,
                                                   pos_5p - self.context_width, pos_5p)
            context_3p = self.get_sequence_context(self.ref_genome_dict[contig_name], contig_name,
                                                   pos_3p + 1, pos_3p + 1 + self.context_width)

            context_termini.append({'pair_id': str(pair_id),
                                   'seq_5p': context_5p + seq_5p,
                                   'seq_3p': seq_3p + context_3p})

        return context_termini


    def get_sequence_context(self, contig, contig_name, start, end):

        add_start_n = 0
        add_end_n = 0
        if start < 0:
            add_start_n = abs(start)
            start = 0
        if end > len(contig):
            add_end_n = end - len(contig)
            end = len(contig)

        reference_sequence = 'N' * add_start_n + contig[start:end] + 'N' * add_end_n
        sequence_context_dict = initialize_sequence_context(reference_sequence, start, end)

        for read in self.ref_bam.fetch(contig_name, start, end):

            ref_positions = read.get_reference_positions(full_length=True)
            read_query = read.query_sequence
            read_qualities = read.query_qualities

            i = 0
            for pos in ref_positions:
                if pos is not None and pos in sequence_context_dict:
                    sequence_context_dict[pos][read_query[i]] += read_qualities[i]
                i += 1

        consensus = self.get_consensus_context(sequence_context_dict)
        return consensus


    def get_consensus_context(self, sequence_context_dict):
        start, end = min(sequence_context_dict.keys()), max(sequence_context_dict.keys()) + 1
        consensus = ''
        for pos in range(start, end):
            best_qual = 0
            best_nuc = ''
            for nuc in sequence_context_dict[pos]:
                qual = sequence_context_dict[pos][nuc]
                if qual > best_qual:
                    best_qual = qual
                    best_nuc = nuc
            consensus += best_nuc
        return consensus


    def get_inferred_sequence(self, forward_read, reverse_read, is_reverse):
        contig = forward_read.reference_name
        start = forward_read.reference_start
        end = reverse_read.reference_end

        inferred_sequence = ''.join(self.genome_dict[contig][start:end])

        inferred_sequence = sctools.left_softclipped_sequence_strict(forward_read) + \
                            inferred_sequence + \
                            sctools.right_softclipped_sequence_strict(reverse_read)

        inferred_sequence = inferred_sequence[self.context_width:-self.context_width]

        if is_reverse:
            inferred_sequence = misc.revcomp(inferred_sequence)

        contig_edge = False
        if sctools.is_left_softclipped_strict(forward_read) and \
                        sctools.left_softclipped_position(forward_read) < 0:
            contig_edge = True
        elif sctools.is_right_softclipped_strict(reverse_read) and \
                        sctools.right_softclipped_position(reverse_read) >= len(self.genome_dict[contig]):
            contig_edge = True


        return inferred_sequence, contig_edge

    def make_dataframe(self):

        outdict = OrderedDict([("pair_id", []), ("method", []), ("loc", []),
                               ("inferred_seq_length", []), ("inferred_seq", [])])

        for aligned_pairs in self.all_aligned_pairs:
            for pair in self.all_aligned_pairs[aligned_pairs].pairs:

                inferred_seq, contig_edge = self.get_inferred_sequence(
                    pair.forward_read, pair.reverse_read, pair.is_reverse()
                )
                outdict['pair_id'].append(int(pair.get_pair_id()))

                if contig_edge == True:
                    outdict['method'].append('inferred_assembly_with_half_context')
                else:
                    outdict['method'].append('inferred_assembly_with_full_context')

                outdict['loc'].append(pair.get_location())
                outdict['inferred_seq_length'].append(len(inferred_seq))
                outdict['inferred_seq'].append(inferred_seq)

        outdf = pd.DataFrame.from_dict(outdict).sort_values(['pair_id', 'loc']).reset_index(drop=True)
        return outdf


class AlignedPairsContext(AlignedPairs):

    context_width = None


    def __init__(self, context_width):
        AlignedPairs.__init__(self)
        self.context_width = context_width


    def add_read(self, read):

        if read.is_reverse:
            if read.query_name.split('_')[-1] == '1':
                self.reverse_reads_mate1.append(read)
                self.reverse_reads_mate1_positions[read.reference_name][read.reference_end-self.context_width] = read
            else:
                self.reverse_reads_mate2.append(read)
                self.reverse_reads_mate2_positions[read.reference_name][read.reference_end-self.context_width] = read
        else:
            if read.query_name.split('_')[-1] == '1':
                self.forward_reads_mate1.append(read)
                self.forward_reads_mate1_positions[read.reference_name][read.reference_start+self.context_width] = read
            else:
                self.forward_reads_mate2.append(read)
                self.forward_reads_mate2_positions[read.reference_name][read.reference_start+self.context_width] = read


    def filter_pairs_max_internal_softclip_prop(self, max_internal_softclip_prop):
        keep_pairs = list()
        for p in self.pairs:

            if sctools.is_left_softclipped_strict(p.forward_read) and \
                sctools.get_left_softclip_length(p.forward_read) > 1 and \
                sctools.is_right_softclipped_strict(p.reverse_read) and \
                sctools.get_right_softclip_length(p.reverse_read) > 1:
                continue

            if sctools.is_right_softclipped_strict(p.forward_read) and \
                p.forward_read.reference_end < p.reverse_read.reference_end and \
                sctools.right_softclip_proportion(p.forward_read) > max_internal_softclip_prop:
                continue

            if sctools.is_left_softclipped_strict(p.reverse_read) and \
                p.reverse_read.reference_start > p.forward_read.reference_start and \
                sctools.left_softclip_proportion(p.reverse_read) > max_internal_softclip_prop:
                continue

            keep_pairs.append(p)

        self.pairs = keep_pairs


    def get_closest_reverse_read(self, read):

        closest_reverse_read = None
        if read.query_name.split('_')[-1] == '1':
            positions = sorted(list(self.reverse_reads_mate2_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestLarger(positions, read.reference_start+self.context_width)
            if closest_pos is not None:
                closest_reverse_read = self.reverse_reads_mate2_positions[read.reference_name][closest_pos]
        else:
            positions = sorted(list(self.reverse_reads_mate1_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestLarger(positions, read.reference_start+self.context_width)
            if closest_pos is not None:
                closest_reverse_read = self.reverse_reads_mate1_positions[read.reference_name][closest_pos]

        print(read.query_name)
        return closest_reverse_read


    def get_closest_forward_read(self, read):

        closest_forward_read = None
        if read.query_name.split('_')[-1] == '1':
            positions = sorted(list(self.forward_reads_mate2_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestSmaller(positions, read.reference_end-self.context_width)
            if closest_pos is not None:
                closest_forward_read = self.forward_reads_mate2_positions[read.reference_name][closest_pos]
        else:
            positions = sorted(list(self.forward_reads_mate1_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestSmaller(positions, read.reference_end-self.context_width)
            if closest_pos is not None:
                closest_forward_read = self.forward_reads_mate1_positions[read.reference_name][closest_pos]

        return closest_forward_read



def get_inferred_sequences(pairs, genome_dict, add_softclipped_bases=False):

    inferred_sequences = []
    for read1, read2 in pairs:
        if read1.query_name.count('_') == 2:
            context_width = int(read1.query_name.split('_')[-2])
            name = read1.reference_name + ':' + str(read1.reference_start+context_width) + '-' + str(read2.reference_end-context_width)

            inferred_sequence = genome_dict[read1.reference_name][read1.reference_start:read2.reference_end]

            if add_softclipped_bases:
                inferred_sequence = sctools.left_softclipped_sequence_strict(read1) + inferred_sequence + sctools.right_softclipped_sequence_strict(read2)

            inferred_sequence = inferred_sequence[context_width:-context_width]

            if read1.query_name.split('_')[-1] == '2':
                inferred_sequence = misc.revcomp(inferred_sequence)

            contig_edge = False
            if sctools.is_left_softclipped_strict(read1) and \
                sctools.left_softclipped_position(read1) < 0:
                contig_edge = True
            elif sctools.is_right_softclipped_strict(read2) and \
                sctools.right_softclipped_position(read2) >= len(genome_dict[read2.reference_name]):
                contig_edge = True

        else:
            name = read1.reference_name + ':' + str(read1.reference_start) + '-' + str(read2.reference_end)
            inferred_sequence = genome_dict[read1.reference_name][read1.reference_start:read2.reference_end]

            if add_softclipped_bases:
                inferred_sequence = sctools.left_softclipped_sequence_strict(read1) + inferred_sequence + sctools.right_softclipped_sequence_strict(read2)

            if read1.query_name.split('_')[-1] == '2':
                inferred_sequence = misc.revcomp(inferred_sequence)

            contig_edge = False
            if sctools.is_left_softclipped_strict(read1) and \
                            sctools.left_softclipped_position(read1) < 0:
                contig_edge = True
            elif sctools.is_right_softclipped_strict(read2) and \
                            sctools.right_softclipped_position(read2) >= len(genome_dict[read2.reference_name]):
                contig_edge = True

        inferred_sequences.append((name, len(inferred_sequence), contig_edge, inferred_sequence))

    return inferred_sequences


def initialize_sequence_context(target_region, expanded_start, expanded_end):
    target_region_reads = defaultdict(lambda: defaultdict(int))
    i = 0
    for pos in range(expanded_start, expanded_end):
        target_region_reads[pos][target_region[i]] += 1
        i += 1
    return target_region_reads


def index_genome(inferseq_assembly):
    if not bowtie2tools.genome_is_indexed(inferseq_assembly):
        click.echo("Indexing inferseq assembly...")
        bowtie2tools.index_genome(inferseq_assembly)
    click.echo("Genome has been indexed...")


def handle_empty_pairsfile(pairs, output_file):
    if pairs.shape[0] == 0:
        outfile = pd.DataFrame(columns=['pair_id', 'method', 'loc', 'inferred_seq_length', 'inferred_seq'])

        if not output_file:
            output_file = 'mgefinder.inferseq_assembly.tsv'

        outfile.to_csv(output_file, sep='\t', index=False)
        click.echo("Empty pairs file, exiting...")
        sys.exit()
