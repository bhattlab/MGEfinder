import sys
import pysam
from collections import defaultdict, OrderedDict
from os.path import join
from random import randint
from mgefinder import fastatools, bowtie2tools, pysamtools, sctools, misc
import pygogo as gogo
from Bio import SeqIO
import pandas as pd
from snakemake import shell


verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger

class InferSequence:

    pairs = None
    genome_fasta = None
    method_name = None
    tmp_dir = None
    bam = None

    min_perc_identity = None
    max_internal_softclip_prop = None
    max_inferseq_size = None
    min_inferseq_size = None
    keep_intermediate = None

    all_aligned_pairs = None

    def __init__(self, pairs, genome_fasta, min_perc_identity, max_internal_softclip_prop,
                 max_inferseq_size, min_inferseq_size, keep_intermediate, method_name='inferred_sequence',
                 tmp_dir='/tmp'):

        self.pairs  = pairs
        self.genome_fasta = genome_fasta
        self.genome_dict = {rec.id: rec.seq for rec in SeqIO.parse(genome_fasta, 'fasta')}
        self.method_name = method_name
        self.tmp_dir = tmp_dir

        self.min_perc_identity = min_perc_identity
        self.max_internal_softclip_prop = max_internal_softclip_prop
        self.max_inferseq_size = max_inferseq_size
        self.min_inferseq_size = min_inferseq_size
        self.keep_intermediate = keep_intermediate

        self.all_aligned_pairs = defaultdict(AlignedPairs)

    def infer_sequences(self):

        self.__align_pairs_to_fasta()

        if not self.keep_intermediate:
            shell('rm {bam} {bam}.bai'.format(bam=self.bam.filename.decode('utf-8')))

        self.__prefilter_reads()

        for read in self.bam:
            name = read.query_name.split('_')[0]
            self.all_aligned_pairs[name].add_read(read)

        self.__match_forward_reverse_reads()

        self.__filter_smallest_overlapping()
        self.__filter_pairs_max_internal_softclip_prop()
        self.__filter_best_alignments()
        self.__filter_pairs_size()

        return self.make_dataframe()

    def __align_pairs_to_fasta(self):
        fasta_prefix = join(self.tmp_dir, 'mgefinder.inferseq.' + str(randint(0, 1e20)))

        writetermini = self.get_termini()

        fastatools.write_termini_to_unpaired_fasta(writetermini, fasta_prefix)

        logger.info("Aligning pairs to genome...")
        pair_bam_path = join(self.tmp_dir, 'mgefinder.inferseq.' + str(randint(0, 1e20)) + '.bam')
        bowtie2tools.align_fasta_to_genome(
            fasta_prefix + '.fasta',
            self.genome_fasta, pair_bam_path, silence=True,
            additional_flags='--all --score-min G,1,9'
        )
        self.bam = pysam.AlignmentFile(pair_bam_path, 'rb')


        if not self.keep_intermediate:
            shell('rm {fasta}'.format(fasta=fasta_prefix+'.fasta'))


    def get_termini(self):
        termini = []
        for index, p in self.pairs.iterrows():
            pair_id, seq_5p, seq_3p = p['pair_id'], p['seq_5p'], p['seq_3p']
            seq_5p, seq_3p = seq_5p.upper().rstrip('N'), seq_3p.upper().lstrip('N')
            termini.append({'pair_id': str(pair_id),
                           'seq_5p': seq_5p,
                           'seq_3p': seq_3p})
        return termini

    def __prefilter_reads(self):

        filtered_bam = []

        for read in self.bam:

            if pysamtools.get_perc_identity(read) < self.min_perc_identity:
                continue

            if not read.is_reverse:
                if sctools.is_left_softclipped_strict(read) and \
                    sctools.get_left_softclip_length(read) > 1 and \
                    sctools.left_softclipped_position(read) >= 0:
                    continue

            if read.is_reverse:
                if sctools.is_right_softclipped_strict(read) and \
                                sctools.get_right_softclip_length(read) > 1 and \
                                sctools.right_softclipped_position(read) < len(self.genome_dict[read.reference_name]):
                    continue

            filtered_bam.append(read)

        self.bam = filtered_bam

    def make_dataframe(self):

        outdict = OrderedDict([("pair_id", []), ("method", []), ("loc", []),
                               ("inferred_seq_length", []), ("inferred_seq", [])])

        for aligned_pairs in self.all_aligned_pairs:
            for pair in self.all_aligned_pairs[aligned_pairs].pairs:

                inferred_seq = self.get_inferred_sequence(
                    pair.forward_read, pair.reverse_read, pair.is_reverse()
                )
                outdict['pair_id'].append(int(pair.get_pair_id()))
                outdict['method'].append(self.method_name)
                outdict['loc'].append(pair.get_location())
                outdict['inferred_seq_length'].append(len(inferred_seq))
                outdict['inferred_seq'].append(inferred_seq)

        outdf = pd.DataFrame.from_dict(outdict).sort_values(['pair_id', 'method', 'loc']).reset_index(drop=True)
        return outdf

    def get_inferred_sequence(self, forward_read, reverse_read, is_reverse):
        contig, start, end = forward_read.reference_name, forward_read.reference_start, reverse_read.reference_end
        inferred_sequence = ''.join(self.genome_dict[contig][start:end])

        inferred_sequence = sctools.left_softclipped_sequence_strict(forward_read) + \
                            inferred_sequence + \
                            sctools.right_softclipped_sequence_strict(reverse_read)
        if is_reverse:
            inferred_sequence = misc.revcomp(inferred_sequence)

        return inferred_sequence

    def __match_forward_reverse_reads(self):
        for aligned_pairs in self.all_aligned_pairs:
            self.all_aligned_pairs[aligned_pairs].match_forward_reverse()

    def print_all_pairs(self):
        for aligned_pairs in self.all_aligned_pairs:
            self.all_aligned_pairs[aligned_pairs].print_pairs()

    def __filter_best_alignments(self):
        for aligned_pairs in self.all_aligned_pairs:
            self.all_aligned_pairs[aligned_pairs].filter_best_alignments()

    def __filter_smallest_overlapping(self):
        for aligned_pairs in self.all_aligned_pairs:
            self.all_aligned_pairs[aligned_pairs].filter_smallest_overlapping()

    def __filter_pairs_max_internal_softclip_prop(self):
        for aligned_pairs in self.all_aligned_pairs:
            self.all_aligned_pairs[aligned_pairs].filter_pairs_max_internal_softclip_prop(self.max_internal_softclip_prop)

    def __filter_pairs_size(self):
        for aligned_pairs in self.all_aligned_pairs:
            self.all_aligned_pairs[aligned_pairs].filter_pairs_size(self.max_inferseq_size, self.min_inferseq_size)


class AlignedPairs:


    forward_reads_mate1 = None
    forward_reads_mate2 = None
    forward_reads_mate1_positions = None
    forward_reads_mate2_positions = None

    reverse_reads_mate1 = None
    reverse_reads_mate2 = None
    reverse_reads_mate1_positions = None
    reverse_reads_mate2_positions = None

    pairs = None

    def __init__(self):

        self.forward_reads_mate1 = []
        self.forward_reads_mate2 = []
        self.forward_reads_mate1_positions = defaultdict(dict)
        self.forward_reads_mate2_positions = defaultdict(dict)

        self.reverse_reads_mate1 = []
        self.reverse_reads_mate2 = []
        self.reverse_reads_mate1_positions = defaultdict(dict)
        self.reverse_reads_mate2_positions = defaultdict(dict)

        self.pairs = set()

    def add_read(self, read):

        if read.is_reverse:
            if read.query_name.split('_')[-1] == '1':
                self.reverse_reads_mate1.append(read)
                self.reverse_reads_mate1_positions[read.reference_name][read.reference_end] = read
            else:
                self.reverse_reads_mate2.append(read)
                self.reverse_reads_mate2_positions[read.reference_name][read.reference_end] = read
        else:
            if read.query_name.split('_')[-1] == '1':
                self.forward_reads_mate1.append(read)
                self.forward_reads_mate1_positions[read.reference_name][read.reference_start] = read
            else:
                self.forward_reads_mate2.append(read)
                self.forward_reads_mate2_positions[read.reference_name][read.reference_start] = read


    def match_forward_reverse(self):

        added = set()
        for read in self.forward_reads_mate1 + self.forward_reads_mate2:
            closest_reverse_read = self.get_closest_reverse_read(read)
            if closest_reverse_read is not None:

                name = (read.query_name, read.reference_name, read.reference_start, read.reference_end,
                        closest_reverse_read.query_name, closest_reverse_read.reference_name,
                        closest_reverse_read.reference_start, closest_reverse_read.reference_end)

                if name not in added:
                    self.pairs.add(AlignedPair(read, closest_reverse_read))
                    added.add(name)

        for read in self.reverse_reads_mate1 + self.reverse_reads_mate2:
            closest_forward_read = self.get_closest_forward_read(read)
            if closest_forward_read is not None:
                name = (closest_forward_read.query_name, closest_forward_read.reference_name,
                        closest_forward_read.reference_start, closest_forward_read.reference_end,
                        read.query_name, read.reference_name, read.reference_start, read.reference_end,)

                if name not in added:
                    self.pairs.add(AlignedPair(closest_forward_read, read))
                    added.add(name)

    def get_closest_reverse_read(self, read):

        closest_reverse_read = None
        if read.query_name.split('_')[-1] == '1':
            positions = sorted(list(self.reverse_reads_mate2_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestLarger(positions, read.reference_start)
            if closest_pos is not None:
                closest_reverse_read = self.reverse_reads_mate2_positions[read.reference_name][closest_pos]
        else:
            positions = sorted(list(self.reverse_reads_mate1_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestLarger(positions, read.reference_start)
            if closest_pos is not None:
                closest_reverse_read = self.reverse_reads_mate1_positions[read.reference_name][closest_pos]

        return closest_reverse_read


    def get_closest_forward_read(self, read):

        closest_forward_read = None
        if read.query_name.split('_')[-1] == '1':
            positions = sorted(list(self.forward_reads_mate2_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestSmaller(positions, read.reference_end)
            if closest_pos is not None:
                closest_forward_read = self.forward_reads_mate2_positions[read.reference_name][closest_pos]
        else:
            positions = sorted(list(self.forward_reads_mate1_positions[read.reference_name].keys()))
            closest_pos = misc.takeClosestSmaller(positions, read.reference_end)
            if closest_pos is not None:
                closest_forward_read = self.forward_reads_mate1_positions[read.reference_name][closest_pos]

        return closest_forward_read


    def has_both_mates(self, read1, read2):

        has_one = False
        has_two = False

        if read1.query_name.split('_')[-1] == '1' or read2.query_name.split('_')[-1] == '1':
            has_one = True
        if read1.query_name.split('_')[-1] == '2' or read2.query_name.split('_')[-1] == '2':
            has_two = True

        return has_one and has_two


    def filter_best_alignments(self):
        best_alignment_score = -sys.maxsize
        for p in self.pairs:
            align_score = p.forward_read.get_tag('AS') + p.reverse_read.get_tag('AS')
            if align_score > best_alignment_score:
                best_alignment_score = align_score

        keep_pairs = list()
        for p in self.pairs:
            align_score = p.forward_read.get_tag('AS') + p.reverse_read.get_tag('AS')
            if align_score == best_alignment_score:
                keep_pairs.append(p)

        self.pairs = keep_pairs


    def filter_smallest_overlapping(self):
        forward_sites = defaultdict(lambda: defaultdict(set))
        reverse_sites = defaultdict(lambda: defaultdict(set))
        for p in self.pairs:
            forward_sites[p.forward_read.query_name][p.forward_read.reference_start].add(p)
            reverse_sites[p.reverse_read.query_name][p.reverse_read.reference_end].add(p)

        keep_pairs = set()
        for p in self.pairs:

            forward_sites_pairs = forward_sites[p.forward_read.query_name][p.forward_read.reference_start]
            reverse_sites_pairs = reverse_sites[p.reverse_read.query_name][p.reverse_read.reference_end]

            if len(forward_sites_pairs) == 1 and len(reverse_sites_pairs) == 1:
                keep_pairs.add(p)

            elif len(forward_sites_pairs) > 1 and len(reverse_sites_pairs) == 1:
                smallest_window_p = self.get_smallest_window_pair(forward_sites_pairs)
                keep_pairs.add(smallest_window_p)

            elif len(forward_sites_pairs) == 1 and len(reverse_sites_pairs) > 1:
                smallest_window_p = self.get_smallest_window_pair(reverse_sites_pairs)
                keep_pairs.add(smallest_window_p)

            elif len(forward_sites_pairs) > 1 and len(reverse_sites_pairs) > 1:
                smallest_window_p = self.get_smallest_window_pair(forward_sites_pairs.union(reverse_sites_pairs))
                keep_pairs.add(smallest_window_p)
            else:
                print()
                print("MISSED SOMETHING - DEBUG")
                sys.exit()

        self.pairs = list(keep_pairs)

    def get_smallest_window_pair(self, pairs):
        smallest_pair = None
        minlength = sys.maxsize
        for p in pairs:
            if p.get_aligned_pair_length() < minlength:
                minlength = p.get_aligned_pair_length()
                smallest_pair = p
        return smallest_pair

    def filter_pairs_max_internal_softclip_prop(self, max_internal_softclip_prop):
        keep_pairs = list()
        for p in self.pairs:

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

    def filter_pairs_size(self, max_inferseq_size, min_inferseq_size):
        keep_pairs = list()
        for p in self.pairs:

            if p.get_aligned_pair_length() < min_inferseq_size or p.get_aligned_pair_length() > max_inferseq_size:
                continue
            keep_pairs.append(p)

        self.pairs = keep_pairs

    def print_pairs(self):

        for p in self.pairs:
            print(p)

class AlignedPair:
    forward_read = None
    reverse_read = None


    def __init__(self, forward_read, reverse_read):
        self.forward_read = forward_read
        self.reverse_read = reverse_read

    def get_alignment_score(self):
        return self.forward_read.get_tag('AS') + self.reverse_read.get_tag('AS')

    def get_aligned_pair_length(self):
        return self.reverse_read.reference_end - self.forward_read.reference_start

    def get_pair_id(self):
        return self.forward_read.query_name.split('_')[0]

    def get_location(self):
        return '{contig}:{start}-{end}'.format(
            contig=self.forward_read.reference_name,
            start=self.forward_read.reference_start,
            end=self.reverse_read.reference_end
        )

    def get_inferred_sequence_contig(self):
        return self.forward_read.reference_name

    def get_inferred_sequence_start(self):
        return self.forward_read.reference_start

    def get_inferred_sequence_end(self):
        return self.reverse_read.reference_end

    def is_reverse(self):
        return self.forward_read.query_name.split('_')[-1] == '2'

    def __str__(self):
        out = '{fname} - {fcontig} {fstart} {fend}, {rname} - {rcontig} {rstart} {rend}, LENGTH={length}'.format(
            fname = self.forward_read.query_name, fcontig = self.forward_read.reference_name,
            fstart = self.forward_read.reference_start, fend = self.forward_read.reference_end,
            rname = self.reverse_read.query_name, rcontig = self.reverse_read.reference_name,
            rstart = self.reverse_read.reference_start, rend = self.reverse_read.reference_end,
            length = self.get_aligned_pair_length()
        )

        return out


