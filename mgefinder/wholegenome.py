from mgefinder import bwatools
from mgefinder import bowtie2tools
from mgefinder import samtools
from mgefinder.find import _find
from mgefinder.pair import _pair
from mgefinder.inferseqassembly import _inferseq_assembly
from Bio import SeqIO
from snakemake import shell
import click


def _wholegenome(reference, query, read_length, read_depth, min_alignment_quality, max_direct_repeat_length,
                 large_insertion_cutoff, query_id, output_prefix):

    min_alignment_inner_length = max_direct_repeat_length + 1

    if not bwatools.genome_is_indexed(reference):
        click.echo("Indexing reference...")
        bwatools.index_genome(reference)
    else:
        click.echo("Reference already indexed...")

    if not bowtie2tools.genome_is_indexed(query):
        click.echo("Indexing query...")
        bowtie2tools.index_genome(query)
    else:
        click.echo("Query already indexed...")

    click.echo("Making query reads...")
    query_reads_path = output_prefix + '.query.tmp.fq'
    make_reads(query, query_reads_path, read_length, read_depth)

    click.echo("Aligning query reads to reference...")
    out_sam = output_prefix + '.query.reference.tmp.sam'
    bwatools.align_to_genome_se(query_reads_path, reference, out_sam, threads=1, verbose=False)

    click.echo("Sorting and indexing alignment file...")
    out_bam = output_prefix + '.query.reference.tmp.bam'
    samtools.sort_coordinate(out_sam, out_bam, delete_in_bam=True)
    samtools.index(out_bam)

    find_file = output_prefix + '.find.tsv'
    _find(out_bam, min_softclip_length=8, min_softclip_count=1, min_alignment_quality=min_alignment_quality,
          min_alignment_inner_length=min_alignment_inner_length, min_distance_to_mate=max_direct_repeat_length + 2,
          min_softclip_ratio=0.01, max_indel_ratio=0.0, large_insertion_cutoff=large_insertion_cutoff,
          min_count_consensus=1, sample_id=query_id, output_file=find_file)

    pair_file = output_prefix + '.pair.tsv'
    _pair(find_file, out_bam, reference, max_direct_repeat_length=max_direct_repeat_length,
          min_alignment_quality=min_alignment_quality, min_alignment_inner_length=min_alignment_inner_length,
          max_junction_spanning_prop=0.01, large_insertion_cutoff=large_insertion_cutoff, output_file=pair_file)

    inferseq_file = output_prefix + '.inferseq.tsv'
    _inferseq_assembly(pair_file, out_bam, query, reference, min_perc_identity=0.95,
                       max_internal_softclip_prop=0.01, max_inferseq_size=500000,
                       min_inferseq_size=30, keep_intermediate=False, output_file=inferseq_file)

    shell('rm %s %s %s' % (query_reads_path, out_bam, out_bam+'.bai'))


def make_reads(genome, out_fastq, readlength, depth):

    with open(out_fastq, 'w') as out:
        step = int(readlength / depth)
        recnum = 0
        for rec in SeqIO.parse(genome, 'fasta'):
            for i in range(0, len(rec.seq)-readlength, step):
                recnum += 1
                outread = rec.seq[i:i + readlength]
                out.write('@SEQ' + str(recnum) + '\n')
                out.write(str(outread) + '\n+\n')
                out.write('F'*len(outread) + '\n')
    return out_fastq




