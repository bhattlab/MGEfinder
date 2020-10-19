import warnings
warnings.filterwarnings("ignore")
import sys
from mgefinder import __version__
from mgefinder.workflow import _workflow
from mgefinder.find import _find
from mgefinder.pair import _pair
from mgefinder.inferseqassembly import _inferseq_assembly
from mgefinder.inferseqoverlap import _inferseq_overlap
from mgefinder.inferseqreference import _inferseq_reference
from mgefinder.inferseqdatabase import _inferseq_database
from mgefinder.makedatabase import _makedatabase
from mgefinder.formatbam import _formatbam
from mgefinder.recall import _recall
from mgefinder.help import CustomHelp
from mgefinder.clusterseq import _clusterseq
from mgefinder.genotype import _genotype
from mgefinder.summarize import _summarize
from mgefinder.makefasta import _makefasta
from mgefinder.wholegenome import _wholegenome
from mgefinder.dependencies import check_dependencies
from mgefinder.misc import aligned_bwa, BWACheckError

import click
from os.path import join, dirname, abspath, basename
from shutil import copyfile

WORKFLOW_DENOVO_ORIGINAL_SNAKEFILE = join(dirname(__file__), 'workflow/denovo.original.Snakefile')
WORKFLOW_DENOVO_ORIGINAL_CONFIG = join(dirname(__file__), 'workflow/denovo.original.config.yml')
WORKFLOW_DENOVO_SENSITIVE_SNAKEFILE = join(dirname(__file__), 'workflow/denovo.sensitive.Snakefile')
WORKFLOW_DENOVO_SENSITIVE_CONFIG = join(dirname(__file__), 'workflow/denovo.sensitive.config.yml')

WORKFLOW_DATABASE_ORIGINAL_SNAKEFILE = join(dirname(__file__), 'workflow/database.original.Snakefile')
WORKFLOW_DATABASE_ORIGINAL_CONFIG = join(dirname(__file__), 'workflow/database.original.config.yml')
WORKFLOW_DATABASE_SENSITIVE_SNAKEFILE = join(dirname(__file__), 'workflow/database.sensitive.Snakefile')
WORKFLOW_DATABASE_SENSITIVE_CONFIG = join(dirname(__file__), 'workflow/database.sensitive.config.yml')

check_dependencies()

def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('MGEfinder version: ' + __version__)
    ctx.exit()

@click.group(cls=CustomHelp, invoke_without_command=True)
@click.option('--version', is_flag=True, callback=print_version, expose_value=False, is_eager=True, help="Get MGEfinder version.")
def cli():
    """Command-line tools to identify mobile genetic element insertions from short-read sequencing data."""
    click.echo("Get help documentation with --help.")
    click.echo("Get version with --version.")
    pass

@cli.group(cls=CustomHelp)
def workflow():
    """Choose a workflow for running mgefinder"""
    pass

@workflow.command(short_help='Run mgefinder on a working directory using the de novo workflow.', help_priority=1)
@click.argument('workdir', type=click.Path(exists=True))
@click.option('--cores', '-t', default=1, help="The number of processors to run while finding flank extensions. default=1")
@click.option('--memory', '-m', default=16000, help="Memory limit in megabytes. default=16000; 0 for unlimited")
@click.option('--unlock/--no-unlock',  default=False, help="Unlock working directory if necessary.")
@click.option('--rerun-incomplete/--no-rerun-incomplete',  default=False, help="Rerun incomplete files in the workflow.")
@click.option('--keep-going/--no-keep-going',  default=False, help="Keep going with independent jobs if one fails.")
@click.option('--sensitive/--original',  default=False, help="Run in original mode (default) or sensitive mode.")
def denovo(workdir, cores, memory, unlock, rerun_incomplete, keep_going, sensitive):
    """A click access point for the workflow module. This is used for creating the command line interface."""

    snakefile = WORKFLOW_DENOVO_ORIGINAL_SNAKEFILE
    configfile = WORKFLOW_DENOVO_ORIGINAL_CONFIG
    if sensitive:
        snakefile = WORKFLOW_DENOVO_SENSITIVE_SNAKEFILE
        configfile = WORKFLOW_DENOVO_SENSITIVE_CONFIG

    log_params(command='workflow', workdir=workdir, cores=cores, memory=memory,
                 unlock=unlock, rerun_incomplete=rerun_incomplete, keep_going=keep_going, sensitive=sensitive)
    _workflow(workdir, snakefile, configfile, cores, memory, unlock, rerun_incomplete, keep_going)

@workflow.command(short_help='Run mgefinder on a working directory using the database workflow.', help_priority=2)
@click.argument('workdir', type=click.Path(exists=True))
@click.option('--cores', '-t', default=1, help="The number of processors to run while finding flank extensions. default=1")
@click.option('--memory', '-m', default=16000, help="Memory limit in megabytes. default=16000; 0 for unlimited")
@click.option('--unlock/--no-unlock',  default=False, help="Unlock working directory if necessary.")
@click.option('--rerun-incomplete/--no-rerun-incomplete',  default=False, help="Rerun incomplete files in the workflow.")
@click.option('--keep-going/--no-keep-going',  default=False, help="Keep going with independent jobs if one fails.")
@click.option('--sensitive/--original',  default=False, help="Run in original mode (default) or sensitive mode.")
def database(workdir, cores, memory, unlock, rerun_incomplete, keep_going, sensitive):
    """A click access point for the workflow module. This is used for creating the command line interface."""

    snakefile = WORKFLOW_DATABASE_ORIGINAL_SNAKEFILE
    configfile = WORKFLOW_DATABASE_ORIGINAL_CONFIG
    if sensitive:
        snakefile = WORKFLOW_DATABASE_SENSITIVE_SNAKEFILE
        configfile = WORKFLOW_DATABASE_SENSITIVE_CONFIG

    log_params(command='workflow', workdir=workdir, cores=cores, memory=memory,
                 unlock=unlock, rerun_incomplete=rerun_incomplete, keep_going=keep_going, sensitive=sensitive)
    _workflow(workdir, snakefile, configfile, cores, memory, unlock, rerun_incomplete, keep_going)

@workflow.command(short_help='Run mgefinder on a working directory using a custom workflow.', help_priority=3)
@click.argument('workdir', type=click.Path(exists=True))
@click.option('--snakefile', '-s', type=click.Path(exists=True), required=True)
@click.option('--configfile', '-s', type=click.Path(exists=True), required=True)
@click.option('--cores', '-t', default=1, help="The number of processors to run while finding flank extensions. default=1")
@click.option('--memory', '-m', default=16000, help="Memory limit in megabytes. default=16000; 0 for unlimited")
@click.option('--unlock/--no-unlock',  default=False, help="Unlock working directory if necessary.")
@click.option('--rerun-incomplete/--no-rerun-incomplete',  default=False, help="Rerun incomplete files in the workflow.")
@click.option('--keep-going/--no-keep-going',  default=False, help="Keep going with independent jobs if one fails.")
def custom(workdir, snakefile, configfile, cores, memory, unlock, rerun_incomplete, keep_going):
    """A click access point for the workflow module. This is used for creating the command line interface."""

    log_params(command='workflow', snakefile=snakefile, configfile=configfile, workdir=workdir, cores=cores, memory=memory,
                 unlock=unlock, rerun_incomplete=rerun_incomplete, keep_going=keep_going)
    _workflow(workdir, snakefile, configfile, cores, memory, unlock, rerun_incomplete, keep_going)


@cli.command(short_help='Get copies of all of the default Snakefile and config files.', help_priority=2)
def getworkflow():
    """A click access point for the getworkflow module. This is used for creating the command line interface."""

    click.echo("Copying Snakefile and config file to current working directory...")
    copyfile(WORKFLOW_DENOVO_ORIGINAL_SNAKEFILE, basename(WORKFLOW_DENOVO_ORIGINAL_SNAKEFILE))
    copyfile(WORKFLOW_DENOVO_ORIGINAL_CONFIG, basename(WORKFLOW_DENOVO_ORIGINAL_CONFIG))
    copyfile(WORKFLOW_DENOVO_SENSITIVE_SNAKEFILE, basename(WORKFLOW_DENOVO_SENSITIVE_SNAKEFILE))
    copyfile(WORKFLOW_DENOVO_SENSITIVE_CONFIG, basename(WORKFLOW_DENOVO_SENSITIVE_CONFIG))

    copyfile(WORKFLOW_DATABASE_ORIGINAL_SNAKEFILE, basename(WORKFLOW_DATABASE_ORIGINAL_SNAKEFILE))
    copyfile(WORKFLOW_DATABASE_ORIGINAL_CONFIG, basename(WORKFLOW_DATABASE_ORIGINAL_CONFIG))
    copyfile(WORKFLOW_DATABASE_SENSITIVE_SNAKEFILE, basename(WORKFLOW_DATABASE_SENSITIVE_SNAKEFILE))
    copyfile(WORKFLOW_DATABASE_SENSITIVE_CONFIG, basename(WORKFLOW_DATABASE_SENSITIVE_CONFIG))
    click.echo("Done.")



@cli.command(short_help='Find insertion sites and reconstruct termini of inserted sequence', help_priority=3)
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--min_softclip_length', '-minlen', default=8, help="For a softclipped site to be considered, there must be at least one softclipped read of this length. default=8")
@click.option('--min_softclip_count', '-mincount', default=2, help="For a softclipped site to be considered, there must be at least this many softclipped reads at the site. default=2")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff. default=20")
@click.option('--min_alignment_inner_length', '-minial', default=21, help="If a read is softclipped on both ends, the aligned portion must be at least this long. Ideally, set this equal to 1 + max_direct_repeat_length. default=21")
@click.option('--min_distance_to_mate', '-mindist', default=22, help="A minimum distance to a potential nearby mate, filters out sites that have no pairs. default=22")
@click.option('--min_softclip_ratio', '-minratio', default=0.10, help="For a softclipped site to be considered, the proportion of softclipped sites + large insertions must not fall below this value. default=0.10")
@click.option('--max_indel_ratio', '-maxir', default=0.03, help="For a softclipped site to be considered, the proportion of small insertions/deletions at this site must not be above this value. default=0.03")
@click.option('--large_insertion_cutoff', '-lins', default=30, help="Keep large insertions if they meet this length. default=30")
@click.option('--min_count_consensus', '-mcc', default=2, help="When building the consensus sequence, stop building consensus if read count drops below this cutoff. default=2")
@click.option('--sample_id', '-id', default=None, help="Specify an ID to use for this sample, default is the absolute path to the bam file")
@click.option('--check-bwa/--no-check-bwa', default=True, help="Check to ensure that BWA MEM was used to generate BAM file. default=True")
@click.option('--output_file', '-o', default='mgefinder.find.tsv', help="The output file to save the results. default=mgefinder.find.tsv")
def find(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, min_alignment_inner_length,
               min_distance_to_mate, min_softclip_ratio, max_indel_ratio, large_insertion_cutoff,
               min_count_consensus, sample_id, check_bwa, output_file):
    """A click access point for the find module. This is used for creating the command line interface."""

    if sample_id is None:
        sample_id = abspath(bamfile)

    if check_bwa and not aligned_bwa(bamfile):
        raise BWACheckError("ERROR: Please use BWA MEM to align your reads to the reference genome. Other aligners may not be compatible. To disable this check, use the --no-check-bwa flag.")
        sys.exit()
    

    log_params(command='find',
                 bamfile=bamfile, min_softclip_length=min_softclip_length, min_softclip_count=min_softclip_count,
                 min_alignment_quality=min_alignment_quality, min_alignment_inner_length=min_alignment_inner_length,
                 min_distance_to_mate=min_distance_to_mate, min_softclip_ratio=min_softclip_ratio,
                 max_indel_ratio=max_indel_ratio, large_insertion_cutoff=large_insertion_cutoff,
                 min_count_consensus=min_count_consensus, sample_id=sample_id, check_bwa=check_bwa, 
                 output_file=output_file)
    _find(bamfile, min_softclip_length, min_softclip_count, min_alignment_quality, min_alignment_inner_length,
                min_distance_to_mate, min_softclip_ratio, max_indel_ratio, large_insertion_cutoff,
                min_count_consensus, sample_id, output_file)


@cli.command(short_help="Pair identified termini with each other to represent 5' and 3' ends of inserted sequence.", help_priority=4)
@click.argument('findfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.argument('genome', type=click.Path(exists=True))
@click.option('--max_direct_repeat_length', '-maxdr', default=20, help="The maximum length of a direct repeat to consider a pair. default=20")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff. default=20")
@click.option('--min_alignment_inner_length', '-minial', default=21, help="If a read is softclipped on both ends, the aligned portion must be at least this long. Ideally, set this equal to 1 + maximum direct repeat length. default=21")
@click.option('--max_junction_spanning_prop', '-maxjsp', default=0.15, help="Removes pairs where this proportion of readsextend across both insertion junctions without softclipping, an indication that the site is a duplicated region. default=0.15")
@click.option('--large_insertion_cutoff', '-lins', default=30, help="Keep large insertions if they meet this length. default=30")
@click.option('--output_file', '-o', default='mgefinder.pairs.tsv', help="The output file to save the results. default=mgefinder.pairs.tsv")
def pair(findfile, bamfile, genome, max_direct_repeat_length, min_alignment_quality,
               min_alignment_inner_length, max_junction_spanning_prop, large_insertion_cutoff, output_file=None):

    log_params(command='pair', findfile=findfile, bamfile=bamfile, genome=genome,
                 max_direct_repeat_length=max_direct_repeat_length, min_alignment_quality=min_alignment_quality,
                 min_alignment_inner_length=min_alignment_inner_length,
                 max_junction_spanning_prop=max_junction_spanning_prop, large_insertion_cutoff=large_insertion_cutoff,
                 output_file=output_file)
    _pair(findfile, bamfile, genome, max_direct_repeat_length, min_alignment_quality,
          min_alignment_inner_length, max_junction_spanning_prop, large_insertion_cutoff, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by aligning terminus pairs to an assembled genome.', help_priority=5)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.argument('inferseq_assembly', type=click.Path(exists=True))
@click.argument('inferseq_reference', type=click.Path(exists=True))
@click.option('--min_perc_identity', '-minident', default=0.95, help="Only consider matches with a percentage identity above this threshold. default=0.95")
@click.option('--max_internal_softclip_prop', '-maxclip', default=0.05, help="Do not consider matches with internal softclipped ends exceeding this proportion of the total read. default=0.05")
@click.option('--max_inferseq_size', '-maxsize', default=200000, help="Do not consider inferred sequences over this size. default=200000")
@click.option('--min_inferseq_size', '-minsize', default=30, help="Do not consider inferred sequences below this size. default=30")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files. default=False")
@click.option('--output_file', '-o', default='mgefinder.inferseq_assembly.tsv', help="The output file to save the results. default=mgefinder.inferseq_assembly.tsv")
def inferseq_assembly(pairsfile, bamfile, inferseq_assembly, inferseq_reference, min_perc_identity,
                      max_internal_softclip_prop, max_inferseq_size, min_inferseq_size, keep_intermediate, output_file=None):
    """Infers the identity of an inserted sequence by aligning terminus pairs to an assembled genome."""

    log_params(command='inferseq_assembly', pairsfile=pairsfile, bamfile=bamfile, inferseq_assembly=inferseq_assembly,
                 inferseq_reference=inferseq_reference, min_perc_identity=min_perc_identity,
                      max_internal_softclip_prop=max_internal_softclip_prop, max_inferseq_size=max_inferseq_size,
                 min_inferseq_size=min_inferseq_size, keep_intermediate=keep_intermediate, output_file=output_file)
    _inferseq_assembly(pairsfile, bamfile, inferseq_assembly, inferseq_reference, min_perc_identity,
                       max_internal_softclip_prop, max_inferseq_size, min_inferseq_size, keep_intermediate, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by aligning terminus pairs to a reference genome. Ideal for re-sequencing experiments where evolved strains are closely related to the reference genome used.',
             help_priority =6)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('inferseq_reference', type=click.Path(exists=True))
@click.option('--min_perc_identity', '-minident', default=0.95, help="Only consider matches with a percentage identity above this threshold. default=0.95")
@click.option('--max_internal_softclip_prop', '-maxclip', default=0.05, help="Do not consider matches with internal softclipped ends exceeding this proportion of the total read. default=0.05")
@click.option('--max_inferseq_size', '-maxsize', default=200000, help="Do not consider inferred sequences over this size. default=200000")
@click.option('--min_inferseq_size', '-minsize', default=30, help="Do not consider inferred sequences below this size. default=30")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files. default=False")
@click.option('--output_file', '-o', default='mgefinder.inferseq_reference.tsv', help="The output file to save the results. default=mgefinder.inferseq_reference.tsv")
def inferseq_reference(pairsfile, inferseq_reference, min_perc_identity, max_internal_softclip_prop,
                       max_inferseq_size, min_inferseq_size, keep_intermediate, output_file=None):
    """
    Infers the identity of an inserted sequence by aligning terminus pairs to a reference genome.
    Ideal for re-sequencing experiments where evolved strains are closely related to the reference genome used.
    """
    log_params(command='inferseq_reference', pairsfile=pairsfile, inferseq_reference=inferseq_reference,
                 min_perc_identity=min_perc_identity, max_internal_softclip_prop=max_internal_softclip_prop,
                 max_inferseq_size=max_inferseq_size, min_inferseq_size=min_inferseq_size,
                 keep_intermediate=keep_intermediate, output_file=output_file)
    _inferseq_reference(pairsfile, inferseq_reference, min_perc_identity, max_internal_softclip_prop,
                        max_inferseq_size, min_inferseq_size, keep_intermediate, output_file)


@cli.command(short_help='Infers the identity of an inserted sequence by checking if they overlap with one another. Only identifies an inserted sequence if the consensus termini are long enough to span the entire insertion.',
             help_priority=7)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.option('--min_overlap_score', '-minscore', default=10, help="The minimum overlap score to keep inferred sequence. default=10")
@click.option('--min_overlap_perc_identity', '-minopi', default=0.9, help="The minimum overlap percent identity to keep inferred sequence. default=0.9")
@click.option('--min_inferseq_size', '-minsize', default=30, help="Do not consider inferred sequences below this size. default=30")
@click.option('--output_file', '-o', default='mgefinder.inferseq_overlap.tsv', help="The output file to save the results. default=mgefinder.inferseq_overlap.tsv")
def inferseq_overlap(pairsfile, min_overlap_score, min_overlap_perc_identity, min_inferseq_size, output_file=None):
    """
    Infers the identity of an inserted sequence by checking if they overlap with one another.
    Only identifies an inserted sequence if the consensus termini are long enough to span the entire insertion.
    """
    log_params(command='inferseq_overlap', pairsfile=pairsfile, min_overlap_score=min_overlap_score,
                 min_overlap_perc_identity=min_overlap_perc_identity, min_inferseq_size=min_inferseq_size,
                 output_file=output_file)
    _inferseq_overlap(pairsfile, min_overlap_score, min_overlap_perc_identity, min_inferseq_size, output_file)


@cli.command(short_help="Make database from inferred sequences.", help_priority=8)
@click.argument('inferseqfiles', type=click.Path(exists=True), nargs=-1)
@click.option('--minimum_size', '-minsize', default=30, help="The minimum size of inferred sequence to use in the database. default=30")
@click.option('--maximum_size', '-maxsize', default=200000, help="The maximum size of inferred sequence to use in the database. default=200000")
@click.option('--threads', '-t', default=1, help="The number of processors to use clustering sequences. default=1")
@click.option('--memory', '-m', default=16000, help="Memory limit in megabytes. default=16000; 0 for unlimited")
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--output_dir', '-o', default='mgefinder.database', help="The output directory to save the results. default=mgefinder.database")
@click.option('--prefix', '-p', default='mgefinder.database', help="The prefix used to name the database files. default=mgefinder.database")
def makedatabase(inferseqfiles, minimum_size, maximum_size, threads, memory, force, output_dir=None, prefix=None):
    """
    Makes a FASTA database of all elements identified in provided inferseq files.
    """
    log_params(command='makedatabase', inferseqfiles=inferseqfiles, minimum_size=minimum_size,
                 maximum_size=maximum_size, threads=threads, memory=memory, force=force, output_dir=output_dir,
                 prefix=prefix)
    _makedatabase(inferseqfiles, minimum_size, maximum_size, threads, memory, force, output_dir, prefix)


@cli.command(short_help='Infers the identity of an inserted sequence by aligning flank pairs to an database of known inserted elements.', help_priority=9)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('inferseq_database', type=click.Path(exists=True))
@click.option('--min_perc_identity', '-minident', default=0.90, help="Only consider matches with a percentage identity above this threshold. default=0.90")
@click.option('--max_internal_softclip_prop', '-maxclip', default=0.05, help="Do not consider matches with internal softclipped ends exceeding this proportion of the total read. default=0.05")
@click.option('--max_edge_distance', '-maxedgedist', default=10, help="Reads must align within this number of bases from the edge of an element to be considered. default=10")
@click.option('--output_file', '-o', default='mgefinder.inferseq_database.tsv', help="The output file to save the results. default=mgefinder.inferseq_database.tsv")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files. default=False")
def inferseq_database(pairsfile, inferseq_database, min_perc_identity,  max_internal_softclip_prop, max_edge_distance, output_file=None, keep_intermediate=False):
    """Infers the identity of an inserted sequence by aligning flank pairs to a database of known inserted elements."""
    log_params(command='inferseq_database', pairsfile=pairsfile, inferseq_database=inferseq_database, 
                 min_perc_identity=min_perc_identity, max_internal_softclip_prop=max_internal_softclip_prop, 
                 max_edge_distance=max_edge_distance, output_file=output_file, keep_intermediate=keep_intermediate)
    _inferseq_database(pairsfile, inferseq_database, min_perc_identity, max_internal_softclip_prop, max_edge_distance, output_file, keep_intermediate)


@cli.command(short_help="Cluster inferred sequences.", help_priority=10)
@click.argument('inferseqfiles', type=click.Path(exists=True), nargs=-1)
@click.option('--minimum_size', '-minsize', default=70, help="The minimum size of inferred sequence to use in the database. default=70")
@click.option('--maximum_size', '-maxsize', default=200000, help="The maximum size of inferred sequence to use in the database. default=200000")
@click.option('--threads', '-t', default=1, help="The number of processors to use clustering sequences. default=1")
@click.option('--memory', '-m', default=16000, help="Memory limit in megabytes. default=16000; 0 for unlimited")
@click.option('--output_file', '-o', default='metamgefinder.clusterseq.tsv', help="The output file to save the results. default=metamgefinder.clusterseq.tsv")
def clusterseq(inferseqfiles, minimum_size, maximum_size, threads, memory, output_file=None):
    """Clusters elements identified in inferseq files using CD-HIT-EST"""
    log_params(command='clusterseq', inferseqfiles=inferseqfiles, minimum_size=minimum_size, maximum_size=maximum_size,
               threads=threads, memory=memory, output_file=output_file)
    _clusterseq(inferseqfiles, minimum_size, maximum_size, threads, memory, output_file)



@cli.command(short_help="Assign final genotypes to samples.", help_priority=11)
@click.argument('clusterseq', type=click.Path(exists=True), nargs=1)
@click.argument('pairfiles', type=click.Path(exists=True), nargs=-1)
@click.option('--filter-clusters-inferred-assembly/--no-filter-clusters-inferred-assembly', default=True, help="Removes sequence clusters that were never inferred from an assembly or by overlap if True. default=True")
@click.option('--output_file', '-o', default='mgefinder.genotype.tsv', help="The output file to save the results. default=mgefinder.genotype.tsv")
def genotype(clusterseq, pairfiles, filter_clusters_inferred_assembly, output_file):
    """Genotypes MGE insertions of input samples with respect to the reference genome."""
    log_params(command='genotype', clusterseq=clusterseq, pairfiles=pairfiles,
               filter_clusters_inferred_assembly=filter_clusters_inferred_assembly,
               output_file=output_file)
    _genotype(clusterseq, pairfiles, filter_clusters_inferred_assembly, output_file)


@cli.command(short_help="Summarize the clusters identified by mgefinder.", help_priority=12)
@click.argument('clusterseq', type=click.Path(exists=True), nargs=1)
@click.argument('genotypes', type=click.Path(exists=True), nargs=1)
@click.option('--output_prefix', '-o', default='mgefinder.summarize', help="The output file to save the results. default=mgefinder.summarize")
def summarize(clusterseq, genotypes, output_prefix):
    """Summarizes the findings of MGEfinder workflow"""
    log_params(command='summarize', clusterseq=clusterseq, genotypes=genotypes, output_prefix=output_prefix)
    _summarize(clusterseq, genotypes,  output_prefix)


@cli.command(short_help="Make FASTA files of the identified sequence clusters.", help_priority=13)
@click.argument('clusterseq', type=click.Path(exists=True), nargs=1)
@click.argument('summarize_clusters', type=click.Path(exists=True), nargs=1)
@click.option('--output_prefix', '-o', default='mgefinder', help="The output file to save the results. default=mgefinder")
def makefasta(clusterseq, summarize_clusters, output_prefix):
    """Makes FASTA files of the putative integrative MGEs identified in the pipeline"""
    log_params(command='makefasta', clusterseq=clusterseq, summarize_clusters=summarize_clusters,
               output_prefix=output_prefix)
    _makefasta(clusterseq, summarize_clusters,  output_prefix)


@cli.command(short_help="Formats a BAM file for use with mgefinder. Usually not necessary, unless using the experiment extendpairs command.", help_priority=14)
@click.argument('in_sam', type=click.Path(exists=True))
@click.argument('out_bam')
@click.option('--single-end', is_flag=True, default=False, help="Add this flag for single-end files. default=False")
@click.option('--keep-tmp-files', is_flag=True, default=False, help="Add this flag if you want to keep intermediate temporary files. default=False")
def formatbam(in_sam, out_bam, single_end, keep_tmp_files):
    """Formats a BAM file for use with mgefinder. Usually not necessary, unless using the experiment extendpairs command."""
    log_params(command='formatbam', in_sam=in_sam, out_bam=out_bam, single_end=single_end, keep_tmp_files=keep_tmp_files)
    _formatbam(in_sam, out_bam, single_end, keep_tmp_files)



@cli.command(short_help='Recall softclip counts and runthrough counts from BAM file at specified insertion sites', help_priority=15)
@click.argument('pairsfile', type=click.Path(exists=True))
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff. default=20")
@click.option('--min_alignment_inner_length', '-minial', default=21, help="If a read is softclipped on both ends, the aligned portion must be at least this long. Ideally, set this equal to 1 + maximum direct repeat length. default=21")
@click.option('--large_insertion_cutoff', '-lins', default=30, help="Same parameter as that used for find command. default=30")
@click.option('--output_file', '-o', default='mgefinder.recall.tsv', help="The output file to save results to. default=mgefinder.recall.tsv")
def recall(pairsfile, bamfile, min_alignment_quality, min_alignment_inner_length, large_insertion_cutoff, output_file):
    """Recall softclip counts and runthrough counts from BAM file at specified insertions sites."""
    log_params(command='recall', pairsfile=pairsfile, bamfile=bamfile, min_alignment_quality=min_alignment_quality,
               min_alignment_inner_length=min_alignment_inner_length, large_insertion_cutoff=large_insertion_cutoff,
               output_file=output_file)
    _recall(pairsfile, bamfile, min_alignment_quality, min_alignment_inner_length, large_insertion_cutoff, output_file)

@cli.command(short_help='Identify insertions using whole-genome alignments of complete genomes.', help_priority=16)
@click.argument('reference', type=click.Path(exists=True))
@click.argument('query', type=click.Path(exists=True))
@click.option('--read_length', '-rl', default=500, help="The length of the reads to simulate. Sets an effective minimum size for identified insertions. default=100")
@click.option('--read_depth', '-rd', default=5, help="The depth of coverage produced to compare genomes. default=5")
@click.option('--min_alignment_quality', '-minq', default=20, help="For a read to be considered, it must meet this alignment quality cutoff. default=20")
@click.option('--max_direct_repeat_length', '-maxdr', default=20, help="Only call insertions where the predicted direct repeat length is at most {max_direct_repeat_length}. default=20.")
@click.option('--large_insertion_cutoff', '-lins', default=30, help="Same parameter as that used for find command. default=30") 
@click.option('--query_id', '-id', default=None, help="Specify an ID to use for this query, default is the absolute path to the query file.")
@click.option('--output_prefix', '-o', default='mgefinder.wholegenome', help="The prefix of saved output files. default=mgefinder.wholegenomev")
def wholegenome(reference, query, read_length, read_depth, min_alignment_quality, max_direct_repeat_length, large_insertion_cutoff, query_id, output_prefix):
    """Identify insertions using whole-genome alignments of complete genomes. Inputs can be draft assemblies or completed genomes"""
    if query_id is None:
        query_id = abspath(query)

    log_params(command='wholegenome', reference=reference, query=query,
               read_length=read_length, read_depth=read_depth,
               min_alignment_quality=min_alignment_quality,
               max_direct_repeat_length=max_direct_repeat_length,
               large_insertion_cutoff=large_insertion_cutoff,
               query_id=query_id, output_prefix=output_prefix)
    _wholegenome(reference, query, read_length, read_depth, min_alignment_quality,
                 max_direct_repeat_length, large_insertion_cutoff, query_id, output_prefix)


def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':

    cli()
