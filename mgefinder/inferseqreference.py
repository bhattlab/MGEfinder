import sys
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from mgefinder import bowtie2tools
from mgefinder.inferseq import InferSequence
import pygogo as gogo
from os.path import dirname

verbose=True
logger = gogo.Gogo(__name__, verbose=verbose).logger


import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def _inferseq_reference(pairsfile, inferseq_reference, min_perc_identity, max_internal_softclip_prop,
                        max_inferseq_size, min_inferseq_size, keep_intermediate, output_file):

    pairs = pd.read_csv(pairsfile, sep='\t', keep_default_na=False, na_values=[
        '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A', 'N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan'])


    handle_empty_pairsfile(pairs, output_file)

    index_genome(inferseq_reference)
    tmp_dir = dirname(output_file)

    inferer = InferSequence(
        pairs, inferseq_reference, min_perc_identity, max_internal_softclip_prop, max_inferseq_size,
        min_inferseq_size, keep_intermediate, 'inferred_reference', tmp_dir
    )

    inferred_sequences = inferer.infer_sequences()

    logger.info("Inferred sequences for %d pairs..." % len(set(list(inferred_sequences['pair_id']))))
    logger.info("Writing results to file %s..." % output_file)
    if pairs.shape[0] > 0:
        sample_id = list(pairs['sample'])[0]
        inferred_sequences.insert(0, 'sample', sample_id)
    else:
        inferred_sequences.insert(0, 'sample', None)

    inferred_sequences.to_csv(output_file, sep='\t', index=False)


def index_genome(inferseq_reference):
    if not bowtie2tools.genome_is_indexed(inferseq_reference):
        logger.info("Indexing inferseq reference genome...")
        bowtie2tools.index_genome(inferseq_reference)
    logger.info("Genome has been indexed...")


def handle_empty_pairsfile(pairs, output_file):
    if pairs.shape[0] == 0:
        outfile = pd.DataFrame(columns=['pair_id', 'method', 'loc', 'inferred_seq_length', 'inferred_seq'])

        if not output_file:
            output_file = 'mgefinder.inferseq_reference.tsv'

        outfile.to_csv(output_file, sep='\t', index=False)
        logger.info("Empty pairs file, exiting...")
        sys.exit()


if __name__ == '__main__':
    _inferseq_reference()