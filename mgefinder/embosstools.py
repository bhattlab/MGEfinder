import warnings
warnings.filterwarnings("ignore")
from snakemake import shell
from mgefinder.fastatools import read_fasta
from pandas import Series

def run_einverted(fasta, gap=12, threshold=15, match=3, mismatch=-4, outfile='einverted.tmp.out', outseq='einverted.tmp.outseq'):
    command = 'einverted -sequence {fasta} -gap {gap} -threshold {threshold} -match {match} ' \
              '-mismatch {mismatch} -outfile {outfile} -outseq {outseq} -auto Y -warning N'.format(
        fasta=fasta, gap=gap, threshold=threshold, match=match, mismatch=mismatch, outfile=outfile, outseq=outseq
    )
    shell(command)


def read_emboss_seq_results(outseq_path):

    pair = []
    for rec in read_fasta(outseq_path):
        d = list(map(int, rec.id.split('_')))
        out = Series({'ir_pos_5p': d[1], 'ir_pos_3p':d[2], 'seq':str(rec.seq)})
        pair.append(out)
        if len(pair) == 2:
            yield pair
            pair = []


