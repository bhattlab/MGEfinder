from glob import glob
from snakemake import shell
from os.path import isfile

def index_genome(genome_path, silence=True):

    if silence:
        shell('makeblastdb -dbtype nucl -in {genome_path} &> /dev/null;'.format(genome_path=genome_path))
    else:
        shell('makeblastdb -dbtype nucl -in {genome_path}'.format(genome_path=genome_path))

    return genome_is_indexed(genome_path)


def genome_is_indexed(genome_path):

    nhr = genome_path + '.nhr'
    nin = genome_path + '.nin'
    nsq = genome_path + '.nsq'

    index_files = [nhr, nin, nsq]

    indexed = True
    files = glob(genome_path+'*')
    for f in index_files:
        if f not in files:
            indexed = False

    return indexed


def align_fasta_to_genome(fasta, genome_path, outfile, threads=1, silence=True, additional_flags=''):

    shell('blastn -query {fasta} -db {genome_path} -outfmt 5 -max_target_seqs 100000 -out {outfile} -parse_deflines '
          '{additional_flags}'.format(
        fasta=fasta, genome_path=genome_path, outfile=outfile, additional_flags=additional_flags
    ))

    if isfile(outfile):
        return True
    else:
        return False