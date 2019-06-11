from snakemake import shell
from glob import glob
from os.path import isfile
import click

def index_genome(genome_path, silence=True):
    if silence:
        shell('bowtie2-build -o 0 -q {genome_path} {genome_path} 2> /dev/null;'.format(genome_path=genome_path))
    else:
        shell('bowtie2-build -o 0 {genome_path} {genome_path}'.format(genome_path=genome_path))

    return genome_is_indexed(genome_path)


def genome_is_indexed(genome_path):

    one = genome_path + '.1.bt2'
    two = genome_path + '.2.bt2'
    three = genome_path + '.3.bt2'
    four = genome_path + '.4.bt2'
    five = genome_path + '.rev.1.bt2'
    six = genome_path + '.rev.2.bt2'

    index_files = [one, two, three, four, five, six]

    indexed = True
    files = glob(genome_path+'*')
    for f in index_files:
        if f not in files:
            indexed = False

    return indexed


def align_fasta_to_genome(fasta, genome_path, out_bam, threads=1, silence=True, additional_flags=''):

    if silence:
        command = "rm -f {out_bam}.sam {out_bam}.sam.tmp; " \
                  "bowtie2 {additional_flags} --local -x {genome_path} -p {threads} -f -U {fasta} -S {out_bam}.sam 2> /dev/null; " \
                  "samtools view -h -q 1 {out_bam}.sam 2> /dev/null 1> {out_bam}.sam.tmp; " \
                  "samtools sort {out_bam}.sam.tmp 2> /dev/null 1> {out_bam}; " \
                  "samtools index {out_bam} 2> /dev/null;" \
                  "rm -f {out_bam}.sam {out_bam}.sam.tmp".format(
            genome_path=genome_path, fasta=fasta, out_bam=out_bam, threads=threads, additional_flags=additional_flags)

    else:
        command = "rm -f {out_bam}.sam {out_bam}.sam.tmp; " \
                  "bowtie2 {additional_flags} --local -x {genome_path} -p {threads} -f -U {fasta} -S {out_bam}.sam; " \
                  "samtools view -h -q 1 {out_bam}.sam > {out_bam}.sam.tmp; " \
                  "samtools sort {out_bam}.sam.tmp > {out_bam}; " \
                  "samtools index {out_bam};" \
                  "rm -f {out_bam}.sam {out_bam}.sam.tmp".format(
            genome_path=genome_path, fasta=fasta, out_bam=out_bam, threads=threads, additional_flags=additional_flags)



    if not silence:
        click.echo("Executing command: %s" % command)

    shell(command, read=silence)

    if isfile(out_bam):
        return True
    else:
        return False

def align_paired_fasta_to_genome(fasta1, fasta2, genome_path, out_bam, threads=1, silence=False, additional_flags=''):

    if silence:
        command = "rm -f {out_bam}.sam {out_bam}.sam.tmp; " \
                  "bowtie2 {additional_flags} --local -x {genome_path} -p {threads} -f -1 {fasta1} -2 {fasta2} -S {out_bam}.sam &> /dev/null; " \
                  "samtools view -h -q 1 {out_bam}.sam 2> /dev/null 1> {out_bam}.sam.tmp; " \
                  "samtools sort {out_bam}.sam.tmp 2> /dev/null 1> {out_bam}; " \
                  "samtools index {out_bam} 2> /dev/null;" \
                  "rm -f {out_bam}.sam {out_bam}.sam.tmp".format(
            genome_path=genome_path, fasta1=fasta1, fasta2=fasta2, out_bam=out_bam, threads=threads, additional_flags=additional_flags)

    else:
        command = "rm -f {out_bam}.sam {out_bam}.sam.tmp; " \
                  "bowtie2 {additional_flags} --local -x {genome_path} -p {threads} -f -1 {fasta1} -2 {fasta2} -S {out_bam}.sam; " \
                  "samtools view -h -q 1 {out_bam}.sam > {out_bam}.sam.tmp; " \
                  "samtools sort {out_bam}.sam.tmp > {out_bam}; " \
                  "samtools index {out_bam};" \
                  "rm -f {out_bam}.sam {out_bam}.sam.tmp".format(
            genome_path=genome_path, fasta1=fasta1, fasta2=fasta2, out_bam=out_bam, threads=threads, additional_flags=additional_flags)



    if not silence:
        click.echo("Executing command: %s" % command)

    shell(command, read=silence)

    if isfile(out_bam):
        return True
    else:
        return False