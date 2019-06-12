from Bio import SeqIO
from snakemake import shell
from mgefinder import misc

def cluster_100p(infile, outfile):

    seqdict1 = dict()
    for rec in SeqIO.parse(infile, 'fasta'):
        seq = str(rec.seq)

        if seq not in seqdict1:
            seqdict1[seq] = rec

    seqdict2 = dict()
    for seq in seqdict1:
        fwd = seq
        rev = misc.revcomp(seq)

        seq_key = tuple(sorted([fwd, rev]))

        if seq_key not in seqdict2:
            seqdict2[seq_key] = seqdict1[seq]

    with open(outfile, "w") as handle:
        SeqIO.write(seqdict2.values(), handle, "fasta")


def cluster_reciprocal_identity(infile, outfile, threads=1, memory=800, alignment_coverage=0.99, perc_identity=0.99):


    shell('cd-hit-est -g 1 -aL {cov} -aS {cov} -d 0 -c {perc_identity} -T {threads} '
          '-M {memory} -i {infile} -o {outfile}'.format(
        cov = alignment_coverage, perc_identity = perc_identity, threads=threads, memory=memory,
        infile = infile, outfile = outfile
    ))


def parse_cdhit_output_100p(unclustered_fasta, cdhit_output):


    final_cluster_mappings = dict()

    seqdict_100p = dict()
    for rec in SeqIO.parse(cdhit_output, 'fasta'):
        fwd = str(rec.seq)
        rev = misc.revcomp(fwd)

        seq_key = tuple(sorted([fwd, rev]))
        seqdict_100p[seq_key] = rec


    for rec in SeqIO.parse(unclustered_fasta, 'fasta'):
        fwd = str(rec.seq)
        rev = misc.revcomp(fwd)

        seq_key = tuple(sorted([fwd, rev]))

        cluster = seqdict_100p[seq_key].id

        final_cluster_mappings[rec.id] = cluster

    return final_cluster_mappings


def parse_cdhit_output_perc_identity(cdhit_output):


    final_cluster_mappings = dict()
    with open(cdhit_output+'.clstr') as infile:
        all_cluster_seqs = []
        repr_seq = None
        for l in infile:
            l = l.strip()
            if l.startswith('>'):
                for seq in all_cluster_seqs:
                    final_cluster_mappings[seq] = repr_seq
                all_cluster_seqs = []
                repr_seq = None
            else:
                seq = l.split()[2].strip('>').strip('.')
                all_cluster_seqs.append(seq)
                if l[-1] == '*':
                    repr_seq = seq

        for seq in all_cluster_seqs:
            final_cluster_mappings[seq] = repr_seq

    return final_cluster_mappings