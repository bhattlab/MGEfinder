import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from mgefinder.misc import revcomp


def write_sequences_to_fasta(sequences, fasta, names=None):
    outseqs = []
    for i in range(len(sequences)):
        if names:
            name = names[i]
        else:
            name = i
        record = SeqRecord(Seq(sequences[i], IUPAC.IUPACAmbiguousDNA),
                           id=str(name), description=str(name))
        outseqs.append(record)
    with open(fasta, 'w') as output_handle:
        SeqIO.write(outseqs, output_handle, 'fasta')

def read_fasta(fasta):
    for rec in SeqIO.parse(fasta, 'fasta'):
        yield rec

def write_termini_to_fasta(termini, fasta):

    outseqs = []
    for index, row in termini.iterrows():
        name5p = str(index) + '_5p'
        name3p = str(index) + '_3p'
        seq5p = row['seq_5p']
        seq3p = row['seq_3p']

        record5p = SeqRecord(Seq(seq5p, IUPAC.IUPACAmbiguousDNA),
                           id=str(name5p), description=str(name5p))
        record3p = SeqRecord(Seq(seq3p, IUPAC.IUPACAmbiguousDNA),
                             id=str(name3p), description=str(name3p))

        outseqs.append(record5p)
        outseqs.append(record3p)

    with open(fasta, 'w') as output_handle:
        SeqIO.write(outseqs, output_handle, 'fasta')

def write_termini_to_paired_fasta(pairs, fasta_prefix):

    fiveprime_outseqs = []
    threeprime_outseqs = []
    for p in pairs:
        name5p = p['pair_id']
        name3p = p['pair_id']
        seq5p = p['seq_5p']
        seq3p = p['seq_3p']

        record5p = SeqRecord(Seq(seq5p, IUPAC.IUPACAmbiguousDNA),
                           id=str(name5p), description=str(name5p)+'_5p')
        record3p = SeqRecord(Seq(revcomp(seq3p), IUPAC.IUPACAmbiguousDNA),
                             id=str(name3p), description=str(name3p)+'_3p')

        fiveprime_outseqs.append(record5p)
        threeprime_outseqs.append(record3p)

    with open(fasta_prefix+'.R1.fasta', 'w') as output_handle:
        SeqIO.write(fiveprime_outseqs, output_handle, 'fasta')

    with open(fasta_prefix + '.R2.fasta', 'w') as output_handle:
        SeqIO.write(threeprime_outseqs, output_handle, 'fasta')

def write_termini_to_unpaired_fasta(pairs, fasta_prefix):

    fiveprime_outseqs = []
    threeprime_outseqs = []
    for p in pairs:
        name5p = p['pair_id']+'_1'
        name3p = p['pair_id']+'_2'
        seq5p = p['seq_5p']
        seq3p = p['seq_3p']

        record5p = SeqRecord(Seq(seq5p, IUPAC.IUPACAmbiguousDNA),
                           id=str(name5p), description=str(name5p)+'_5p')
        record3p = SeqRecord(Seq(revcomp(seq3p), IUPAC.IUPACAmbiguousDNA),
                             id=str(name3p), description=str(name3p)+'_3p')

        fiveprime_outseqs.append(record5p)
        threeprime_outseqs.append(record3p)

    with open(fasta_prefix+'.fasta', 'w') as output_handle:
        SeqIO.write(fiveprime_outseqs, output_handle, 'fasta')

    with open(fasta_prefix + '.fasta', 'a') as output_handle:
        SeqIO.write(threeprime_outseqs, output_handle, 'fasta')

def write_panisa_termini_to_fasta(termini, fasta):

    outseqs = []
    for index, row in termini.iterrows():
        name5p = str(index) + '_5p'
        name3p = str(index) + '_3p'
        seq5p = row['Left sequence']
        seq3p = row['Right sequence']

        record5p = SeqRecord(Seq(seq5p, IUPAC.IUPACAmbiguousDNA),
                           id=str(name5p), description=str(name5p))
        record3p = SeqRecord(Seq(seq3p, IUPAC.IUPACAmbiguousDNA),
                             id=str(name3p), description=str(name3p))

        outseqs.append(record5p)
        outseqs.append(record3p)

    with open(fasta, 'w') as output_handle:
        SeqIO.write(outseqs, output_handle, 'fasta')