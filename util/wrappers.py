import subprocess
import tempfile

from Bio import SeqIO, AlignIO

def mafft(sequences):
    '''
    Call mafft on a set of sequences and return the resulting alignment.

    Args:
        sequences: Sequences to be aligned

    Returns:
        An AlignIO object
    '''

    with tempfile.NamedTemporaryFile() as alignment_input, tempfile.NamedTemporaryFile() as alignment_output:
        SeqIO.write(sequences, alignment_input.name, "fasta")
        subprocess.call(["mafft", "--quiet", alignment_input.name], shell=False, stdout=alignment_output)
        alignment = AlignIO.read(alignment_output.name, "fasta")
        return alignment
