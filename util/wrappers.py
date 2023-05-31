import os
import subprocess
import sys
import uuid

from Bio import SeqIO, AlignIO

def mafft(working_dir, sequences):
    '''
    Call mafft on a set of sequences and return the resulting alignment.

    Args:
        sequences: Sequences to be aligned

    Returns:
        An AlignIO object
    '''

    alignment_input = os.path.join(working_dir, str(uuid.uuid4()) + ".fasta")
    alignment_output = os.path.join(working_dir, str(uuid.uuid4()) + ".fasta")

    SeqIO.write(sequences, alignment_input, "fasta")

    with open(alignment_output, 'w') as output:
        subprocess.call(["mafft", "--quiet", alignment_input], shell=False, stdout=output)

    alignment = AlignIO.read(alignment_output, "fasta")

    os.unlink(alignment_input)
    os.unlink(alignment_output)

    return alignment
