import dataclasses
import subprocess
import tempfile

from Bio import SeqIO, AlignIO

from util.blastrow import BlastRow

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

def blast(alignment_file, input_file, output_file):
    fileformat = "10" # .csv format
    fields = [field.name for field in dataclasses.fields(BlastRow)]
    outfmt = fileformat + ' ' + ' '.join(fields)

    with open(output_file, "w") as output, \
         tempfile.NamedTemporaryFile() as alignment_output:

        subprocess.call(
            ["blastn",
             "-query", input_file,
             "-db", alignment_file,
             "-num_alignments", "1",
             "-reward", "1",
             "-penalty", "-1",
             "-gapopen", "2",
             "-gapextend", "1",
             "-out", alignment_output.name,
             "-outfmt", outfmt],
            shell=False)

        output.write(','.join(fields) + '\n')
        alignment_output.seek(0)
        for line in alignment_output:
            output.write(line.decode('utf-8'))
