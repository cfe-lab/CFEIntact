import dataclasses
import subprocess
import tempfile
from typing import Iterable

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from cfeintact.blastrow import BlastRow
from cfeintact.user_error import UserError


def mafft(sequences: Iterable[SeqRecord]) -> MultipleSeqAlignment:
    '''
    Call mafft on a set of sequences and return the resulting alignment.

    Args:
        sequences: Sequences to be aligned

    Returns:
        An AlignIO object
    '''

    with tempfile.NamedTemporaryFile() as alignment_input, tempfile.NamedTemporaryFile() as alignment_output:
        SeqIO.write(sequences, alignment_input.name, "fasta")

        try:
            subprocess.run(["mafft", "--quiet", alignment_input.name],
                           shell=False, stdout=alignment_output, check=True)
        except Exception as e:
            raise UserError(f"MAFFT run failed with {e}. "
                            f"Please check if the input .FASTA file is correctly formatted.") from e

        alignment: MultipleSeqAlignment = AlignIO.read(alignment_output.name, "fasta")
        return alignment


def blast(alignment_file: str, input_file: str, output_file: str) -> None:
    fileformat = "10"  # .csv format
    fields = [field.name for field in dataclasses.fields(BlastRow)]
    outfmt = fileformat + ' ' + ' '.join(fields)

    with open(output_file, "w") as output, \
            tempfile.NamedTemporaryFile() as alignment_output:

        try:
            subprocess.run(
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
                shell=False,
                check=True,
            )
        except Exception as e:
            raise UserError(f"BLAST run failed with {e}. "
                            f"Please check if the input .FASTA file is correctly formatted.") from e

        output.write(','.join(fields) + '\n')
        alignment_output.seek(0)
        for line in alignment_output:
            output.write(line.decode('utf-8'))
