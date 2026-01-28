import dataclasses
import subprocess
import tempfile
from typing import Iterable

from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cfeintact.blastrow import BlastRow
from cfeintact.user_error import UserError


def global_align(sequences: Iterable[SeqRecord]) -> MultipleSeqAlignment:
    '''
    Perform global pairwise alignment using BioPython's PairwiseAligner.

    Args:
        sequences: Two sequences to be aligned

    Returns:
        A MultipleSeqAlignment object
    '''
    seq_list = list(sequences)
    if len(seq_list) != 2:
        raise UserError("Global alignment requires exactly 2 sequences, got %d", len(seq_list))
    
    reference, query = seq_list
    
    # Use BioPython's PairwiseAligner for global alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = Align.substitution_matrices.load("NUC.4.4")
    
    # Perform alignment and get the first (best) alignment
    try:
        alignments = aligner.align(str(reference.seq), str(query.seq))
        alignment = next(iter(alignments))
        
        # Convert to MultipleSeqAlignment
        aligned_ref = SeqRecord(Seq(alignment[0]), id=reference.id, description=reference.description)
        aligned_query = SeqRecord(Seq(alignment[1]), id=query.id, description=query.description)
        
        return MultipleSeqAlignment([aligned_ref, aligned_query])
    except StopIteration:
        raise UserError("Global alignment failed - no alignment found")
    except Exception as e:
        raise UserError("Global alignment failed with %s", e) from e


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
            raise UserError("BLAST run failed with %s. "
                            "Please check if the input FASTA file is correctly formatted.",
                            e) from e
        except KeyboardInterrupt as e:
            raise UserError("BLAST run got cancelled. Cannot continue the analysis.") from e

        output.write(','.join(fields) + '\n')
        alignment_output.seek(0)
        for line in alignment_output:
            output.write(line.decode('utf-8'))
