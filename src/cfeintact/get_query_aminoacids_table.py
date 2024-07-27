from cfeintact.translate_to_aminoacids import translate_to_aminoacids
from typing import Tuple
from Bio.SeqRecord import SeqRecord


def get_query_aminoacids_table(sequence: SeqRecord) -> Tuple[str, str, str]:
    if not hasattr(sequence.seq, "query_aminoacids_table"):
        ret = (translate_to_aminoacids(sequence.seq, frame=0),
               translate_to_aminoacids(sequence.seq, frame=1),
               translate_to_aminoacids(sequence.seq, frame=2))
        sequence.seq.query_aminoacids_table = ret
    else:
        ret = sequence.seq.query_aminoacids_table

    return ret
