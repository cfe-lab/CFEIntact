from cfeintact.translate_to_aminoacids import translate_to_aminoacids
from typing import Tuple
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord


def get_query_aminoacids_table(sequence: SeqRecord) -> Tuple[str, str, str]:
    if not hasattr(sequence.seq, "query_aminoacids_table"):
        seq = sequence.seq
        assert isinstance(seq, Seq) or isinstance(seq, MutableSeq)
        ret = (translate_to_aminoacids(seq, frame=0),
               translate_to_aminoacids(seq, frame=1),
               translate_to_aminoacids(seq, frame=2))
        sequence.seq.query_aminoacids_table = ret  # type: ignore
    else:
        ret = sequence.seq.query_aminoacids_table  # type: ignore

    return ret
