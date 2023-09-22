import Bio
from util.translate_to_aminoacids import translate_to_aminoacids

def get_query_aminoacids_table(sequence):
    if not hasattr(sequence.seq, "query_aminoacids_table"):
        try:
            sequence.seq.query_aminoacids_table = [translate_to_aminoacids(sequence.seq, i) for i in range(3)]
        except Bio.Data.CodonTable.TranslationError as e:
            sequence.seq.query_aminoacids_table = e

    return sequence.seq.query_aminoacids_table
