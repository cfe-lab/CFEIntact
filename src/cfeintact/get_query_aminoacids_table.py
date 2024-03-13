import Bio
from cfeintact.translate_to_aminoacids import translate_to_aminoacids

def get_query_aminoacids_table(sequence):
    if not hasattr(sequence.seq, "query_aminoacids_table"):
        sequence.seq.query_aminoacids_table = [translate_to_aminoacids(sequence.seq, i) for i in range(3)]

    return sequence.seq.query_aminoacids_table
