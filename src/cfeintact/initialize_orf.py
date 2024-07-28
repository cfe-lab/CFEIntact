from cfeintact.translate_to_aminoacids import translate_to_aminoacids
from cfeintact.find_orf import find_orf
from cfeintact.original_orf import OriginalORF
from cfeintact.aligned_sequence import AlignedSequence


def initialize_orf(aligned_sequence: AlignedSequence, name: str,
                   start: int, end: int,
                   max_deletions: int, max_insertions: int,
                   max_distance: float, is_small: bool) -> OriginalORF:
    nucleotides = aligned_sequence.reference.seq[start:(end + 1)]
    aminoacids = translate_to_aminoacids(nucleotides)
    protein = aminoacids.strip("*")
    reference_orf = \
        OriginalORF(
            name=name,
            start=start,
            end=end,
            max_deletions=abs(max_deletions),
            max_insertions=abs(max_insertions),
            max_distance=max_distance,
            nucleotides=nucleotides,
            aminoacids=aminoacids,
            protein=protein,
            is_small=is_small,
        )

    return find_orf(aligned_sequence, reference_orf).query
