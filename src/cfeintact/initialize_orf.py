from cfeintact.translate_to_aminoacids import translate_to_aminoacids
from cfeintact.find_orf import find_orf
from cfeintact.original_orf import OriginalORF
from cfeintact.aligned_sequence import AlignedSequence


def initialize_orf(aligned_sequence: AlignedSequence, name: str,
                   start: int, end: int, deletion_tolerence: int,
                   is_small: bool) -> OriginalORF:
    nucleotides = aligned_sequence.reference.seq[start:(end + 1)]
    aminoacids = translate_to_aminoacids(nucleotides)
    protein = aminoacids.strip("*")
    reference_orf = \
        OriginalORF(
            name=name,
            start=start,
            end=end,
            deletion_tolerence=deletion_tolerence,
            nucleotides=nucleotides,
            aminoacids=aminoacids,
            protein=protein,
            is_small=is_small,
        )

    return find_orf(aligned_sequence, reference_orf).query
