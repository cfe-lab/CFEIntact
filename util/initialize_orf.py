from util.translate_to_aminoacids import translate_to_aminoacids
from util.find_orf import find_orf
from util.original_orf import OriginalORF


def initialize_orf(aligned_sequence, name, start, end, deletion_tolerence):
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
        )

    return find_orf(aligned_sequence, reference_orf).query
