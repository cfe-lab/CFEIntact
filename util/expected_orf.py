from dataclasses import dataclass
from util.translate_to_aminoacids import translate_to_aminoacids
from util.find_orf import find_orf


@dataclass
class ExpectedORF:
    name: str
    start: int
    end: int
    deletion_tolerence: int
    nucleotides: str
    aminoacids: str
    protein: str

    @staticmethod
    def subtyped(aligned_sequence, name, start, end, deletion_tolerence):
        nucleotides = aligned_sequence.reference.seq[start:(end + 1)]
        aminoacids = translate_to_aminoacids(nucleotides)
        protein = aminoacids.strip("*")
        reference_orf = \
            ExpectedORF(
                name=name,
                start=start,
                end=end,
                deletion_tolerence=deletion_tolerence,
                nucleotides=nucleotides,
                aminoacids=aminoacids,
                protein=protein,
            )

        subtype_orf = find_orf(aligned_sequence, reference_orf)
        subtype_orf.deletion_tolerence = deletion_tolerence

        return subtype_orf
