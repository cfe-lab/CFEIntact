from dataclasses import dataclass
from util.reference_index import ReferenceIndex
from util.translate_to_aminoacids import translate_to_aminoacids
from util.get_biggest_protein import get_biggest_protein


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
        start_s = ReferenceIndex(start).mapto(aligned_sequence)
        end_s = ReferenceIndex(end).mapto(aligned_sequence)

        nucleotides = str(aligned_sequence.this.seq[start_s:(end_s + 1)])
        aminoacids = translate_to_aminoacids(nucleotides)
        has_start_codon = translate_to_aminoacids(aligned_sequence.this.seq[start:(end + 1)]).startswith("M")
        protein = get_biggest_protein(has_start_codon, aminoacids)

        return ExpectedORF(name=name,
                           start=start_s,
                           end=end_s,
                           deletion_tolerence=deletion_tolerence,
                           nucleotides=nucleotides,
                           aminoacids=aminoacids,
                           protein=protein,
                           )

