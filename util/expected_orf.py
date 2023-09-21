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
        vpr_defective_insertion_pos = 5772
        start = start if start < vpr_defective_insertion_pos else start - 1
        end = end if end < vpr_defective_insertion_pos else end - 1

        start_s = ReferenceIndex(start - 1).mapto(aligned_sequence) # decrement is needed because original "start" is 1-based.
        end_s = ReferenceIndex(end).mapto(aligned_sequence)

        nucleotides = str(aligned_sequence.this.seq[start_s:end_s])
        aminoacids = translate_to_aminoacids(nucleotides)
        has_start_codon = translate_to_aminoacids(aligned_sequence.this.seq[(start - 1):end]).startswith("M")
        protein = get_biggest_protein(has_start_codon, aminoacids)

        return ExpectedORF(name=name,
                           start=start_s,
                           end=end_s,
                           deletion_tolerence=deletion_tolerence,
                           nucleotides=nucleotides,
                           aminoacids=aminoacids,
                           protein=protein,
                           )

