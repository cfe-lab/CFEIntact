
from Bio import Align
from dataclasses import dataclass

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2


@dataclass(frozen=True)
class Alignment:
    reference: str
    query: str
    score: float

    def __getitem__(self, index: int) -> str:
        if index == 0:
            return self.reference
        elif index == 1:
            return self.query
        else:
            raise IndexError("Index out of range")

    def distance(self) -> float:
        denominator = max(1, len(self.reference))
        shift: int = aligner.match_score
        absolute = (-1 * self.score) / denominator + shift
        return absolute


def align(reference: str, query: str) -> Alignment:
    if reference and query:
        x = aligner.align(reference, query)[0]
        return Alignment(str(x[0]), str(x[1]), x.score)
    elif reference or query:
        return align(reference or "-" * len(query), query or "-" * len(reference))
    else:
        return Alignment(reference, query, 0)
