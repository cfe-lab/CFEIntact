from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class OriginalORF:
    name: str
    start: int
    end: int
    max_deletions: int
    max_insertions: int
    max_distance: Fraction
    nucleotides: str
    aminoacids: str
    protein: str
    is_small: bool
