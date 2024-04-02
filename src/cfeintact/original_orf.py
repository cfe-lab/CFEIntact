from dataclasses import dataclass


@dataclass(frozen=True)
class OriginalORF:
    name: str
    start: int
    end: int
    max_deletions: int
    max_insertions: int
    max_distance: float
    nucleotides: str
    aminoacids: str
    protein: str
    is_small: bool
