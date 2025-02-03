from dataclasses import dataclass
from Bio.Seq import Seq
from functools import cached_property


@dataclass(frozen=True)
class OriginalORF:
    name: str
    start: int
    end: int
    max_deletions: int
    max_insertions: int
    max_distance: float
    nucleotides: Seq
    aminoacids: str
    protein: str
    is_small: bool

    @cached_property
    def has_start_codon(self) -> bool:
        return self.aminoacids.startswith("M")

    @cached_property
    def has_stop_codon(self) -> bool:
        return self.aminoacids.endswith("*")
