from dataclasses import dataclass
from functools import cached_property


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
    region_nucleotides: str
    region_aminoacids: str
    is_small: bool

    @cached_property
    def has_start_codon(self) -> bool:
        return self.region_aminoacids.startswith("M")

    @cached_property
    def has_stop_codon(self) -> bool:
        return self.region_aminoacids.endswith("*")
