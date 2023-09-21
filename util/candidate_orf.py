from dataclasses import dataclass

@dataclass
class CandidateORF:
    name: str
    start: int
    end: int
    subtype_start: int
    subtype_end: int
    orientation: str
    distance: float
    protein: str
    aminoacids: str
