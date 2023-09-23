from dataclasses import dataclass

@dataclass
class OriginalORF:
    name: str
    start: int
    end: int
    deletion_tolerence: int
    nucleotides: str
    aminoacids: str
    protein: str
