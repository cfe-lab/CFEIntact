from dataclasses import dataclass

from cfeintact.original_orf import OriginalORF

@dataclass
class MappedORF:
    reference: OriginalORF
    query: OriginalORF
    orientation: str
    distance: float
