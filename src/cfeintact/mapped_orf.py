from dataclasses import dataclass
from fractions import Fraction

from cfeintact.original_orf import OriginalORF


@dataclass
class MappedORF:
    reference: OriginalORF
    query: OriginalORF
    orientation: str
    distance: Fraction
