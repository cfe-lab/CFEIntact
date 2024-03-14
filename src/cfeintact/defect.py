
from dataclasses import dataclass
from typing import Union

from cfeintact.original_orf import OriginalORF


@dataclass(frozen=True)
class LongDeletion:
    def __str__(self):
        return "Query sequence contains a long deletion."


@dataclass(frozen=True)
class DeletionInOrf:
    e: OriginalORF
    q: OriginalORF
    is_small: bool
    deletions: int

    def __str__(self):
        return (f"{'Smaller ' if self.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start}-{self.q.end}"
                f" can have maximum deletions "
                f"{self.e.deletion_tolerence}, got {self.deletions}")


@dataclass(frozen=True)
class InsertionInOrf:
    q: OriginalORF
    e: OriginalORF
    is_small: bool
    insertions: int

    def __str__(self):
        return (f"{'Smaller ' if self.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start}-{self.q.end} can have maximum insertions "
                f"{3 * self.e.deletion_tolerence}, got {self.insertions}")


@dataclass(frozen=True)
class InternalStopInOrf:
    e: OriginalORF
    q: OriginalORF
    is_small: bool
    position: int

    def __str__(self):
        return (f"{'Smaller ' if self.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start}-{self.q.end}"
                f" contains an internal stop codon at {self.position}")


@dataclass(frozen=True)
class FrameshiftInOrf:
    e: OriginalORF
    q: OriginalORF
    is_small: bool
    impacted_positions: int

    def __str__(self):
        return (f"{'Smaller ' if self.is_small else ''}"
                f"ORF {self.e.name} at {self.e.start}-{self.e.end}"  # FIXME: change this to q.{start,end}
                f" contains out of frame indels that impact {self.impacted_positions} positions.")


@dataclass(frozen=True)
class MajorSpliceDonorSiteMutated:
    splice_site: str

    def __str__(self):
        adj = 'missing' if all(x == '-' for x in self.splice_site) else 'mutated'
        return f"Query sequence has a {adj} splice donor site, {self.splice_site}."


@dataclass(frozen=True)
class PackagingSignalDeletion:
    deletions: int
    tolerance: int

    def __str__(self):
        return (f"Query Sequence exceeds maximum deletion tolerance in PSI. "
                f"Contains {self.deletions} deletions with max tolerance of {self.tolerance} deletions.")


@dataclass(frozen=True)
class PackagingSignalNotComplete:
    position: int

    def __str__(self):
        return ("Query does not encompass the complete PSI region. "
                f"PSI starts at reference position {self.position}.")


@dataclass(frozen=True)
class RevResponseElementDeletion:
    deletions: int
    tolerance: int

    def __str__(self):
        return (f"Query Sequence exceeds maximum deletion tolerance in RRE. "
                f"Contains {self.deletions} deletions with max tolerance of {self.tolerance} deletions.")


@dataclass(frozen=True)
class APOBECHypermutationDetected:
    p_value: float

    def __str__(self):
        return (f"Query sequence shows evidence of APOBEC3F/G-mediated hypermutation "
                f"(p = {self.p_value}).")


@dataclass(frozen=True)
class NonHIV:
    def __str__(self):
        return "Sequence contains unrecognized parts. It is probably a Human/HIV Chimera sequence."


@dataclass(frozen=True)
class Scramble:
    direction: str

    def __str__(self):
        return f"Sequence is {self.direction}-scrambled."


@dataclass(frozen=True)
class InternalInversion:
    def __str__(self):
        return "Sequence contains an internal inversion."


@dataclass(frozen=True)
class UnknownNucleotide:
    details: str

    def __str__(self):
        return f"Sequence contains invalid nucleotides: {self.details}"


# The final exported Defect type will encompass all defined classes
Defect = Union[
    LongDeletion, DeletionInOrf, InsertionInOrf, InternalStopInOrf,
    FrameshiftInOrf, MajorSpliceDonorSiteMutated, PackagingSignalDeletion,
    PackagingSignalNotComplete, RevResponseElementDeletion, APOBECHypermutationDetected,
    NonHIV, Scramble, InternalInversion, UnknownNucleotide
]
