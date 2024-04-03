
from dataclasses import dataclass
from typing import Union
from fractions import Fraction

from cfeintact.original_orf import OriginalORF


@dataclass(frozen=True)
class ORFDefect:
    q: OriginalORF


@dataclass(frozen=True)
class LongDeletion:
    def __str__(self) -> str:
        return "Query sequence contains a long deletion."


@dataclass(frozen=True)
class DeletionInOrf(ORFDefect):
    e: OriginalORF
    deletions: int

    def __str__(self) -> str:
        return (f"{'Smaller ' if self.q.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start + 1}-{self.q.end + 1}"
                f" can have maximum deletions "
                f"{self.e.max_deletions}, got {self.deletions}.")


@dataclass(frozen=True)
class InsertionInOrf(ORFDefect):
    e: OriginalORF
    insertions: int

    def __str__(self) -> str:
        return (f"{'Smaller ' if self.q.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start + 1}-{self.q.end + 1} can have maximum insertions "
                f"{self.e.max_insertions}, got {self.insertions}.")


@dataclass(frozen=True)
class InternalStopInOrf(ORFDefect):
    e: OriginalORF
    position: int

    def __str__(self) -> str:
        return (f"{'Smaller ' if self.q.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start + 1}-{self.q.end + 1}"
                f" contains an internal stop codon at {self.position + 1}.")


@dataclass(frozen=True)
class FrameshiftInOrf(ORFDefect):
    e: OriginalORF
    impacted_positions: int

    def __str__(self) -> str:
        return (f"{'Smaller ' if self.q.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start + 1}-{self.q.end + 1}"
                f" contains out of frame indels that impact {self.impacted_positions} positions.")


@dataclass(frozen=True)
class SequenceDivergence(ORFDefect):
    e: OriginalORF
    distance: Fraction

    def __str__(self) -> str:
        ex_dist = float(round(self.e.max_distance, 5))
        it_dist = float(round(self.distance, 5))
        return (f"{'Smaller ' if self.q.is_small else ''}"
                f"ORF {self.e.name} at {self.q.start + 1}-{self.q.end + 1} can have maximum distance of "
                f"{ex_dist} from its subtype ORF's aminoacid sequence, got {it_dist}.")


@dataclass(frozen=True)
class MajorSpliceDonorSiteMutated:
    splice_site: str
    context: str

    def __str__(self) -> str:
        adj = 'missing' if all(x == '-' for x in self.splice_site) else 'mutated'
        return f"Query sequence has a {adj} splice donor site, {self.splice_site}, around {self.context}."


@dataclass(frozen=True)
class PackagingSignalDeletion:
    deletions: int
    tolerance: int

    def __str__(self) -> str:
        return (f"Query Sequence exceeds maximum deletion tolerance in PSI. "
                f"Contains {self.deletions} deletions with max tolerance of {self.tolerance} deletions.")


@dataclass(frozen=True)
class PackagingSignalNotComplete:
    position: int

    def __str__(self) -> str:
        return ("Query does not encompass the complete PSI region. "
                f"PSI starts at reference position {self.position + 1}.")


@dataclass(frozen=True)
class RevResponseElementDeletion:
    deletions: int
    tolerance: int

    def __str__(self) -> str:
        return (f"Query Sequence exceeds maximum deletion tolerance in RRE. "
                f"Contains {self.deletions} deletions with max tolerance of {self.tolerance} deletions.")


@dataclass(frozen=True)
class APOBECHypermutationDetected:
    p_value: float

    def __str__(self) -> str:
        return (f"Query sequence shows evidence of APOBEC3F/G-mediated hypermutation "
                f"(p = {self.p_value}).")


@dataclass(frozen=True)
class NonHIV:
    def __str__(self) -> str:
        return "Sequence contains unrecognized parts. It is probably a Human/HIV Chimera sequence."


@dataclass(frozen=True)
class Scramble:
    direction: str

    def __str__(self) -> str:
        return f"Sequence is {self.direction}-scrambled."


@dataclass(frozen=True)
class InternalInversion:
    def __str__(self) -> str:
        return "Sequence contains an internal inversion."


@dataclass(frozen=True)
class UnknownNucleotide:
    details: str

    def __str__(self) -> str:
        return f"Sequence contains invalid nucleotides: {self.details}."


# The final exported Defect type encompasses all defined classes
DefectType = Union[
    LongDeletion, DeletionInOrf, InsertionInOrf, InternalStopInOrf,
    FrameshiftInOrf, SequenceDivergence, MajorSpliceDonorSiteMutated, PackagingSignalDeletion,
    PackagingSignalNotComplete, RevResponseElementDeletion, APOBECHypermutationDetected,
    NonHIV, Scramble, InternalInversion, UnknownNucleotide
]


@dataclass(frozen=True)
class Defect:
    qseqid: str
    error: DefectType
