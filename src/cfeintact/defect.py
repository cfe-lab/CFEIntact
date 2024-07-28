
from dataclasses import dataclass
from typing import Union

from cfeintact.original_orf import OriginalORF


@dataclass(frozen=True)
class ORFDefect:
    q: OriginalORF


@dataclass(frozen=True)
class LongDeletion:
    def __str__(self) -> str:
        return "Query sequence contains a very large deletion."


@dataclass(frozen=True)
class Deletion(ORFDefect):
    e: OriginalORF
    deletions: int

    def __str__(self) -> str:
        return (f"ORF {self.e.name!r} exceeds maximum deletion tolerance."
                f" Contains {self.deletions} deletions with max tolerance of {self.e.max_deletions} deletions.")


@dataclass(frozen=True)
class Insertion(ORFDefect):
    e: OriginalORF
    insertions: int

    def __str__(self) -> str:
        return (f"ORF {self.e.name!r} exceeds maximum insertion tolerance."
                f" Contains {self.insertions} insertions with max tolerance of {self.e.max_insertions} insertions.")


@dataclass(frozen=True)
class MutatedStopCodon(ORFDefect):
    e: OriginalORF
    codon: str

    def __str__(self) -> str:
        codon = self.codon.ljust(3, "-")
        return f"ORF {self.e.name!r} has a mutated stop codon: {codon!r}."


@dataclass(frozen=True)
class MutatedStartCodon(ORFDefect):
    e: OriginalORF
    codon: str

    def __str__(self) -> str:
        codon = self.codon.rjust(3, "-")
        return f"ORF {self.e.name!r} has a mutated start codon: {codon!r}."


@dataclass(frozen=True)
class InternalStop(ORFDefect):
    e: OriginalORF
    position: int

    def __str__(self) -> str:
        return (f"ORF {self.e.name!r} at {self.q.start + 1}-{self.q.end + 1}"
                f" contains an internal stop codon at {self.position + 1}.")


@dataclass(frozen=True)
class Frameshift(ORFDefect):
    e: OriginalORF
    impacted_positions: float

    def __str__(self) -> str:
        return (f"ORF {self.e.name!r} at {self.q.start + 1}-{self.q.end + 1}"
                f" contains out of frame indels that impact {self.impacted_positions} positions.")


@dataclass(frozen=True)
class SequenceDivergence(ORFDefect):
    e: OriginalORF
    distance: float

    def __str__(self) -> str:
        ex_dist = float(round(self.e.max_distance, 5))
        it_dist = float(round(self.distance, 5))
        return (f"ORF {self.e.name!r} exceeds maximum distance tolerance. "
                f"It is {it_dist} units of distance away from its reference ORF's aminoacid sequence"
                f" with max tolerance of {ex_dist}.")


@dataclass(frozen=True)
class MajorSpliceDonorSiteMutated:
    splice_site: str
    context: str

    def __str__(self) -> str:
        adj = 'missing' if all(x == '-' for x in self.splice_site) else 'mutated'
        return f"Query sequence has a {adj} splice donor site: {self.splice_site}. The context is {self.context}."


@dataclass(frozen=True)
class PackagingSignalDeletion:
    deletions: int
    tolerance: int

    def __str__(self) -> str:
        return (f"Query sequence exceeds maximum deletion tolerance in PSI."
                f" Contains {self.deletions} deletions with max tolerance of {self.tolerance} deletions.")


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
        return (f"Query Sequence exceeds maximum deletion tolerance in RRE."
                f" Contains {self.deletions} deletions with max tolerance of {self.tolerance} deletions.")


@dataclass(frozen=True)
class APOBECHypermutation:
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
    LongDeletion, Deletion, Insertion, MutatedStartCodon, MutatedStopCodon, InternalStop,
    Frameshift, SequenceDivergence, MajorSpliceDonorSiteMutated, PackagingSignalDeletion,
    PackagingSignalNotComplete, RevResponseElementDeletion, APOBECHypermutation,
    NonHIV, Scramble, InternalInversion, UnknownNucleotide
]


@dataclass(frozen=True)
class Defect:
    qseqid: str
    error: DefectType
