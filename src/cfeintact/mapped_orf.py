from dataclasses import dataclass
from functools import cached_property

import cfeintact.detailed_aligner as detailed_aligner
from cfeintact.original_orf import OriginalORF
from cfeintact.detailed_aligner import Alignment
from cfeintact.get_indel_impact import get_indel_impact


@dataclass(frozen=True)
class MappedORF:
    reference: OriginalORF
    query: OriginalORF
    orientation: str

    @cached_property
    def amino_alignment(self) -> Alignment:
        if self.reference.has_start_codon and self.reference.has_stop_codon:
            return detailed_aligner.align(self.reference.protein, self.query.protein)
        else:
            return detailed_aligner.align(self.reference.aminoacids, self.query.aminoacids)

    @cached_property
    def nuc_alignment(self) -> Alignment:
        got_nucleotides = self.query.nucleotides[:(len(self.query.protein) * 3) + 1].upper()
        exp_nucleotides = self.reference.nucleotides.upper()
        return detailed_aligner.align(exp_nucleotides, got_nucleotides)

    @cached_property
    def distance(self) -> float:
        return self.amino_alignment.distance()

    @cached_property
    def indel_impact(self) -> int:
        return get_indel_impact(self.nuc_alignment[0], self.nuc_alignment[1])

    @cached_property
    def starts_properly(self) -> bool:
        if self.reference.has_start_codon:
            return self.query.has_start_codon
        else:
            return True

    @cached_property
    def ends_properly(self) -> bool:
        if self.reference.has_stop_codon:
            return self.query.has_stop_codon
        else:
            return True
