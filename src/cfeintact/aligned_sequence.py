from dataclasses import dataclass
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from aligntools import Cigar
from functools import cached_property
from typing import Sequence

import cfeintact.wrappers as wrappers
from Bio.Align import MultipleSeqAlignment


@dataclass(frozen=True)
class AlignedSequence:
    this: SeqRecord
    reference: SeqRecord

    @cached_property
    def alignment(self) -> MultipleSeqAlignment:
        return wrappers.global_align([self.reference, self.this])

    @property
    def aligned_reference(self) -> SeqRecord:
        ret: SeqRecord = self.alignment[0]
        return ret

    @property
    def aligned_this(self) -> SeqRecord:
        ret: SeqRecord = self.alignment[1]
        return ret

    @cached_property
    def cigar(self) -> Cigar:
        reference = self.aligned_reference
        query = self.aligned_this
        seq_reference: Sequence[str] = reference.seq or ""  # type: ignore
        seq_query: Sequence[str] = query.seq or ""  # type: ignore
        return Cigar.from_msa(reference=seq_reference, query=seq_query)  # type: ignore

    @cached_property
    def coordinate_mapping(self):
        return self.cigar.coordinate_mapping

    def reverse(self):
        seq = self.this.seq
        assert seq is not None
        newthis = SeqRecord(Seq.reverse_complement(seq),
                            id=self.this.id,
                            name=self.this.name
                            )

        return AlignedSequence(this=newthis, reference=self.reference)

    def alignment_score(self) -> int:
        left = self.aligned_reference.seq or ""
        right = self.aligned_this.seq or ""
        return sum([a == b for a, b in zip(left, right)])
