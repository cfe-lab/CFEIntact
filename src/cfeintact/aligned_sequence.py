from dataclasses import dataclass
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from aligntools import Cigar
from functools import cached_property

import cfeintact.wrappers as wrappers
from Bio.Align import MultipleSeqAlignment


@dataclass(frozen=True)
class AlignedSequence:
    this: SeqRecord
    reference: SeqRecord

    @cached_property
    def alignment(self) -> MultipleSeqAlignment:
        return wrappers.mafft([self.reference, self.this])

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
        return Cigar.from_msa(reference=reference.seq, query=query.seq)

    @cached_property
    def coordinate_mapping(self):
        return self.cigar.coordinate_mapping

    def reverse(self):
        newthis = SeqRecord(Seq.reverse_complement(self.this.seq),
                            id=self.this.id,
                            name=self.this.name
                            )

        return AlignedSequence(this=newthis, reference=self.reference)

    def alignment_score(self) -> int:
        return sum([a == b for a, b in zip(self.aligned_reference.seq, self.aligned_this.seq)])
