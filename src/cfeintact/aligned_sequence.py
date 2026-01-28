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
        seq_reference: Sequence[str] = str(reference.seq or "")
        seq_query: Sequence[str] = str(query.seq or "")
        return Cigar.from_msa(reference=seq_reference, query=seq_query)

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


def create_aligned_sequence(this: SeqRecord, reference: SeqRecord, use_mappy: bool = False, preset: str = "asm20"):
    """
    Factory function to create an aligned sequence using either MAFFT or mappy.
    
    Args:
        this: Query sequence to align
        reference: Reference sequence
        use_mappy: If True, use mappy/minimap2; if False, use MAFFT (default: False)
        preset: Minimap2 preset to use when use_mappy=True (default: "asm20")
    
    Returns:
        AlignedSequence or MappyAlignedSequence instance
    
    Example:
        >>> # Use MAFFT (default)
        >>> aligned = create_aligned_sequence(query, reference)
        >>> 
        >>> # Use mappy
        >>> aligned = create_aligned_sequence(query, reference, use_mappy=True)
    """
    if use_mappy:
        from cfeintact.mappy_aligned_sequence import MappyAlignedSequence
        return MappyAlignedSequence(this=this, reference=reference, preset=preset)
    else:
        return AlignedSequence(this=this, reference=reference)
