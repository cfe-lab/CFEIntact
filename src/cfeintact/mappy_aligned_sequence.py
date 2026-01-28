"""
Mappy-based implementation of AlignedSequence.

This module provides a drop-in replacement for the MAFFT-based AlignedSequence class,
using mappy (minimap2) for alignment.
"""

from dataclasses import dataclass
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from functools import cached_property
from typing import Optional

from cfeintact.mappy_aligner import MappyAligner, MappyHit
from cfeintact.cigar_parser import CigarCoordinateMapper, AlignToolsCompatibleMapping
from cfeintact.user_error import UserError


@dataclass(frozen=True)
class MappyAlignedSequence:
    """
    Represents a pairwise alignment between a query and reference sequence using mappy.
    
    This class provides a compatible interface with the original AlignedSequence class
    but uses mappy/minimap2 instead of MAFFT for alignment.
    
    Attributes:
        this: Query sequence (the sequence being analyzed)
        reference: Reference sequence (subtype or HXB2)
        preset: Minimap2 preset for alignment (default: "asm20")
    """
    this: SeqRecord
    reference: SeqRecord
    preset: str = "asm20"
    
    @cached_property
    def _mappy_hit(self) -> Optional[MappyHit]:
        """
        Perform alignment and cache the best hit.
        
        Returns:
            Best MappyHit or None if alignment fails
        """
        aligner = MappyAligner(preset=self.preset)
        aligner.index_reference(self.reference)
        return aligner.align_best(self.this)
    
    @cached_property
    def alignment(self) -> MultipleSeqAlignment:
        """
        Generate a MultipleSeqAlignment object from mappy alignment.
        
        This reconstructs a gap-aligned MSA from the CIGAR string for compatibility
        with existing code that expects MAFFT-style alignment output.
        
        Returns:
            BioPython MultipleSeqAlignment object
        """
        hit = self._mappy_hit
        
        if hit is None:
            # No alignment found - return unaligned sequences
            return MultipleSeqAlignment([self.reference, self.this])
        
        # Reconstruct aligned sequences from CIGAR string
        ref_aligned, query_aligned = self._reconstruct_alignment(hit)
        
        # Create SeqRecord objects for the alignment
        ref_record = SeqRecord(
            Seq(ref_aligned),
            id=self.reference.id or "reference",
            name=self.reference.name or "",
            description=""
        )
        
        query_record = SeqRecord(
            Seq(query_aligned),
            id=self.this.id or "query",
            name=self.this.name or "",
            description=""
        )
        
        return MultipleSeqAlignment([ref_record, query_record])
    
    def _reconstruct_alignment(self, hit: MappyHit) -> tuple[str, str]:
        """
        Reconstruct aligned sequences with gaps from CIGAR string.
        
        Args:
            hit: MappyHit containing CIGAR string
        
        Returns:
            Tuple of (reference_aligned, query_aligned) strings with gaps
        """
        from cfeintact.cigar_parser import parse_cigar
        
        ref_seq = str(self.reference.seq)
        query_seq = str(self.this.seq)
        
        ref_aligned = []
        query_aligned = []
        
        ref_pos = hit.ref_start
        query_pos = hit.query_start
        
        cigar_ops = parse_cigar(hit.cigar_str)
        
        for length, op in cigar_ops:
            if op in 'M=X':  # Match or mismatch
                ref_aligned.append(ref_seq[ref_pos:ref_pos + length])
                query_aligned.append(query_seq[query_pos:query_pos + length])
                ref_pos += length
                query_pos += length
            elif op == 'D':  # Deletion from reference (gap in query)
                ref_aligned.append(ref_seq[ref_pos:ref_pos + length])
                query_aligned.append('-' * length)
                ref_pos += length
            elif op == 'I':  # Insertion to reference (gap in reference)
                ref_aligned.append('-' * length)
                query_aligned.append(query_seq[query_pos:query_pos + length])
                query_pos += length
            elif op == 'N':  # Skipped region
                ref_aligned.append(ref_seq[ref_pos:ref_pos + length])
                query_aligned.append('-' * length)
                ref_pos += length
            elif op == 'S':  # Soft clipping (query sequence present but not aligned)
                query_pos += length
            elif op == 'H':  # Hard clipping (query sequence not present)
                pass
        
        return ''.join(ref_aligned), ''.join(query_aligned)
    
    @property
    def aligned_reference(self) -> SeqRecord:
        """
        Get the aligned reference sequence.
        
        Returns:
            Reference sequence from the alignment (with gaps)
        """
        result: SeqRecord = self.alignment[0]
        return result
    
    @property
    def aligned_this(self) -> SeqRecord:
        """
        Get the aligned query sequence.
        
        Returns:
            Query sequence from the alignment (with gaps)
        """
        result: SeqRecord = self.alignment[1]
        return result
    
    @cached_property
    def cigar(self):
        """
        Get CIGAR representation compatible with aligntools.
        
        Returns:
            Cigar-like object (for compatibility, but from mappy)
        """
        # This returns a CIGAR string directly rather than aligntools.Cigar
        # since we're using mappy's CIGAR output
        hit = self._mappy_hit
        if hit is None:
            raise UserError("No alignment found for CIGAR generation")
        
        # Return a simple object that provides the cigar_str attribute
        # for compatibility with code that might access it
        class CigarWrapper:
            def __init__(self, cigar_str, ref_aligned, query_aligned):
                self.cigar_str = cigar_str
                self._ref = ref_aligned
                self._query = query_aligned
            
            @property
            def coordinate_mapping(self):
                # This is accessed by downstream code
                return None  # Will use self.coordinate_mapping directly instead
        
        return CigarWrapper(hit.cigar_str, self.aligned_reference, self.aligned_this)
    
    @cached_property
    def coordinate_mapping(self):
        """
        Get coordinate mapping for position translation.
        
        Returns:
            AlignToolsCompatibleMapping object
        """
        hit = self._mappy_hit
        if hit is None:
            raise UserError("No alignment found for coordinate mapping")
        
        mapper = CigarCoordinateMapper(
            hit.cigar_str,
            ref_start=hit.ref_start,
            query_start=hit.query_start
        )
        
        return AlignToolsCompatibleMapping(mapper)
    
    def reverse(self):
        """
        Create a reversed (reverse complement) version of this alignment.
        
        Returns:
            New MappyAlignedSequence with reverse complemented query
        """
        seq = self.this.seq
        assert seq is not None
        newthis = SeqRecord(
            Seq.reverse_complement(seq),
            id=self.this.id,
            name=self.this.name
        )
        
        return MappyAlignedSequence(
            this=newthis,
            reference=self.reference,
            preset=self.preset
        )
    
    def alignment_score(self) -> int:
        """
        Calculate alignment score based on matching positions.
        
        Returns:
            Number of matching positions in the alignment
        """
        left = self.aligned_reference.seq or ""
        right = self.aligned_this.seq or ""
        return sum([a == b for a, b in zip(left, right)])
