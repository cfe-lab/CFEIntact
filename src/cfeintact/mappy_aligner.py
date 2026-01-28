"""
Mappy (minimap2) wrapper for sequence alignment in CFEIntact.

This module provides a high-performance alternative to MAFFT for sequence alignment
using minimap2 through its Python bindings (mappy).
"""

import mappy  # type: ignore
from typing import Optional, Tuple, List, Any
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass

from cfeintact.user_error import UserError


@dataclass
class MappyHit:
    """
    Represents a single alignment hit from mappy.
    
    Attributes:
        query_start: Start position in query (0-based)
        query_end: End position in query (0-based, exclusive)
        ref_start: Start position in reference (0-based)
        ref_end: End position in reference (0-based, exclusive)
        strand: 1 for forward, -1 for reverse
        mapq: Mapping quality
        cigar_str: CIGAR string representation
        match_length: Number of matching bases
        alignment_length: Total alignment length including gaps
    """
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    strand: int
    mapq: int
    cigar_str: str
    match_length: int
    alignment_length: int
    
    @classmethod
    def from_mappy_alignment(cls, hit: mappy.Alignment) -> "MappyHit":
        """
        Create MappyHit from a mappy.Alignment object.
        
        Args:
            hit: mappy.Alignment object
        
        Returns:
            MappyHit instance
        """
        return cls(
            query_start=hit.q_st,
            query_end=hit.q_en,
            ref_start=hit.r_st,
            ref_end=hit.r_en,
            strand=hit.strand,
            mapq=hit.mapq,
            cigar_str=hit.cigar_str,
            match_length=hit.mlen,
            alignment_length=hit.blen,
        )


class MappyAligner:
    """
    Wrapper for mappy alignment operations.
    
    This class provides a convenient interface for aligning sequences using minimap2
    through the mappy Python bindings.
    """
    
    # Preset configurations optimized for different sequence types
    PRESETS = {
        "asm5": "Assembly alignment with ~5% divergence",
        "asm10": "Assembly alignment with ~10% divergence",
        "asm20": "Assembly alignment with ~20% divergence",
        "map-ont": "Oxford Nanopore long reads",
        "map-pb": "PacBio long reads",
        "sr": "Short reads",
        "splice": "Splice-aware alignment",
    }
    
    def __init__(self, preset: str = "asm20", k: Optional[int] = None, 
                 w: Optional[int] = None, scoring: Optional[Tuple[int, int, int, int]] = None):
        """
        Initialize the mappy aligner.
        
        Args:
            preset: Minimap2 preset (default: "asm20" for ~20% divergence)
            k: k-mer size (None = auto)
            w: Minimizer window size (None = auto)
            scoring: Tuple of (match, mismatch, gap_open, gap_extend) scores
        """
        if preset not in self.PRESETS and preset is not None:
            raise ValueError(f"Unknown preset '{preset}'. Available: {list(self.PRESETS.keys())}")
        
        self.preset = preset
        self.k = k
        self.w = w
        self.scoring = scoring
        self._aligner = None
    
    def index_reference(self, reference: SeqRecord) -> None:
        """
        Index the reference sequence for alignment.
        
        Args:
            reference: BioPython SeqRecord object
        """
        if reference.seq is None:
            raise ValueError("Reference sequence is empty")
        
        # Create mappy aligner with the reference sequence
        ref_seq = str(reference.seq)
        
        # Build kwargs for Aligner
        kwargs: dict[str, Any] = {}
        if self.preset:
            kwargs['preset'] = self.preset
        if self.k is not None:
            kwargs['k'] = self.k
        if self.w is not None:
            kwargs['w'] = self.w
        if self.scoring is not None:
            kwargs['scoring'] = self.scoring
        
        try:
            self._aligner = mappy.Aligner(seq=ref_seq, **kwargs)
        except Exception as e:
            raise UserError(f"Failed to create mappy aligner: {e}") from e
        
        if self._aligner is None:
            raise UserError("Failed to index reference sequence with mappy")
    
    def align(self, query: SeqRecord) -> List[MappyHit]:
        """
        Align a query sequence to the indexed reference.
        
        Args:
            query: BioPython SeqRecord object to align
        
        Returns:
            List of MappyHit objects (may be empty if no alignment found)
        
        Raises:
            ValueError: If reference has not been indexed
        """
        if self._aligner is None:
            raise ValueError("Reference must be indexed before alignment")
        
        if query.seq is None:
            return []
        
        query_seq = str(query.seq)
        
        try:
            hits = []
            for hit in self._aligner.map(query_seq):
                hits.append(MappyHit.from_mappy_alignment(hit))
            return hits
        except Exception as e:
            raise UserError(f"Mappy alignment failed: {e}") from e
    
    def align_best(self, query: SeqRecord) -> Optional[MappyHit]:
        """
        Align a query sequence and return the best hit.
        
        The "best" hit is determined by the number of matching bases.
        
        Args:
            query: BioPython SeqRecord object to align
        
        Returns:
            Best MappyHit or None if no alignment found
        """
        hits = self.align(query)
        if not hits:
            return None
        
        # Return hit with most matching bases
        return max(hits, key=lambda h: h.match_length)


def align_pair(reference: SeqRecord, query: SeqRecord, preset: str = "asm20") -> Optional[MappyHit]:
    """
    Convenience function to align two sequences.
    
    Args:
        reference: Reference sequence
        query: Query sequence
        preset: Minimap2 preset
    
    Returns:
        Best alignment hit or None if no alignment found
    
    Example:
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Seq import Seq
        >>> ref = SeqRecord(Seq("ACGTACGTACGT"), id="ref")
        >>> qry = SeqRecord(Seq("ACGTACGTACGT"), id="qry")
        >>> hit = align_pair(ref, qry)
        >>> print(hit.cigar_str)
    """
    aligner = MappyAligner(preset=preset)
    aligner.index_reference(reference)
    return aligner.align_best(query)
