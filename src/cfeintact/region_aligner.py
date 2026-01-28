"""
Global alignment of small regions using BioPython PairwiseAligner.

This module provides functionality to align small critical regions (like MSD, PSI, RRE)
using BioPython's PairwiseAligner. These regions are too small for local alignment
and require precise global alignment.
"""

from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass
from typing import Optional


@dataclass
class RegionAlignment:
    """
    Result of aligning a small region globally.
    
    Attributes:
        reference_aligned: Aligned reference sequence (with gaps)
        query_aligned: Aligned query sequence (with gaps)
        score: Alignment score
        query_start: Start position in query where alignment was searched
        query_end: End position in query where alignment was searched
    """
    reference_aligned: str
    query_aligned: str
    score: float
    query_start: int
    query_end: int


class RegionAligner:
    """
    Aligns small critical regions using BioPython PairwiseAligner.
    
    This class performs global alignment on small regions where precise
    positional information is needed (e.g., Major Splice Donor site).
    """
    
    def __init__(self):
        """
        Initialize the region aligner with default scoring parameters.
        
        Scoring scheme:
        - Match: +2
        - Mismatch: -1
        - Gap open: -1.5
        - Gap extend: -0.2
        """
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -1.5
        self.aligner.extend_gap_score = -0.2
    
    def align_region(self,
                    reference_seq: str,
                    query_seq: str,
                    query_search_start: int = 0,
                    query_search_end: Optional[int] = None) -> RegionAlignment:
        """
        Globally align a small region.
        
        Args:
            reference_seq: Small reference sequence (e.g., "GT" for MSD)
            query_seq: Full query sequence or substring
            query_search_start: Where to start looking in query (0-based)
            query_search_end: Where to stop looking in query (exclusive, None = end)
        
        Returns:
            RegionAlignment with aligned sequences
        
        Example:
            >>> aligner = RegionAligner()
            >>> # Look for MSD site near position 740
            >>> result = aligner.align_region("GT", query_seq, 700, 800)
            >>> if result.query_aligned.upper().replace("-", "") != "GT":
            ...     print("MSD is mutated!")
        """
        if query_search_end is None:
            query_search_end = len(query_seq)
        
        # Extract search region from query
        search_region = query_seq[query_search_start:query_search_end]
        
        if not reference_seq or not search_region:
            # Return empty alignment
            return RegionAlignment(
                reference_aligned=reference_seq,
                query_aligned="-" * len(reference_seq),
                score=0.0,
                query_start=query_search_start,
                query_end=query_search_start
            )
        
        # Perform global alignment
        try:
            alignments = list(self.aligner.align(reference_seq, search_region))
            
            if not alignments:
                # No alignment found
                return RegionAlignment(
                    reference_aligned=reference_seq,
                    query_aligned="-" * len(reference_seq),
                    score=0.0,
                    query_start=query_search_start,
                    query_end=query_search_end
                )
            
            # Get best alignment
            best = alignments[0]
            
            return RegionAlignment(
                reference_aligned=str(best[0]),
                query_aligned=str(best[1]),
                score=float(best.score),
                query_start=query_search_start,
                query_end=query_search_end
            )
        
        except Exception:
            # If alignment fails, return empty result
            return RegionAlignment(
                reference_aligned=reference_seq,
                query_aligned="-" * len(reference_seq),
                score=0.0,
                query_start=query_search_start,
                query_end=query_search_end
            )
    
    def find_in_sequence(self,
                        pattern: str,
                        sequence: str,
                        context_before: int = 50,
                        context_after: int = 50) -> Optional[RegionAlignment]:
        """
        Find a pattern in a sequence with context.
        
        This is useful for finding small motifs like MSD site where we
        want to search in a window around the expected position.
        
        Args:
            pattern: Pattern to find (e.g., "GT" for MSD)
            sequence: Sequence to search in
            context_before: Bases before pattern to include
            context_after: Bases after pattern to include
        
        Returns:
            RegionAlignment if found, None if not found with high confidence
        """
        # Simple substring search first
        pos = sequence.upper().find(pattern.upper())
        
        if pos >= 0:
            # Found exact match
            start = max(0, pos - context_before)
            end = min(len(sequence), pos + len(pattern) + context_after)
            
            # Align the region to confirm
            result = self.align_region(pattern, sequence, start, end)
            return result
        
        # If no exact match, could do fuzzy search here
        # For now, return None
        return None
