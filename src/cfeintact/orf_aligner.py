"""
Local alignment of ORFs using minimap2/mappy.

This module provides functionality to align individual Open Reading Frames (ORFs)
from a reference sequence to a query sequence using minimap2 through its Python
bindings (mappy). This approach is more efficient than whole-sequence alignment
for ORF detection.
"""

import mappy  # type: ignore
from Bio.SeqRecord import SeqRecord
from aligntools.cigar_hit import CigarHit
from aligntools.cigar import Cigar
from typing import Optional, List
from dataclasses import dataclass

from cfeintact.original_orf import OriginalORF
from cfeintact.user_error import UserError


@dataclass
class OrfAlignment:
    """
    Result of aligning an ORF between reference and query sequences.
    
    Attributes:
        cigar_hit: aligntools CigarHit object with coordinate mapping
        reference_orf: Original ORF definition from reference
        query_start: Start position in full query sequence (0-based)
        query_end: End position in full query sequence (0-based, inclusive)
        mapq: Mapping quality from minimap2
        strand: 1 for forward, -1 for reverse
        match_length: Number of matching bases
    """
    cigar_hit: CigarHit
    reference_orf: OriginalORF
    query_start: int
    query_end: int
    mapq: int
    strand: int
    match_length: int


class OrfAligner:
    """
    Aligns ORFs using minimap2 local alignment.
    
    This class uses minimap2's strength as a local aligner to find where
    specific ORFs map from a reference sequence to a query sequence.
    """
    
    def __init__(self, preset: str = "asm20"):
        """
        Initialize the ORF aligner.
        
        Args:
            preset: Minimap2 preset. Options:
                - "asm5": ~5% divergence (stricter)
                - "asm10": ~10% divergence
                - "asm20": ~20% divergence (default, more permissive)
        """
        self.preset = preset
    
    def align_orf(self, 
                  reference: SeqRecord,
                  query: SeqRecord,
                  orf: OriginalORF) -> Optional[OrfAlignment]:
        """
        Align a specific ORF region from reference to query using local alignment.
        
        This performs local alignment of the ORF sequence, allowing it to find
        where the ORF maps in the query even if there are large insertions/
        deletions elsewhere in the sequence.
        
        Args:
            reference: Reference sequence (e.g., HXB2 or subtype)
            query: Query sequence to analyze
            orf: ORF definition with start/end positions in reference
        
        Returns:
            OrfAlignment if alignment found, None otherwise
        
        Example:
            >>> aligner = OrfAligner()
            >>> alignment = aligner.align_orf(hxb2, query, gag_orf)
            >>> if alignment:
            ...     print(f"Found at query position {alignment.query_start}")
            ...     mapping = alignment.cigar_hit.coordinate_mapping
        """
        if reference.seq is None or query.seq is None:
            return None
        
        # Extract ORF sequence from reference
        ref_orf_seq = str(reference.seq[orf.start:orf.end + 1])
        query_seq = str(query.seq)
        
        if not ref_orf_seq:
            raise UserError(f"ORF {orf.name} has empty reference sequence")
        
        try:
            # Create aligner for this specific ORF
            aligner = mappy.Aligner(seq=ref_orf_seq, preset=self.preset)
            
            if aligner is None:
                raise UserError(f"Failed to create mappy aligner for ORF {orf.name}")
            
            # Find local alignments
            hits = list(aligner.map(query_seq))
            
            if not hits:
                return None
            
            # Get best hit (most matching bases)
            best_hit = max(hits, key=lambda h: h.mlen)
            
            # Convert minimap2 CIGAR to aligntools Cigar
            cigar = Cigar.coerce(best_hit.cigar_str)
            
            # Create aligntools CigarHit
            # Note: minimap2 uses exclusive end (r_en, q_en)
            # CigarHit uses inclusive end (r_ei, q_ei)
            cigar_hit = CigarHit(
                cigar=cigar,
                r_st=best_hit.r_st,
                r_ei=best_hit.r_en - 1,
                q_st=best_hit.q_st,
                q_ei=best_hit.q_en - 1
            )
            
            return OrfAlignment(
                cigar_hit=cigar_hit,
                reference_orf=orf,
                query_start=best_hit.q_st,
                query_end=best_hit.q_en - 1,
                mapq=best_hit.mapq,
                strand=best_hit.strand,
                match_length=best_hit.mlen
            )
        
        except Exception as e:
            raise UserError(f"Failed to align ORF {orf.name}: {e}") from e
    
    def align_all_orfs(self,
                       reference: SeqRecord,
                       query: SeqRecord,
                       orfs: List[OriginalORF]) -> List[Optional[OrfAlignment]]:
        """
        Align multiple ORFs from reference to query.
        
        Args:
            reference: Reference sequence
            query: Query sequence
            orfs: List of ORF definitions
        
        Returns:
            List of OrfAlignment (or None if no alignment found for that ORF)
        """
        results = []
        for orf in orfs:
            alignment = self.align_orf(reference, query, orf)
            results.append(alignment)
        return results
