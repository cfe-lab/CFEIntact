"""
CIGAR string parsing utilities for converting mappy alignment output
to coordinate mappings compatible with the rest of CFEIntact.
"""

import re
from typing import List, Tuple, Optional


def parse_cigar(cigar_str: str) -> List[Tuple[int, str]]:
    """
    Parse a CIGAR string into a list of (length, operation) tuples.
    
    Args:
        cigar_str: CIGAR string (e.g., "100M2D50M3I75M")
    
    Returns:
        List of (length, operation) tuples
        
    Example:
        >>> parse_cigar("100M2D50M")
        [(100, 'M'), (2, 'D'), (50, 'M')]
    """
    if not cigar_str:
        return []
    
    pattern = r'(\d+)([MIDNSHP=X])'
    matches = re.findall(pattern, cigar_str)
    return [(int(length), op) for length, op in matches]


class CigarCoordinateMapper:
    """
    Maps coordinates between reference and query sequences using CIGAR operations.
    
    This class provides coordinate translation similar to aligntools' CoordinateMapping
    but derived from mappy's CIGAR output.
    """
    
    def __init__(self, cigar_str: str, ref_start: int = 0, query_start: int = 0):
        """
        Initialize the coordinate mapper.
        
        Args:
            cigar_str: CIGAR string from mappy alignment
            ref_start: Starting position in reference (0-based)
            query_start: Starting position in query (0-based)
        """
        self.cigar = parse_cigar(cigar_str)
        self.ref_start = ref_start
        self.query_start = query_start
        self._build_mappings()
    
    def _build_mappings(self):
        """
        Build forward and reverse coordinate mapping tables.
        
        This creates lookup tables for efficient position translation.
        """
        self.ref_to_query_map: dict[int, Optional[int]] = {}
        self.query_to_ref_map: dict[int, Optional[int]] = {}
        
        ref_pos = self.ref_start
        query_pos = self.query_start
        
        for length, op in self.cigar:
            if op in 'M=X':  # Match or mismatch
                for i in range(length):
                    self.ref_to_query_map[ref_pos + i] = query_pos + i
                    self.query_to_ref_map[query_pos + i] = ref_pos + i
                ref_pos += length
                query_pos += length
            elif op == 'D':  # Deletion from reference (gap in query)
                # Reference advances, query does not
                for i in range(length):
                    self.ref_to_query_map[ref_pos + i] = None
                ref_pos += length
            elif op == 'I':  # Insertion to reference (gap in reference)
                # Query advances, reference does not
                for i in range(length):
                    self.query_to_ref_map[query_pos + i] = None
                query_pos += length
            elif op == 'N':  # Skipped region (for RNA-seq)
                ref_pos += length
            elif op == 'S':  # Soft clipping
                query_pos += length
            elif op == 'H':  # Hard clipping
                # No position change
                pass
    
    def ref_to_query_position(self, ref_pos: int) -> Optional[int]:
        """
        Map a reference position to query position.
        
        Args:
            ref_pos: Position in reference sequence (0-based)
        
        Returns:
            Corresponding query position, or None if position is in a deletion
        """
        return self.ref_to_query_map.get(ref_pos)
    
    def query_to_ref_position(self, query_pos: int) -> Optional[int]:
        """
        Map a query position to reference position.
        
        Args:
            query_pos: Position in query sequence (0-based)
        
        Returns:
            Corresponding reference position, or None if position is in an insertion
        """
        return self.query_to_ref_map.get(query_pos)
    
    def right_min(self, ref_pos: int) -> Optional[int]:
        """
        Find the minimum query position that maps to or after the reference position.
        
        This handles cases where the exact reference position is in a deletion.
        
        Args:
            ref_pos: Reference position
        
        Returns:
            Minimum query position >= mapped position, or None if not found
        """
        # First try exact mapping
        query_pos = self.ref_to_query_position(ref_pos)
        if query_pos is not None:
            return query_pos
        
        # If exact position is in deletion, scan forward to find next match
        for r_pos in range(ref_pos, max(self.ref_to_query_map.keys()) + 1):
            q_pos = self.ref_to_query_map.get(r_pos)
            if q_pos is not None:
                return q_pos
        
        return None
    
    def left_max(self, ref_pos: int) -> Optional[int]:
        """
        Find the maximum query position that maps to or before the reference position.
        
        This handles cases where the exact reference position is in a deletion.
        
        Args:
            ref_pos: Reference position
        
        Returns:
            Maximum query position <= mapped position, or None if not found
        """
        # First try exact mapping
        query_pos = self.ref_to_query_position(ref_pos)
        if query_pos is not None:
            return query_pos
        
        # If exact position is in deletion, scan backward to find previous match
        for r_pos in range(ref_pos, min(self.ref_to_query_map.keys()) - 1, -1):
            q_pos = self.ref_to_query_map.get(r_pos)
            if q_pos is not None:
                return q_pos
        
        return None


class AlignToolsCompatibleMapping:
    """
    Adapter to make CigarCoordinateMapper compatible with aligntools' CoordinateMapping interface.
    
    This allows mappy-based alignments to work with existing CFEIntact code that expects
    aligntools coordinate mappings.
    """
    
    def __init__(self, mapper: CigarCoordinateMapper):
        """
        Initialize with a CIGAR-based coordinate mapper.
        
        Args:
            mapper: CigarCoordinateMapper instance
        """
        self.mapper = mapper
        self.ref_to_query = RefToQueryMapper(mapper)
        self.query_to_ref = QueryToRefMapper(mapper)


class RefToQueryMapper:
    """Adapter for reference-to-query mapping."""
    
    def __init__(self, mapper: CigarCoordinateMapper):
        self.mapper = mapper
    
    def right_min(self, ref_pos: int) -> Optional[int]:
        return self.mapper.right_min(ref_pos)
    
    def left_max(self, ref_pos: int) -> Optional[int]:
        return self.mapper.left_max(ref_pos)
    
    def __call__(self, ref_pos: int) -> Optional[int]:
        return self.mapper.ref_to_query_position(ref_pos)


class QueryToRefMapper:
    """Adapter for query-to-reference mapping."""
    
    def __init__(self, mapper: CigarCoordinateMapper):
        self.mapper = mapper
    
    def __call__(self, query_pos: int) -> Optional[int]:
        return self.mapper.query_to_ref_position(query_pos)
