"""
Unit tests for CIGAR parsing and coordinate mapping.
"""

import pytest
from cfeintact.cigar_parser import (
    parse_cigar,
    CigarCoordinateMapper,
    AlignToolsCompatibleMapping
)


def test_parse_cigar_basic():
    """Test basic CIGAR string parsing."""
    cigar = "100M2D50M3I75M"
    ops = parse_cigar(cigar)
    
    assert len(ops) == 5
    assert ops[0] == (100, 'M')
    assert ops[1] == (2, 'D')
    assert ops[2] == (50, 'M')
    assert ops[3] == (3, 'I')
    assert ops[4] == (75, 'M')


def test_parse_cigar_simple():
    """Test simple CIGAR string."""
    cigar = "100M"
    ops = parse_cigar(cigar)
    
    assert len(ops) == 1
    assert ops[0] == (100, 'M')


def test_parse_cigar_empty():
    """Test empty CIGAR string."""
    cigar = ""
    ops = parse_cigar(cigar)
    
    assert len(ops) == 0


def test_coordinate_mapper_match_only():
    """Test coordinate mapping with match operations only."""
    cigar = "100M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    
    # Direct position mapping
    assert mapper.ref_to_query_position(0) == 0
    assert mapper.ref_to_query_position(50) == 50
    assert mapper.ref_to_query_position(99) == 99
    
    # Reverse mapping
    assert mapper.query_to_ref_position(0) == 0
    assert mapper.query_to_ref_position(50) == 50
    assert mapper.query_to_ref_position(99) == 99


def test_coordinate_mapper_with_deletion():
    """Test coordinate mapping with deletion."""
    # 50M2D50M means:
    # - First 50 positions: direct mapping
    # - Positions 50-51 in ref: deletion (no query mapping)
    # - Next 50 positions: ref 52-101 -> query 50-99
    cigar = "50M2D50M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    
    # Before deletion
    assert mapper.ref_to_query_position(0) == 0
    assert mapper.ref_to_query_position(49) == 49
    
    # In deletion
    assert mapper.ref_to_query_position(50) is None
    assert mapper.ref_to_query_position(51) is None
    
    # After deletion
    assert mapper.ref_to_query_position(52) == 50
    assert mapper.ref_to_query_position(101) == 99


def test_coordinate_mapper_with_insertion():
    """Test coordinate mapping with insertion."""
    # 50M3I50M means:
    # - First 50 positions: direct mapping
    # - Positions 50-52 in query: insertion (no ref mapping)
    # - Next 50 positions: ref 50-99 -> query 53-102
    cigar = "50M3I50M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    
    # Before insertion
    assert mapper.ref_to_query_position(0) == 0
    assert mapper.ref_to_query_position(49) == 49
    
    # After insertion (query positions 50-52 are insertion)
    assert mapper.ref_to_query_position(50) == 53
    assert mapper.ref_to_query_position(99) == 102
    
    # Reverse mapping in insertion
    assert mapper.query_to_ref_position(50) is None
    assert mapper.query_to_ref_position(51) is None
    assert mapper.query_to_ref_position(52) is None


def test_coordinate_mapper_right_min():
    """Test right_min method for finding next valid position."""
    cigar = "50M2D50M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    
    # Normal position
    assert mapper.right_min(0) == 0
    
    # Position in deletion - should return next valid query position
    assert mapper.right_min(50) == 50
    assert mapper.right_min(51) == 50


def test_coordinate_mapper_left_max():
    """Test left_max method for finding previous valid position."""
    cigar = "50M2D50M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    
    # Normal position
    assert mapper.left_max(0) == 0
    
    # Position in deletion - should return previous valid query position
    assert mapper.left_max(50) == 49
    assert mapper.left_max(51) == 49


def test_aligntools_compatible_mapping():
    """Test compatibility wrapper for aligntools."""
    cigar = "100M2D50M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    compat = AlignToolsCompatibleMapping(mapper)
    
    # Test ref_to_query interface
    assert compat.ref_to_query(0) == 0
    assert compat.ref_to_query.right_min(100) == 98
    assert compat.ref_to_query.left_max(100) == 98
    
    # Test query_to_ref interface
    assert compat.query_to_ref(0) == 0


def test_coordinate_mapper_with_offset():
    """Test coordinate mapping with non-zero start positions."""
    cigar = "100M"
    mapper = CigarCoordinateMapper(cigar, ref_start=100, query_start=200)
    
    assert mapper.ref_to_query_position(100) == 200
    assert mapper.ref_to_query_position(150) == 250
    assert mapper.query_to_ref_position(200) == 100


def test_complex_cigar():
    """Test complex CIGAR string with multiple operations."""
    cigar = "10M2I20M5D30M3I15M"
    mapper = CigarCoordinateMapper(cigar, ref_start=0, query_start=0)
    
    # Verify basic functionality works with complex CIGAR
    assert mapper.ref_to_query_position(0) == 0
    assert mapper.ref_to_query_position(9) == 9
    # After 10M2I: ref 10 -> query 12
    assert mapper.ref_to_query_position(10) == 12
    

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
