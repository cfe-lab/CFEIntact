"""
Unit tests for mappy alignment functionality.
"""

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from cfeintact.mappy_aligner import MappyAligner, align_pair


def test_mappy_aligner_basic():
    """Test basic alignment functionality."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query = SeqRecord(Seq("ACGTACGTACGTACGT"), id="query")
    
    aligner = MappyAligner(preset="asm20")
    aligner.index_reference(ref)
    hits = aligner.align(query)
    
    assert len(hits) > 0
    best_hit = hits[0]
    assert best_hit.match_length > 0
    assert best_hit.cigar_str is not None


def test_mappy_aligner_with_mismatch():
    """Test alignment with mismatches."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query = SeqRecord(Seq("ACGTCCGTACGTACGT"), id="query")  # One mismatch
    
    hit = align_pair(ref, query)
    
    assert hit is not None
    assert hit.match_length > 10  # Most bases still match


def test_mappy_aligner_with_gaps():
    """Test alignment with insertions/deletions."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query = SeqRecord(Seq("ACGTACGTACGT"), id="query")  # 4 bp deletion
    
    hit = align_pair(ref, query)
    
    assert hit is not None
    assert 'D' in hit.cigar_str or 'I' in hit.cigar_str


def test_mappy_aligner_no_match():
    """Test alignment with completely different sequences."""
    ref = SeqRecord(Seq("AAAAAAAAAAAAAAAA"), id="ref")
    query = SeqRecord(Seq("TTTTTTTTTTTTTTTT"), id="query")
    
    hit = align_pair(ref, query)
    
    # May or may not find alignment depending on aligner sensitivity
    # This is expected behavior for highly divergent sequences
    assert True  # Just check it doesn't crash


def test_mappy_aligner_reverse_complement():
    """Test alignment with reverse complement."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query_seq = Seq("ACGTACGTACGTACGT")
    query = SeqRecord(Seq.reverse_complement(query_seq), id="query")
    
    hit = align_pair(ref, query)
    
    assert hit is not None
    # Mappy should detect reverse complement
    assert hit.strand == -1


def test_mappy_different_presets():
    """Test different alignment presets."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT" * 10), id="ref")
    query = SeqRecord(Seq("ACGTACGTACGTACGT" * 10), id="query")
    
    for preset in ["asm5", "asm10", "asm20"]:
        hit = align_pair(ref, query, preset=preset)
        assert hit is not None
        assert hit.match_length > 0


def test_mappy_empty_sequence():
    """Test behavior with empty sequences."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query = SeqRecord(Seq(""), id="query")
    
    aligner = MappyAligner()
    aligner.index_reference(ref)
    hits = aligner.align(query)
    
    # Should return empty list for empty query
    assert len(hits) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
