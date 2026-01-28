"""
Unit tests for MappyAlignedSequence class.
"""

import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from cfeintact.mappy_aligned_sequence import MappyAlignedSequence
from cfeintact.aligned_sequence import AlignedSequence


def create_test_sequences():
    """Create test sequences for alignment."""
    ref = SeqRecord(
        Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        id="reference",
        name="ref"
    )
    query = SeqRecord(
        Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        id="query",
        name="qry"
    )
    return ref, query


def test_mappy_aligned_sequence_basic():
    """Test basic MappyAlignedSequence functionality."""
    ref, query = create_test_sequences()
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    
    # Test alignment property
    alignment = aligned.alignment
    assert len(alignment) == 2
    assert alignment[0].id == "reference"
    assert alignment[1].id == "query"


def test_mappy_aligned_sequence_properties():
    """Test aligned sequence properties."""
    ref, query = create_test_sequences()
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    
    # Test aligned_reference property
    aligned_ref = aligned.aligned_reference
    assert aligned_ref.id == "reference"
    
    # Test aligned_this property
    aligned_qry = aligned.aligned_this
    assert aligned_qry.id == "query"


def test_mappy_aligned_sequence_coordinate_mapping():
    """Test coordinate mapping functionality."""
    ref, query = create_test_sequences()
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    
    # Test coordinate mapping
    coord_mapping = aligned.coordinate_mapping
    assert coord_mapping is not None
    assert hasattr(coord_mapping, 'ref_to_query')
    assert hasattr(coord_mapping, 'query_to_ref')


def test_mappy_aligned_sequence_reverse():
    """Test reverse complement functionality."""
    ref, query = create_test_sequences()
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    reversed_aligned = aligned.reverse()
    
    # Check that it returns a new MappyAlignedSequence
    assert isinstance(reversed_aligned, MappyAlignedSequence)
    assert reversed_aligned.reference.id == ref.id
    
    # Check that query is reverse complemented
    original_seq = str(query.seq)
    reversed_seq = str(reversed_aligned.this.seq)
    assert reversed_seq == str(Seq(original_seq).reverse_complement())


def test_mappy_aligned_sequence_alignment_score():
    """Test alignment scoring."""
    ref, query = create_test_sequences()
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    score = aligned.alignment_score()
    
    # For identical sequences, score should be equal to sequence length
    assert score > 0
    assert score == len(ref.seq)


def test_mappy_aligned_sequence_with_mismatch():
    """Test alignment with mismatches."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query = SeqRecord(Seq("ACGTCCGTACGTACGT"), id="qry")  # One mismatch
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    score = aligned.alignment_score()
    
    # Score should be less than perfect match
    assert score < len(ref.seq)
    assert score > len(ref.seq) - 5  # Most positions still match


def test_mappy_aligned_sequence_with_gap():
    """Test alignment with gaps (insertions/deletions)."""
    ref = SeqRecord(Seq("ACGTACGTACGTACGT"), id="ref")
    query = SeqRecord(Seq("ACGTACGTACGT"), id="qry")  # 4 bp deletion
    
    aligned = MappyAlignedSequence(this=query, reference=ref)
    alignment = aligned.alignment
    
    # Check that alignment is created
    assert len(alignment) == 2
    
    # One of the sequences should have gaps
    ref_aligned_str = str(alignment[0].seq)
    qry_aligned_str = str(alignment[1].seq)
    assert '-' in ref_aligned_str or '-' in qry_aligned_str


def test_mappy_vs_mafft_compatibility():
    """Test that MappyAlignedSequence has same interface as AlignedSequence."""
    ref, query = create_test_sequences()
    
    # Both should have the same methods and properties
    mappy_aligned = MappyAlignedSequence(this=query, reference=ref)
    
    # Check all required attributes exist
    assert hasattr(mappy_aligned, 'alignment')
    assert hasattr(mappy_aligned, 'aligned_reference')
    assert hasattr(mappy_aligned, 'aligned_this')
    assert hasattr(mappy_aligned, 'coordinate_mapping')
    assert hasattr(mappy_aligned, 'reverse')
    assert hasattr(mappy_aligned, 'alignment_score')
    
    # Check they're callable/accessible
    assert mappy_aligned.alignment is not None
    assert mappy_aligned.aligned_reference is not None
    assert mappy_aligned.aligned_this is not None
    assert mappy_aligned.coordinate_mapping is not None
    assert callable(mappy_aligned.reverse)
    assert callable(mappy_aligned.alignment_score)


def test_different_presets():
    """Test different minimap2 presets."""
    ref, query = create_test_sequences()
    
    for preset in ["asm5", "asm10", "asm20"]:
        aligned = MappyAlignedSequence(this=query, reference=ref, preset=preset)
        alignment = aligned.alignment
        
        assert len(alignment) == 2
        assert aligned.alignment_score() > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
