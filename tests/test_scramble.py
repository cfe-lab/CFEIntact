
import pytest
import os

import intact.intact as intact
from intact.intact import is_scrambled, is_nonhiv, IntactnessError

@pytest.mark.parametrize("lst, expected", [
    ([1, 2, 3, 4, 5], True),
    ([1, 3, 2, 4, 5], False),
    ([], True),
    ([5], True),
    ([1, 1, 2, 2, 3, 3], True),
])
def test_is_sorted(lst, expected):
    assert intact.is_sorted(lst) == expected

# def test_blast(tmp_path):
#     input_file = "tests/data-large.fasta"
#     output_file = os.path.join(tmp_path, "out.tsv")

#     # intact.blast(
#     #     subtype="B",
#     #     input_file=input_file,
#     #     output_file=output_file,
#     # )

#     it = intact.blast_interate(subtype="B", input_file=input_file)
#     for v in it:
#         print("\n\n\n")
#         print("---------------------------")
#         print(v[0].qseqid)
#         print(intact.is_scrambled(v))

#     assert False

class BlastRow:
    def __init__(self, sstart, send, qstart, sstrand):
        self.sstart = sstart
        self.send = send
        self.qstart = qstart
        self.sstrand = sstrand

def test_is_scrambled_no_alignment():
    # Test case where there is no alignment
    blast_rows = []
    assert is_scrambled("id", blast_rows) is None

def test_is_scrambled_internal_inversion():
    # Test case where some parts of the sequence are aligned in forward direction
    # and some in reverse, indicating an internal inversion error.
    blast_rows = [
        BlastRow(sstart=900, send=910, qstart=5, sstrand="plus"),
        BlastRow(sstart=930, send=940, qstart=25, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.INTERNALINVERSION_ERROR

def test_is_scrambled_plus_strand_sorted():
    # Test case where the direction is "plus" and the sstart values are sorted.
    blast_rows = [
        BlastRow(sstart=900, send=910, qstart=5, sstrand="plus"),
        BlastRow(sstart=930, send=940, qstart=25, sstrand="plus"),
    ]
    assert is_scrambled("id", blast_rows) is None

def test_is_scrambled_minus_strand_sorted():
    # Test case where the direction is "minus" and the send values are sorted in reverse order.
    blast_rows = [
        BlastRow(sstart=910, send=900, qstart=25, sstrand="minus"),
        BlastRow(sstart=940, send=930, qstart=5, sstrand="minus"),
    ]
    assert is_scrambled("id", blast_rows) is None

def test_is_scrambled_plus_strand_unsorted():
    # Test case with mixed directions and inversions
    blast_rows = [
        BlastRow(sstart=700, send=700, qstart=99, sstrand="plus"),
        BlastRow(sstart=880, send=890, qstart=25, sstrand="plus"),
        BlastRow(sstart=920, send=930, qstart=45, sstrand="plus"),
        BlastRow(sstart=950, send=980, qstart=85, sstrand="plus"),
    ]
    result = is_scrambled("id", blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.SCRAMBLE_ERROR

def test_is_scrambled_plus_strand_unsorted_5prime():
    # Test case with mixed directions and inversions
    blast_rows = [
        BlastRow(sstart=100, send=110, qstart=99, sstrand="plus"),
        BlastRow(sstart=880, send=890, qstart=25, sstrand="plus"),
        BlastRow(sstart=920, send=930, qstart=45, sstrand="plus"),
        BlastRow(sstart=950, send=980, qstart=85, sstrand="plus"),
    ]
    assert is_scrambled("id", blast_rows) == None

def test_is_scrambled_minus_strand_unsorted():
    # Test case where the direction is "minus" but the send values are not sorted in reverse order.
    blast_rows = [
        BlastRow(sstart=910, send=900, qstart=5, sstrand="minus"),
        BlastRow(sstart=940, send=930, qstart=25, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.SCRAMBLE_ERROR

def test_is_scrambled_mixed_direction():
    # Test case where some rows have "plus" direction and some have "minus" direction.
    blast_rows = [
        BlastRow(sstart=900, send=910, qstart=5, sstrand="plus"),
        BlastRow(sstart=930, send=940, qstart=25, sstrand="minus"),
        BlastRow(sstart=950, send=990, qstart=45, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.INTERNALINVERSION_ERROR

def test_is_scrambled_single_row_plus_strand():
    # Test case with a single row aligned in the "plus" strand
    blast_rows = [
        BlastRow(sstart=700, send=710, qstart=5, sstrand="plus"),
    ]
    assert is_scrambled("id", blast_rows) is None

def test_is_scrambled_single_row_minus_strand():
    # Test case with a single row aligned in the "minus" strand
    blast_rows = [
        BlastRow(sstart=900, send=910, qstart=5, sstrand="minus"),
    ]
    assert is_scrambled("id", blast_rows) is None

def test_is_scrambled_multiple_directions_and_inversions():
    # Test case with mixed directions and inversions
    blast_rows = [
        BlastRow(sstart=900, send=910, qstart=5, sstrand="plus"),
        BlastRow(sstart=880, send=890, qstart=25, sstrand="minus"),
        BlastRow(sstart=920, send=930, qstart=45, sstrand="minus"),
        BlastRow(sstart=850, send=880, qstart=85, sstrand="plus"),
        BlastRow(sstart=890, send=880, qstart=85, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.INTERNALINVERSION_ERROR
