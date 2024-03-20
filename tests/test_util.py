
import pytest

# from cfeintact.coordinates import map_positions, map_nonaligned_to_aligned_positions
from cfeintact.coordinates import map_nonaligned_to_aligned_positions
import cfeintact.subtypes as st

def map_positions_hxb2_style(reference, query):
    return [st.convert_from_aligned_to_reference(position, [reference, query])
            for position in range(len(reference) - reference.count("-"))]

# @pytest.mark.parametrize("reference, query, expected_mapping", [
#     ("ACGTTA", "ACGTTA", [0, 1, 2, 3, 4, 5]),  # No gaps or variations
#     ("A-GTTA", "ACGTTA", [0, 2, 3, 4, 5]),  # Gap 1 in reference
#     ("ACG-TA", "ACGTTA", [0, 1, 2, 4, 5]),  # Gap 1 in reference
#     ("ACGTTA", "ACG-TA", [0, 1, 2, 3, 3, 4]),  # Gap 1 in query
#     ("AC--TA", "ACGTTA", [0, 1, 4, 5]),  # Gap 2 in reference
#     ("ACTTTA", "AC--TA", [0, 1, 2, 2, 2, 3]),  # Gap 2 in query
#     ("AC-TAC", "A-TTAC", [0, 1, 2, 3, 4]),  # Gap 1 in both
#     ("AC--TA", "A--TTA", [0, 1, 2, 3]),  # Gap 2 in both
#     ("AC-TTA", "ACGT-A", [0, 1, 3, 4, 4]),  # Multiple gaps
#     ("AC-TTA", "AC-TTA", [0, 1, 2, 3, 4]),  # Gap 1 in both
#     ("AC-TTA", "A-GTTA", [0, 1, 2, 3, 4]),  # Multiple gap 1 in both
#     ("AC--TA", "A-G-TA", [0, 1, 2, 3]),  # Variation
#     ("ACGTT-", "ACGTTA", [0, 1, 2, 3, 4]),  # Reference ends in -
#     ("-CGTTA", "ACGTTA", [1, 2, 3, 4, 5]),  # Reference starts with -
#     ("ACGTTA", "ACGTT-", [0, 1, 2, 3, 4, 5]),  # Query ends in -
#     ("ACGTTA", "-CGTTA", [0, 0, 1, 2, 3, 4]),  # Query starts with -
#     ("-----A", "ACGTTA", [5]),  # Almost empty reference 2
#     ("A-----", "ACGTTA", [0]),  # Almost empty reference 1
#     ("ACGTTA", "-----A", [0, 0, 0, 0, 0, 0]),  # Almost empty query 1
#     ("ACGTTA", "A-----", [0, 1, 1, 1, 1, 1]),  # Almost empty query 2
#     ("ACGTTA", "------", [0, 0, 0, 0, 0, 0]),  # Empty query
#     ("------", "ACGTTA", []),  # Empty reference
#     ("------", "------", []),  # Empty both
# ])
# def test_map_positions(reference, query, expected_mapping):
#     assert map_positions(reference, query) == expected_mapping

#     if "-" not in reference:
#         assert map_positions_hxb2_style(reference, query) \
#             == map_positions(reference, query)


@pytest.mark.parametrize("reference, aligned, expected_mapping", [
    ("ACGTTA", "ACGTTA", [0, 1, 2, 3, 4, 5]),  # No gaps or variations
    ("ACGTTA", "AC-GTTA", [0, 1, 3, 4, 5, 6]),  # One gap
    ("ACGTTA", "ACG--TTA", [0, 1, 2, 5, 6, 7]),  # Two gaps
    ("ACGTTA", "ACGTT-A", [0, 1, 2, 3, 4, 6]),  # Gap at the end
    ("ACGTTA", "-ACGTTA", [1, 2, 3, 4, 5, 6]),  # Gap at the beginning
    ("ACGTTA", "-AC-GT-TA", [1, 2, 4, 5, 7, 8]),  # Multiple gaps
])
def test_map_nonaligned_to_aligned_positions(reference, aligned, expected_mapping):
    assert map_nonaligned_to_aligned_positions(reference, aligned) == expected_mapping
