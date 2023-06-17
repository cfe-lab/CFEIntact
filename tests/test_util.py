
import pytest

from util.coordinates import map_positions
import util.subtypes as st

def map_positions_hxb2_style(reference, query):
    return [st.convert_from_aligned_to_reference(position, [reference, query])
            for position in range(len(reference) - reference.count("-"))]

@pytest.mark.parametrize("reference, query, expected_mapping", [
    ("ACGTTA", "ACGTTA", [0, 1, 2, 3, 4, 5]),  # No gaps or variations
    ("A-GTTA", "ACGTTA", [0, 2, 3, 4, 5]),  # Gap 1 in reference
    ("ACG-TA", "ACGTTA", [0, 1, 2, 4, 5]),  # Gap 1 in reference
    ("ACGTTA", "ACG-TA", [0, 1, 2, 3, 3, 4]),  # Gap 1 in query
    ("AC--TA", "ACGTTA", [0, 1, 4, 5]),  # Gap 2 in reference
    ("ACTTTA", "AC--TA", [0, 1, 2, 2, 2, 3]),  # Gap 2 in query
    ("AC--TA", "A--TTA", [0, 1, 2, 3]),  # Gap 2 in both
    ("AC-TTA", "ACGT-A", [0, 1, 3, 4, 4]),  # Multiple gaps
    ("AC-TTA", "AC-TTA", [0, 1, 2, 3, 4]),  # Gap 1 in both
    ("AC-TTA", "A-GTTA", [0, 1, 2, 3, 4]),  # Multiple gap 1 in both
    ("AC--TA", "A-G-TA", [0, 1, 2, 3]),  # Variation
])
def test_map_positions(reference, query, expected_mapping):
    assert map_positions(reference, query) == expected_mapping

    if "-" not in reference:
        assert map_positions_hxb2_style(reference, query) \
            == map_positions(reference, query)
