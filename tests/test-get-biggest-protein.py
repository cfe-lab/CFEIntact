
import pytest

from cfeintact.intact import get_biggest_protein

@pytest.mark.parametrize("aminoseq, expected_output", [
    ("MXX*MXXX*MX*", "MXXX"),
    ("MXX*MXXX*MX", "MXXX"),
    ("MXX*MXXX*MXXXXX", "MXXXXX"),
    ("MXX**MX*", "MXX"),
    ("MXX***MX***", "MXX"),
    ("MXXX", "MXXX"),
    ("XMXX*MX*", "MXX"),
    ("XXXXXXMXX*MX*", "MXX"),
    ("XMXX*XMX*", "MXX"),
    ("XMXX*XXXXXMXXX*", "MXXX"),
    ("XMXXX*XXXXXXMX*", "MXXX"),
    ("This is a test", ""),
    ("", ""),
])
def test_find_biggest_protein_with_start_codon(aminoseq, expected_output):
    assert get_biggest_protein(True, aminoseq) == expected_output

@pytest.mark.parametrize("aminoseq, expected_output", [
    ("MXX*MXXX*MX*", "MXXX"),
    ("MXX*MXXX*MX", "MXXX"),
    ("MXX*MXXX*MXXXXX", "MXXXXX"),
    ("MXX**MX*", "MXX"),
    ("MXX***MX***", "MXX"),
    ("MXXX", "MXXX"),
    ("XMXX*MX*", "XMXX"),
    ("XXXXXXMXX*MX*", "XXXXXXMXX"),
    ("XMXX*XMX*", "XMXX"),
    ("XMXX*XXXXXMXXX*", "XXXXXMXXX"),
    ("XMXXX*XXXXXXMX*", "XXXXXXMX"),
    ("This is a test", "This is a test"),
    ("", ""),
])
def test_find_biggest_protein_without_start_codon(aminoseq, expected_output):
    assert get_biggest_protein(False, aminoseq) == expected_output
