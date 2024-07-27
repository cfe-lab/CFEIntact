
import pytest

from cfeintact.get_biggest_protein import get_biggest_protein

@pytest.mark.parametrize("aminoseq, expected_output", [
    ("MACTG*", "MACTG"),
    ("ACTG*", None),  # No start codon.
    ("MACTG", "MACTG"),  # No stop codon is fine.
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
    ("This is a test", None),
    ("", None),
])
def test_find_biggest_protein_with_start_codon(aminoseq, expected_output):
    result = get_biggest_protein(True, aminoseq)
    if expected_output is not None:
        assert result is not None
        protein, start, end = result
        assert protein == expected_output
        assert aminoseq[start:end + 1] == protein
    else:
        assert result is None


@pytest.mark.parametrize("aminoseq, expected_output", [
    ("MACTG*", "MACTG"),
    ("ACTG*", "ACTG"),
    ("MACTG", "MACTG"),  # No stop codon.
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
    ("ACTG*", "ACTG"),
    ("MACTG*", "MACTG"),
    ("This is a test", "This is a test"),
    ("", None),
])
def test_find_biggest_protein_without_start_codon(aminoseq, expected_output):
    result = get_biggest_protein(False, aminoseq)
    if expected_output is not None:
        assert result is not None
        protein, start, end = result
        assert protein == expected_output
        assert aminoseq[start:end + 1] == protein
    else:
        assert result is None
