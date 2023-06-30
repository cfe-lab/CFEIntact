
import pytest
import os

import intact.intact as intact
from intact.intact import is_nonhiv, IntactnessError

class BlastRow:
    def __init__(self, qstart, qend, qlen):
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen


def test_is_nonhiv_empty_blast_rows():
    blast_rows = []
    result = is_nonhiv("id", 5, blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.NONHIV_ERROR


def test_is_nonhiv_low_ratio():
    blast_rows = [
        BlastRow(0, 100, 1000),
        BlastRow(200, 300, 1000),
        BlastRow(400, 500, 1000)
    ]
    result = is_nonhiv("id", 1000, blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.NONHIV_ERROR


def test_is_nonhiv_high_ratio():
    blast_rows = [
        BlastRow(0, 200, 1000),
        BlastRow(200, 400, 1000),
        BlastRow(400, 900, 1000)
    ]
    assert is_nonhiv("id", 1000, blast_rows) is None


def test_is_nonhiv_high_ratio_reversed():
    blast_rows = [
        BlastRow(0, 200, 1000),
        BlastRow(200, 400, 1000),
        BlastRow(900, 400, 1000)
    ]
    assert is_nonhiv("id", 1000, blast_rows) is None


def test_is_nonhiv_single_blast_row():
    blast_rows = [BlastRow(0, 200, 1000)]
    result = is_nonhiv("id", 1000, blast_rows)
    assert isinstance(result, IntactnessError)
    assert result.error == intact.NONHIV_ERROR
