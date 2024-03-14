
from typing import List

import cfeintact.defect as defect
from cfeintact.intact import is_nonhiv, HolisticInfo
from cfeintact.defect import Defect


class BlastRow:
    def __init__(self, qstart, qend, qlen):
        self.qstart = qstart
        self.qend = qend
        self.qlen = qlen


def test_is_nonhiv_empty_blast_rows():
    blast_rows: List[BlastRow] = []
    result = is_nonhiv(HolisticInfo(), "id", 5, blast_rows)
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.NonHIV)


def test_is_nonhiv_low_ratio():
    blast_rows = [
        BlastRow(0, 100, 1000),
        BlastRow(200, 300, 1000),
        BlastRow(400, 500, 1000)
    ]
    result = is_nonhiv(HolisticInfo(), "id", 1000, blast_rows)
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.NonHIV)


def test_is_nonhiv_high_ratio():
    blast_rows = [
        BlastRow(0, 200, 1000),
        BlastRow(200, 400, 1000),
        BlastRow(400, 900, 1000)
    ]
    assert is_nonhiv(HolisticInfo(), "id", 1000, blast_rows) is None


def test_is_nonhiv_high_ratio_reversed():
    blast_rows = [
        BlastRow(0, 200, 1000),
        BlastRow(200, 400, 1000),
        BlastRow(900, 400, 1000)
    ]
    assert is_nonhiv(HolisticInfo(), "id", 1000, blast_rows) is None


def test_is_nonhiv_single_blast_row():
    blast_rows = [BlastRow(0, 200, 1000)]
    result = is_nonhiv(HolisticInfo(), "id", 1000, blast_rows)
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.NonHIV)
