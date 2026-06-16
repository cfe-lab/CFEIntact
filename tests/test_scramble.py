
import pytest
from typing import List, Dict, Tuple

import cfeintact.intact as intact
import cfeintact.defect as defect
from cfeintact.intact import (
    is_scrambled,
    contains_internal_inversion,
    most_frequent_element,
    SubjectPlacement,
    normalized_subject_interval,
    subject_segment,
    subject_segment_length,
    find_exact_occurrences,
    candidate_subject_placements,
    is_subject_location_informative,
    has_monotone_subject_assignment,
    has_monotone_subject_assignment_for_group,
    row_has_enough_order_evidence,
    MIN_ORDER_EVIDENCE_LENGTH,
    MIN_INTERNAL_INVERSION_EVIDENCE_LENGTH,
)
from cfeintact.defect import Defect
from cfeintact.blastrow import BlastRow


# ---------------------------------------------------------------------------
# Helper: build a subject string where each row gets a unique segment at its
# interval.  Rows MUST have non-overlapping subject intervals.
# ---------------------------------------------------------------------------

def subject_from_rows(blast_rows, length=10000):
    seq = bytearray(b"N" * length)
    for i, row in enumerate(blast_rows):
        lo = min(row.sstart, row.send) - 1
        hi = max(row.sstart, row.send)
        tag = f"SEGM{i:04d}X"
        reps = (hi - lo) // len(tag) + 1
        val = (tag * reps)[: hi - lo]
        seq[lo:hi] = val.encode()
    return seq.decode()


# ---------------------------------------------------------------------------
# Existing parametrized tests (unchanged)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lst, expected", [
    ([1, 2, 3, 4, 2, 2, 3, 1, 4, 4, 4], 4),
    ([1, 2, 3, 4, 2, 2, 3, 1], 2),
    ([], None),
    ([5], 5),
    ([1, 2, 3, 4, 5], 1),
    ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 1),
    ([1, 2, 3, 2, 4, 5, 4, 6, 6, 6], 6)
])
def test_most_frequent_element(lst, expected):
    assert most_frequent_element(lst) == expected


@pytest.mark.parametrize("lst, expected", [
    ([1, 2, 3, 4, 5], True),
    ([1, 3, 2, 4, 5], False),
    ([], True),
    ([5], True),
    ([1, 1, 2, 2, 3, 3], True),
])
def test_is_sorted(lst, expected):
    assert intact.is_sorted(lst) == expected


# ---------------------------------------------------------------------------
# Mock row and factory
# ---------------------------------------------------------------------------

class MockBlastRow:
    def __init__(self, sstart, send, qstart, sstrand, qend=None, sseqid="mock_subject"):
        self.sstart = sstart
        self.send = send
        self.qstart = qstart
        self.sstrand = sstrand
        self.qend = qend if qend is not None else qstart
        self.qlen = 9000
        self.slen = 9000
        self.sseqid = sseqid


def mk_blast_row(sstart, send, qstart, sstrand, qend=None, sseqid="mock_subject"):
    return MockBlastRow(sstart, send, qstart, sstrand, qend=qend, sseqid=sseqid)


# ---------------------------------------------------------------------------
# Unit tests for helpers
# ---------------------------------------------------------------------------

class TestNormalizedSubjectInterval:
    def test_plus_strand(self):
        sp = normalized_subject_interval(mk_blast_row(900, 950, 10, "plus"))
        assert sp == SubjectPlacement(900, 950)

    def test_minus_strand(self):
        sp = normalized_subject_interval(mk_blast_row(950, 900, 10, "minus"))
        assert sp == SubjectPlacement(900, 950)

    def test_zero_length(self):
        sp = normalized_subject_interval(mk_blast_row(500, 500, 10, "plus"))
        assert sp == SubjectPlacement(500, 500)


class TestSubjectSegment:
    def test_extracts_correct_slice(self):
        row = mk_blast_row(5, 9, 10, "plus")
        subj = "ABCDEFGHIJKLMNOP"
        assert subject_segment(row, subj) == "EFGHI"

    def test_minus_strand(self):
        row = mk_blast_row(9, 5, 10, "minus")
        subj = "ABCDEFGHIJKLMNOP"
        assert subject_segment(row, subj) == "EFGHI"

    def test_raises_on_overflow(self):
        row = mk_blast_row(95, 105, 10, "plus")
        subj = "A" * 100
        with pytest.raises(ValueError, match="exceeds"):
            subject_segment(row, subj)


class TestSubjectSegmentLength:
    def test_normal(self):
        assert subject_segment_length(mk_blast_row(100, 200, 10, "plus")) == 101

    def test_reversed(self):
        assert subject_segment_length(mk_blast_row(200, 100, 10, "plus")) == 101

    def test_single_base(self):
        assert subject_segment_length(mk_blast_row(50, 50, 10, "plus")) == 1


class TestFindExactOccurrences:
    def test_single_occurrence(self):
        result = find_exact_occurrences("ABC", "XXABCYY")
        assert result == [SubjectPlacement(3, 5)]

    def test_multiple_non_overlapping(self):
        result = find_exact_occurrences("ABA", "ABAXXABA")
        assert result == [SubjectPlacement(1, 3), SubjectPlacement(6, 8)]

    def test_overlapping(self):
        result = find_exact_occurrences("AAA", "AAAAA")
        assert result == [
            SubjectPlacement(1, 3),
            SubjectPlacement(2, 4),
            SubjectPlacement(3, 5),
        ]

    def test_empty_needle(self):
        assert find_exact_occurrences("", "ABC") == []

    def test_no_match(self):
        assert find_exact_occurrences("XYZ", "ABC") == []

    def test_whole_string(self):
        result = find_exact_occurrences("HELLO", "HELLO")
        assert result == [SubjectPlacement(1, 5)]


class TestCandidateSubjectPlacements:
    def test_includes_original(self):
        subj = "NNNABCDEABCDFNNN"
        row = mk_blast_row(4, 8, 10, "plus")
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        placements = candidate_subject_placements(row, subj, cache, row.sseqid)
        assert SubjectPlacement(4, 8) in placements

    def test_multiple_copies(self):
        subj = "NNNABCDEXXXABCDENNN"
        row = mk_blast_row(4, 8, 10, "plus")
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        placements = candidate_subject_placements(row, subj, cache, row.sseqid)
        assert len(placements) == 2
        assert SubjectPlacement(4, 8) in placements
        assert SubjectPlacement(12, 16) in placements

    def test_caches_results(self):
        subj = "ABCABC"
        row_a = mk_blast_row(1, 3, 10, "plus")
        row_b = mk_blast_row(4, 6, 20, "plus")
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        p1 = candidate_subject_placements(row_a, subj, cache, row_a.sseqid)
        p2 = candidate_subject_placements(row_b, subj, cache, row_b.sseqid)
        assert p1 is p2


class TestIsSubjectLocationInformative:
    def test_short_segment_not_informative(self):
        subj = "A" * 10000
        row = mk_blast_row(500, 500 + MIN_ORDER_EVIDENCE_LENGTH - 1, 10, "plus")
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        assert not is_subject_location_informative(row, subj, cache, row.sseqid)

    def test_long_enough_unique_is_informative(self):
        subj = "N" * 200 + "X" * MIN_ORDER_EVIDENCE_LENGTH + "N" * 9700
        row = mk_blast_row(201, 200 + MIN_ORDER_EVIDENCE_LENGTH, 10, "plus")
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        assert is_subject_location_informative(row, subj, cache, row.sseqid)

    def test_long_enough_duplicate_not_informative(self):
        seg = "X" * MIN_ORDER_EVIDENCE_LENGTH
        subj = "N" * 100 + seg + "N" * 500 + seg + "N" * 9300
        row = mk_blast_row(101, 100 + MIN_ORDER_EVIDENCE_LENGTH, 10, "plus")
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        assert not is_subject_location_informative(row, subj, cache, row.sseqid)


class TestRowHasEnoughOrderEvidence:
    def test_short_row(self):
        row = mk_blast_row(1, MIN_ORDER_EVIDENCE_LENGTH - 1, 10, "plus")
        assert not row_has_enough_order_evidence(row, "N" * 10000)

    def test_long_enough(self):
        row = mk_blast_row(1, MIN_ORDER_EVIDENCE_LENGTH, 10, "plus")
        assert row_has_enough_order_evidence(row, "N" * 10000)


class TestHasMonotoneSubjectAssignmentForGroup:
    def test_empty_rows(self):
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        assert has_monotone_subject_assignment_for_group([], "A" * 100, "plus", cache) is True

    def test_plus_sorted(self):
        rows = [
            mk_blast_row(100, 200, 10, "plus"),
            mk_blast_row(300, 400, 20, "plus"),
        ]
        subj = subject_from_rows(rows)
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        assert has_monotone_subject_assignment_for_group(rows, subj, "plus", cache) is True

    def test_plus_unsorted_unique(self):
        rows = [
            mk_blast_row(100, 200, 10, "plus"),
            mk_blast_row(50, 150, 20, "plus"),
        ]
        subj = subject_from_rows(rows)
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        assert has_monotone_subject_assignment_for_group(rows, subj, "plus", cache) is False

    def test_minus_uses_end(self):
        """For minus: use end (high coord) for monotonicity check."""
        rows = [
            mk_blast_row(1000, 900, 10, "minus"),   # end=1000
            mk_blast_row(800, 700, 20, "minus"),    # end=800
        ]
        subj = subject_from_rows(rows)
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        # end: 1000 -> 800, nonincreasing -> monotone
        assert has_monotone_subject_assignment_for_group(rows, subj, "minus", cache) is True

    def test_minus_unsorted_unique(self):
        rows = [
            mk_blast_row(800, 700, 10, "minus"),    # end=800
            mk_blast_row(1000, 900, 20, "minus"),   # end=1000
        ]
        subj = subject_from_rows(rows)
        cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}
        # end: 800 -> 1000, 1000 > 800 -> not nonincreasing -> not monotone
        assert has_monotone_subject_assignment_for_group(rows, subj, "minus", cache) is False


# ---------------------------------------------------------------------------
# Behavioral tests — edge cases
# ---------------------------------------------------------------------------

def test_is_scrambled_no_alignment():
    blast_rows: List[BlastRow] = []
    assert is_scrambled("id", blast_rows, {"mock_subject": ""}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": ""}) is None


# ---------------------------------------------------------------------------
# Monotone unique rows (not scrambled)
# ---------------------------------------------------------------------------

def test_plus_strand_sorted():
    blast_rows = [
        mk_blast_row(sstart=700, send=810, qstart=50, sstrand="plus"),
        mk_blast_row(sstart=900, send=1010, qstart=200, sstrand="plus"),
    ]
    subj = subject_from_rows(blast_rows)
    assert is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subj}) is None


def test_minus_strand_sorted():
    blast_rows = [
        mk_blast_row(sstart=600, send=500, qstart=200, sstrand="minus"),
        mk_blast_row(sstart=400, send=300, qstart=50, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    assert is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subj}) is None


def test_single_row_plus_strand():
    blast_rows = [mk_blast_row(sstart=700, send=810, qstart=5, sstrand="plus")]
    subj = subject_from_rows(blast_rows)
    assert is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subj}) is None


def test_single_row_minus_strand():
    blast_rows = [mk_blast_row(sstart=900, send=910, qstart=5, sstrand="minus")]
    subj = subject_from_rows(blast_rows)
    assert is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subj}) is None


# ---------------------------------------------------------------------------
# Unique unsorted regions -> Scramble
# ---------------------------------------------------------------------------

def test_plus_strand_unsorted_unique():
    blast_rows = [
        mk_blast_row(sstart=880, send=990, qstart=700, sstrand="plus"),
        mk_blast_row(sstart=100, send=210, qstart=1100, sstrand="plus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


def test_minus_strand_unsorted_unique():
    """For minus: monotone if end (high coord) nonincreasing.
    Row A (q=700, end=980) -> Row B (q=800, end=990): 980 > 990? No -> scramble."""
    blast_rows = [
        mk_blast_row(sstart=990, send=980, qstart=700, sstrand="minus"),
        mk_blast_row(sstart=1000, send=990, qstart=800, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Test A: LTR-style terminal cross-hit should not be scramble
# ---------------------------------------------------------------------------

def test_ltr_terminal_cross_hit_not_scrambled():
    """LTR-like terminal repeats: subject has identical prefix and suffix.
    BLAST maps query 5' to subject 3' and vice versa. Exact-match solver
    resolves the ambiguity -> no Scramble, no InternalInversion."""
    segment_len = 200
    prefix = "X" * segment_len
    suffix = "X" * segment_len  # identical
    middle = "Y" * 1000
    subject = prefix + middle + suffix
    query_len = len(prefix) + len(middle) + segment_len  # matching length

    blast_rows = [
        # Full-length alignment
        mk_blast_row(sstart=1, send=query_len, qstart=1, qend=query_len, sstrand="plus"),
        # Query 5' -> subject 3' (LTR cross-hit)
        mk_blast_row(sstart=len(prefix) + len(middle) + 1, send=query_len,
                      qstart=1, qend=segment_len, sstrand="plus"),
        # Query 3' -> subject 5' (LTR cross-hit)
        mk_blast_row(sstart=1, send=segment_len,
                      qstart=len(prefix) + len(middle) + 1, qend=query_len, sstrand="plus"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert result is None
    assert contains_internal_inversion("id", blast_rows, {"mock_subject": subject}) is None


# ---------------------------------------------------------------------------
# Test B: Real reorder of unique regions should still be scramble
# ---------------------------------------------------------------------------

def test_real_reorder_unique_regions():
    """Three unique regions A, B, C; query order maps to subject order A, C, B."""
    subject = "N" * 2000
    # Manually build with distinct segments
    seg_a = "AAAA" * 50   # 200 bp
    seg_b = "BBBB" * 50   # 200 bp
    seg_c = "CCCC" * 50   # 200 bp
    subject = seg_a + seg_b + seg_c + "N" * 9400

    blast_rows = [
        mk_blast_row(sstart=1, send=200, qstart=100, qend=299, sstrand="plus"),     # A at 100
        mk_blast_row(sstart=401, send=600, qstart=500, qend=699, sstrand="plus"),   # C at 500 (but B is at 201-400)
        mk_blast_row(sstart=201, send=400, qstart=300, qend=499, sstrand="plus"),   # B at 300
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Test C: Reorder explainable by exact duplicated region should not be scramble
# ---------------------------------------------------------------------------

def test_exact_duplicate_explains_reorder():
    """Earlier query row maps to a late copy; the segment also occurs earlier."""
    segment = "EXACT_DUPLICATE_" * 10
    subject = (
        "N" * 99 +
        segment +          # positions 100-259
        "N" * 740 +
        segment +          # positions 1000-1159 (same text)
        "N" * (10000 - 1160)
    )
    blast_rows = [
        mk_blast_row(sstart=1000, send=1159, qstart=100, qend=259, sstrand="plus"),
        mk_blast_row(sstart=400, send=499, qstart=300, qend=399, sstrand="plus"),
    ]
    assert is_scrambled("id", blast_rows, {"mock_subject": subject}) is None


# ---------------------------------------------------------------------------
# Duplicate not automatically forgiven
# ---------------------------------------------------------------------------

def test_scramble_duplicate_not_automatically_forgiven():
    """Unique rows constrain the assignment; the duplicate cannot fix order."""
    dup = "DUPBLOCK" * 20
    subject = (
        "N" * 49 +
        dup +              # positions 50-209 (early copy)
        "N" * 90 +
        "UNIQUE_A" * 20 +  # positions 300-459
        "N" * 40 +
        "UNIQUE_D" * 20 +  # positions 500-659
        "N" * 40 +
        dup +              # positions 700-859 (late copy)
        "N" * (10000 - 860)
    )
    blast_rows = [
        mk_blast_row(sstart=300, send=459, qstart=100, qend=259, sstrand="plus"),
        mk_blast_row(sstart=700, send=859, qstart=200, qend=359, sstrand="plus"),
        mk_blast_row(sstart=500, send=659, qstart=300, qend=459, sstrand="plus"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Test D: Short HSPs should not create order evidence
# ---------------------------------------------------------------------------

def test_short_hsp_does_not_create_scramble():
    """A row below MIN_ORDER_EVIDENCE_LENGTH should not cause a scramble."""
    short_len = MIN_ORDER_EVIDENCE_LENGTH - 1
    blast_rows = [
        mk_blast_row(sstart=500, send=500 + short_len - 1, qstart=100, qend=100 + short_len - 1, sstrand="plus"),
        mk_blast_row(sstart=100, send=100 + short_len - 1, qstart=200, qend=200 + short_len - 1, sstrand="plus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj})
    assert result is None  # no evidence -> no decision -> not scrambled


def test_short_hsp_does_not_hide_real_scramble():
    """Short rows are ignored; long unique rows still detect scramble."""
    blast_rows = [
        # Short row (< 100 bp) - should be ignored
        mk_blast_row(sstart=50, send=100, qstart=50, qend=100, sstrand="plus"),
        # Long rows (> 100 bp) - provide order evidence
        mk_blast_row(sstart=900, send=1010, qstart=500, qend=610, sstrand="plus"),
        mk_blast_row(sstart=200, send=310, qstart=700, qend=810, sstrand="plus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Test E: Cache must be subject-specific
# ---------------------------------------------------------------------------

def test_cache_is_subject_specific():
    """Same segment string at different coordinates in different subjects
    must produce separate cache entries."""
    subj_a = "ABCABC"
    subj_b = "XXXXABCXXXX"
    row = mk_blast_row(1, 3, 10, "plus")  # segment = "ABC"

    cache: Dict[Tuple[str, str], List[SubjectPlacement]] = {}

    # Look up in subject A (sseqid="subj_A")
    placements_a = candidate_subject_placements(row, subj_a, cache, "subj_A")
    assert placements_a == [SubjectPlacement(1, 3), SubjectPlacement(4, 6)]

    # Look up in subject B (sseqid="subj_B") - should not reuse cache
    placements_b = candidate_subject_placements(row, subj_b, cache, "subj_B")
    assert placements_b == [SubjectPlacement(5, 7)]
    assert placements_a is not placements_b  # different cached objects


# ---------------------------------------------------------------------------
# Test F: Minus-strand monotone assignment must use placement.end
# ---------------------------------------------------------------------------

def test_minus_uses_placement_end():
    """For minus, the greedy assignment uses end (high coord)."""
    rows = [
        mk_blast_row(sstart=1000, send=900, qstart=10, sstrand="minus"),   # end=1000
        mk_blast_row(sstart=800, send=700, qstart=20, sstrand="minus"),    # end=800
    ]
    subj = subject_from_rows(rows)
    result = is_scrambled("id", rows, {"mock_subject": subj})
    assert result is None  # end: 1000 -> 800, nonincreasing -> not scrambled


def test_minus_unsorted_uses_end():
    """Minus with end increasing should be scramble."""
    rows = [
        mk_blast_row(sstart=700, send=800, qstart=10, sstrand="minus"),   # end=800
        mk_blast_row(sstart=900, send=1000, qstart=20, sstrand="minus"),   # end=1000
    ]
    subj = subject_from_rows(rows)
    result = is_scrambled("id", rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Test G: Internal inversion ignores repeated regions
# ---------------------------------------------------------------------------

def test_inversion_ignores_duplicate():
    """Plus and minus rows where the minus row maps to a duplicated segment."""
    seg = "DUP" * 50
    subject = (
        "N" * 99 +
        seg +              # positions 100-248
        "N" * 751 +
        seg +              # positions 1000-1148 (duplicate)
        "N" * (10000 - 1149)
    )
    blast_rows = [
        mk_blast_row(sstart=100, send=248, qstart=100, qend=248, sstrand="plus"),
        mk_blast_row(sstart=1148, send=1000, qstart=500, qend=648, sstrand="minus"),
    ]
    assert contains_internal_inversion("id", blast_rows, {"mock_subject": subject}) is None


# ---------------------------------------------------------------------------
# Test H: Internal inversion requires sufficient minority-strand evidence
# ---------------------------------------------------------------------------

def test_inversion_insufficient_minority_evidence():
    """Tiny unique minority-strand row below threshold -> no inversion."""
    short_len = MIN_INTERNAL_INVERSION_EVIDENCE_LENGTH - 1
    subject = "N" * 200 + "U" * 200 + "N" * 9600
    blast_rows = [
        mk_blast_row(sstart=201, send=400, qstart=100, qend=299, sstrand="plus"),  # 200bp, unique
        mk_blast_row(sstart=201, send=200 + short_len, qstart=500, qend=500 + short_len - 1, sstrand="minus"),  # too short
    ]
    assert contains_internal_inversion("id", blast_rows, {"mock_subject": subject}) is None


def test_inversion_sufficient_minority_evidence():
    """Long enough unique minority-strand row -> inversion."""
    subject = "N" * 200 + "U" * 200 + "N" * 200 + "V" * 200 + "N" * 9200
    blast_rows = [
        mk_blast_row(sstart=201, send=400, qstart=100, qend=299, sstrand="plus"),  # 200bp, unique
        mk_blast_row(sstart=601, send=800, qstart=500, qend=699, sstrand="minus"),  # 200bp, unique
    ]
    result = contains_internal_inversion("id", blast_rows, {"mock_subject": subject})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.InternalInversion)


# ---------------------------------------------------------------------------
# Internal inversion: full pipeline (is_scrambled + contains_internal_inversion)
# ---------------------------------------------------------------------------

def test_internal_inversion_detected():
    """Mixed-strand evidence in informative rows -> Inversion."""
    blast_rows = [
        mk_blast_row(sstart=900, send=1010, qstart=700, qend=810, sstrand="plus"),
        mk_blast_row(sstart=1200, send=1310, qstart=800, qend=900, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.InternalInversion)


def test_internal_inversion_not_reported_for_duplicate_only():
    """Mixed strand only in duplicate regions -> no inversion."""
    seg = "DUPSTRAND" * 50
    subject = (
        "N" * 99 +
        seg +              # positions 100-498
        "N" * 501 +
        seg +              # positions 1000-1398
        "N" * (10000 - 1399)
    )
    blast_rows = [
        mk_blast_row(sstart=100, send=498, qstart=100, qend=498, sstrand="plus"),
        mk_blast_row(sstart=1398, send=1000, qstart=500, qend=898, sstrand="minus"),
    ]
    assert is_scrambled("id", blast_rows, {"mock_subject": subject}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subject}) is None


def test_internal_inversion_not_hidden_by_duplicate():
    """Informative rows with mixed strands -> inversion, even with duplicates."""
    dup = "X" * 180
    unique_plus = "P" * 160
    unique_minus = "M" * 160
    subject = (
        "N" * 99 +
        dup +              # positions 100-279 (copy A)
        "N" * 220 +
        unique_plus +      # positions 500-659
        "N" * 340 +
        unique_minus +     # positions 1000-1159
        "N" * 840 +
        dup +              # positions 2000-2179 (copy B)
        "N" * 7821
    )
    blast_rows = [
        mk_blast_row(sstart=100, send=279, qstart=100, qend=279, sstrand="plus"),
        mk_blast_row(sstart=500, send=659, qstart=300, qend=459, sstrand="plus"),
        mk_blast_row(sstart=1159, send=1000, qstart=500, qend=659, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subject})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.InternalInversion)


# ---------------------------------------------------------------------------
# Mixed direction tests
# ---------------------------------------------------------------------------

def test_mixed_direction_returns_scramble():
    blast_rows = [
        mk_blast_row(sstart=900, send=1010, qstart=700, qend=810, sstrand="plus"),
        mk_blast_row(sstart=1200, send=1310, qstart=800, qend=900, sstrand="minus"),
        mk_blast_row(sstart=1400, send=1510, qstart=900, qend=1040, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


def test_multiple_directions_and_inversions():
    blast_rows = [
        mk_blast_row(sstart=900, send=1010, qstart=700, qend=810, sstrand="plus"),
        mk_blast_row(sstart=1200, send=1310, qstart=800, qend=900, sstrand="minus"),
        mk_blast_row(sstart=1400, send=1510, qstart=900, qend=1040, sstrand="minus"),
        mk_blast_row(sstart=1600, send=1710, qstart=1000, qend=1140, sstrand="plus"),
        mk_blast_row(sstart=1800, send=1910, qstart=1000, qend=1140, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Unknown sseqid
# ---------------------------------------------------------------------------

def test_unknown_sseqid_skipped():
    """If no sseqid is in the subject dict, row is skipped -> no scramble."""
    blast_rows = [
        mk_blast_row(sstart=900, send=1010, qstart=50, qend=150, sstrand="plus",
                     sseqid="unknown_ref"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": "A" * 10000})
    assert result is None  # skipped, not scrambled


def test_unknown_sseqid_does_not_affect_known():
    """Rows with known sseqid still evaluated; unknown rows are skipped."""
    blast_rows = [
        mk_blast_row(sstart=900, send=1010, qstart=50, qend=150, sstrand="plus",
                     sseqid="known_ref"),
        mk_blast_row(sstart=100, send=210, qstart=200, qend=310, sstrand="plus",
                     sseqid="unknown_ref"),
    ]
    subj = subject_from_rows([blast_rows[0]])  # only for known_ref
    result = is_scrambled("id", blast_rows, {"known_ref": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# has_monotone_subject_assignment direct tests
# ---------------------------------------------------------------------------

def test_has_monotone_subject_assignment_no_evidence_rows():
    """No rows with enough evidence -> returns False."""
    assert has_monotone_subject_assignment([], {"s": "A" * 10000}, "plus") is False

    # Rows exist but below threshold
    short = MIN_ORDER_EVIDENCE_LENGTH - 1
    rows = [mk_blast_row(1, short, 10, "plus")]
    assert has_monotone_subject_assignment(rows, {"mock_subject": "N" * 10000}, "plus") is False


def test_has_monotone_subject_assignment_single():
    rows = [mk_blast_row(100, 200, 10, "plus")]
    subj = subject_from_rows(rows)
    assert has_monotone_subject_assignment(rows, {"mock_subject": subj}, "plus") is True


def test_has_monotone_subject_assignment_plus_unsorted_unique():
    rows = [
        mk_blast_row(100, 200, 10, "plus"),
        mk_blast_row(50, 150, 20, "plus"),
    ]
    subj = subject_from_rows(rows)
    assert has_monotone_subject_assignment(rows, {"mock_subject": subj}, "plus") is False


def test_has_monotone_subject_assignment_minus_unsorted_unique():
    """For minus: end nonincreasing. Row A end=200, Row B end=150 -> nonincreasing? No, 200->150 is ok.
    Row A end=150, Row B end=200 -> 150->200 not nonincreasing -> fails."""
    rows = [
        mk_blast_row(200, 100, 10, "minus"),   # end=200
        mk_blast_row(250, 150, 20, "minus"),   # end=250, 250 > 200 -> not nonincreasing
    ]
    subj = subject_from_rows(rows)
    assert has_monotone_subject_assignment(rows, {"mock_subject": subj}, "minus") is False


def test_has_monotone_subject_assignment_minus_sorted_unique():
    rows = [
        mk_blast_row(250, 150, 10, "minus"),   # end=250
        mk_blast_row(200, 100, 20, "minus"),   # end=200, 200 <= 250 -> nonincreasing
    ]
    subj = subject_from_rows(rows)
    assert has_monotone_subject_assignment(rows, {"mock_subject": subj}, "minus") is True


# ---------------------------------------------------------------------------
# LTR analogue (no hard-coded 622)
# ---------------------------------------------------------------------------

def test_scramble_ltr_analogue():
    """Identical prefix/suffix blocks (like HIV LTRs) do not trigger scramble."""
    prefix = "IDENTICAL_PREFIX_" * 11
    suffix = "IDENTICAL_PREFIX_" * 11
    middle = "MIDDLE_UNIQUE_" * 100
    subject = prefix + middle + suffix
    blast_rows = [
        mk_blast_row(
            sstart=len(prefix) + len(middle) + 1,
            send=len(prefix) + len(middle) + len(suffix),
            qstart=50,
            qend=50 + len(prefix) - 1,
            sstrand="plus",
        ),
        mk_blast_row(
            sstart=len(prefix) + 1,
            send=len(prefix) + len(middle),
            qstart=800,
            qend=800 + len(middle) - 1,
            sstrand="plus",
        ),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert result is None


def test_scramble_minus_duplicate_explains_reorder():
    """Minus-strand rows whose apparent disorder is explained by duplicates."""
    seg = "REPEATED_MINUS_" * 10
    subject = (
        "N" * 499 +
        seg +              # positions 500-659
        "N" * 340 +
        seg +              # positions 1000-1159
        "N" * (10000 - 1160)
    )
    blast_rows = [
        mk_blast_row(sstart=1159, send=1000, qstart=100, qend=259, sstrand="minus"),
        mk_blast_row(sstart=299, send=200, qstart=300, qend=399, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert result is None
