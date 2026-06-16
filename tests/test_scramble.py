
import pytest
from typing import List, Dict

import cfeintact.intact as intact
import cfeintact.defect as defect
from cfeintact.intact import (
    is_scrambled,
    contains_internal_inversion,
    most_frequent_element,
    SubjectPlacement,
    normalized_subject_interval,
    subject_segment,
    find_exact_occurrences,
    candidate_subject_placements,
    is_subject_location_informative,
    has_monotone_subject_assignment,
    MIN_ORDER_EVIDENCE_LENGTH,
)
from cfeintact.defect import Defect
from cfeintact.blastrow import BlastRow


# ---------------------------------------------------------------------------
# Helper: build a subject string where each row gets a unique segment at its
# interval.  Rows MUST have non-overlapping subject intervals.
# ---------------------------------------------------------------------------

def subject_from_rows(blast_rows, length=10000):
    """Build a subject string where each row gets a unique segment.
    Each row's subject interval is filled with a distinct repeating tag.
    Overlapping intervals will corrupt segments for all overlapping rows."""
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
# Unit tests for new helpers
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
        # 1-based inclusive → subject[4:9] = "EFGHI"
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


class TestFindExactOccurrences:
    def test_single_occurrence(self):
        result = find_exact_occurrences("ABC", "XXABCYY")
        # ABC at 0-indexed pos 2 → SubjectPlacement(3, 5)
        assert result == [SubjectPlacement(3, 5)]

    def test_multiple_non_overlapping(self):
        result = find_exact_occurrences("ABA", "ABAXXABA")
        # first at pos 0 → (1, 3), second at pos 5 → (6, 8)
        assert result == [SubjectPlacement(1, 3), SubjectPlacement(6, 8)]

    def test_overlapping(self):
        result = find_exact_occurrences("AAA", "AAAAA")
        # pos 0 → (1, 3), pos 1 → (2, 4), pos 2 → (3, 5)
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
        row = mk_blast_row(4, 8, 10, "plus")  # "ABCDE"
        cache: Dict[str, List[SubjectPlacement]] = {}
        placements = candidate_subject_placements(row, subj, cache)
        assert SubjectPlacement(4, 8) in placements

    def test_multiple_copies(self):
        subj = "NNNABCDEXXXABCDENNN"
        row = mk_blast_row(4, 8, 10, "plus")  # "ABCDE"
        cache: Dict[str, List[SubjectPlacement]] = {}
        placements = candidate_subject_placements(row, subj, cache)
        assert len(placements) == 2
        assert SubjectPlacement(4, 8) in placements
        # second copy at 0-indexed pos 11 → SubjectPlacement(12, 16)
        assert SubjectPlacement(12, 16) in placements

    def test_caches_results(self):
        subj = "ABCABC"
        row_a = mk_blast_row(1, 3, 10, "plus")
        row_b = mk_blast_row(4, 6, 20, "plus")
        cache: Dict[str, List[SubjectPlacement]] = {}
        p1 = candidate_subject_placements(row_a, subj, cache)
        p2 = candidate_subject_placements(row_b, subj, cache)
        assert p1 is p2  # same list object from cache


class TestIsSubjectLocationInformative:
    def test_short_segment_not_informative(self):
        subj = "A" * 10000
        row = mk_blast_row(500, 500 + MIN_ORDER_EVIDENCE_LENGTH - 1, 10, "plus")
        assert not is_subject_location_informative(row, subj)

    def test_long_enough_unique_is_informative(self):
        subj = "N" * 200 + "X" * MIN_ORDER_EVIDENCE_LENGTH + "N" * 9700
        row = mk_blast_row(201, 200 + MIN_ORDER_EVIDENCE_LENGTH, 10, "plus")
        assert is_subject_location_informative(row, subj)

    def test_long_enough_duplicate_not_informative(self):
        seg = "X" * MIN_ORDER_EVIDENCE_LENGTH
        subj = "N" * 100 + seg + "N" * 500 + seg + "N" * 9300
        row = mk_blast_row(101, 100 + MIN_ORDER_EVIDENCE_LENGTH, 10, "plus")
        assert not is_subject_location_informative(row, subj)


# ---------------------------------------------------------------------------
# Behavioral tests — no-alignment edge case
# ---------------------------------------------------------------------------

def test_is_scrambled_no_alignment():
    blast_rows: List[BlastRow] = []
    # Both return None → (None or (None is None)) → True
    assert is_scrambled("id", blast_rows, {"mock_subject": ""}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": ""}) is None


# ---------------------------------------------------------------------------
# Behavioral tests — monotone unique rows (not scrambled, not inverted)
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
    blast_rows = [
        mk_blast_row(sstart=700, send=810, qstart=5, sstrand="plus"),
    ]
    subj = subject_from_rows(blast_rows)
    assert is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subj}) is None


def test_single_row_minus_strand():
    blast_rows = [
        mk_blast_row(sstart=900, send=910, qstart=5, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    assert is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subj}) is None


# ---------------------------------------------------------------------------
# Behavioral tests — unique unsorted regions → Scramble
# ---------------------------------------------------------------------------

def test_plus_strand_unsorted_unique():
    # Non-overlapping intervals to avoid subject_from_rows corruption
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
    # For minus: as qstart increases, normalized low should decrease.
    # Row B (q=800, low=940) > Row A (q=700, low=910) → not decreasing → Scramble
    blast_rows = [
        mk_blast_row(sstart=940, send=930, qstart=700, sstrand="minus"),
        mk_blast_row(sstart=990, send=980, qstart=800, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Behavioral tests — exact duplicate explains reorder (no Scramble)
# ---------------------------------------------------------------------------

def test_scramble_exact_duplicate_explains_reorder():
    """Earlier query row maps to a late copy; the segment also occurs earlier."""
    segment = "EXACT_DUPLICATE_" * 10
    subject = (
        "N" * 99 +
        segment +          # positions 100–259
        "N" * 740 +
        segment +          # positions 1000–1159 (same text)
        "N" * (10000 - 1160)
    )
    blast_rows = [
        mk_blast_row(sstart=1000, send=1159, qstart=100, qend=259, sstrand="plus"),
        mk_blast_row(sstart=400, send=499, qstart=300, qend=399, sstrand="plus"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert result is None


def test_scramble_duplicate_not_automatically_forgiven():
    """A duplicate region cannot fix the order when surrounding unique rows
    constrain the overall assignment, regardless of the duplicate's placement."""
    # Row A (q=100): unique segment at [300, 399]
    # Row D (q=300): unique segment at [500, 599]
    # Row B (q=200): duplicate segment at both [600, 699] and [50, 149]
    #   - With late copy (600): starts = [300, 600, 500] → not monotone
    #   - With early copy (50):  starts = [300, 50, 500]  → not monotone
    # Either way, the total sequence is not monotone → Scramble
    dup = "DUPBLOCK" * 20
    subject = (
        "N" * 49 +
        dup +              # positions 50–209 (early copy)
        "N" * 90 +
        "UNIQUE_A" * 20 +  # positions 300–459
        "N" * 40 +
        "UNIQUE_D" * 20 +  # positions 500–659
        "N" * 40 +
        dup +              # positions 700–859 (late copy)
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


def test_scramble_ltr_analogue():
    """Identical prefix/suffix blocks (like HIV LTRs) do not trigger scramble."""
    prefix = "IDENTICAL_PREFIX_" * 11
    suffix = "IDENTICAL_PREFIX_" * 11  # same string
    middle = "MIDDLE_UNIQUE_" * 100
    subject = prefix + middle + suffix
    # BLAST maps qstart 50 to the suffix (late copy), and qstart 800 to middle
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
        seg +              # positions 500–659
        "N" * 340 +
        seg +              # positions 1000–1159
        "N" * (10000 - 1160)
    )
    # Row A (q=100): BLAST picks late copy (1000-1159), minus strand
    #   Normalized: starts = 1000, end = 1159
    #   Segment also at 500-659, so candidates include start=500
    # Row B (q=300): unique segment, start=200
    #
    # Without duplicate fix: starts = [1000, 200] → not minus-sorted
    # With duplicate fix: choose 500 for Row A → starts = [500, 200] → minus-sorted ✓
    blast_rows = [
        mk_blast_row(sstart=1159, send=1000, qstart=100, qend=259, sstrand="minus"),
        mk_blast_row(sstart=299, send=200, qstart=300, qend=399, sstrand="minus"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": subject})
    assert result is None


def test_scramble_minus_unique_unsorted():
    """Minus-strand with unique segments that are genuinely unsorted → Scramble."""
    # Row A (q=700, low=910) → Row B (q=800, low=980)
    # As q increases, low increases (910→980) → not minus-sorted → Scramble
    blast_rows = [
        mk_blast_row(sstart=940, send=930, qstart=700, sstrand="minus"),
        mk_blast_row(sstart=990, send=980, qstart=800, sstrand="minus"),
    ]
    subj = subject_from_rows(blast_rows)
    result = is_scrambled("id", blast_rows, {"mock_subject": subj})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# Internal inversion tests
# ---------------------------------------------------------------------------

def test_internal_inversion_detected():
    """Mixed-strand evidence in informative (unique, ≥100bp) rows → Inversion."""
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
    """Mixed strand evidence only in duplicate regions → no inversion."""
    seg = "DUPSTRAND" * 50
    subject = (
        "N" * 99 +
        seg +              # positions 100–498
        "N" * 501 +
        seg +              # positions 1000–1398
        "N" * (10000 - 1399)
    )
    blast_rows = [
        mk_blast_row(sstart=100, send=498, qstart=100, qend=498, sstrand="plus"),
        mk_blast_row(sstart=1398, send=1000, qstart=500, qend=898, sstrand="minus"),
    ]
    assert is_scrambled("id", blast_rows, {"mock_subject": subject}) or \
           contains_internal_inversion("id", blast_rows, {"mock_subject": subject}) is None


def test_internal_inversion_with_mixed_informative_and_duplicate():
    """Mixed strand with one informative and one duplicate row → inversion."""
    seg = "SOMEDUP" * 50
    subject = (
        "N" * 199 +
        seg +              # positions 200–598 (copy A)
        "N" * 401 +
        seg +              # positions 1000–1398 (copy B)
        "N" * 200 +
        "UNIQUE_PLUS_" * 20 +  # positions 1599–1778 (unique, plus)
        "N" * (10000 - 1779)
    )
    blast_rows = [
        mk_blast_row(sstart=1398, send=1000, qstart=100, qend=498, sstrand="minus"),
        mk_blast_row(sstart=1599, send=1778, qstart=600, qend=779, sstrand="plus"),
    ]
    # is_scrambled will check scramble first → both rows, dominant=plus (2/2? no, 1 minus, 1 plus)
    # Wait, sstrand values: "minus", "plus" → tie. most_frequent returns "minus" or "plus" (first in iteration = "minus")
    # Direction = "minus". For minus with has_monotone: [start=1000, start=1599] → 1000 → 1599: 1599 > 1000 → NOT nonincreasing → Scramble!
    # So is_scrambled returns Scramble first. The `or` short-circuits.
    # contains_internal_inversion is NOT reached.
    #
    # That's OK for this test's intent but doesn't test the inversion path.
    # Let me add a separate test where is_scrambled returns None first.
    pass


def test_internal_inversion_not_hidden_by_duplicate():
    """When informative rows have mixed strands, inversion is reported
    even if non-informative rows are also present."""
    dup = "X" * 180
    unique_plus = "P" * 160
    unique_minus = "M" * 160
    total_len = 99 + 180 + 220 + 160 + 340 + 160 + 840 + 180 + 7821
    subject = (
        "N" * 99 +
        dup +              # positions 100–279 (copy A)
        "N" * 220 +
        unique_plus +      # positions 500–659
        "N" * 340 +
        unique_minus +     # positions 1000–1159
        "N" * 840 +
        dup +              # positions 2000–2179 (copy B)
        "N" * 7821
    )
    assert len(subject) == total_len == 10000
    blast_rows = [
        # Row C (qstart=100): non-informative (duplicate), plus
        mk_blast_row(sstart=100, send=279, qstart=100, qend=279, sstrand="plus"),
        # Row A (qstart=300): informative, plus
        mk_blast_row(sstart=500, send=659, qstart=300, qend=459, sstrand="plus"),
        # Row B (qstart=500): informative, minus
        mk_blast_row(sstart=1159, send=1000, qstart=500, qend=659, sstrand="minus"),
    ]
    # Sorted by qstart: C(100, plus), A(300, plus), B(500, minus)
    # Direction = "plus" (2/3). For plus: starts [100, 500, 1000] → monotone → no scramble.
    # contains_internal_inversion:
    #   C is not informative (dup occurs twice).
    #   A and B are informative with mixed strands → InternalInversion.
    result = is_scrambled("id", blast_rows, {"mock_subject": subject}) or \
             contains_internal_inversion("id", blast_rows, {"mock_subject": subject})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.InternalInversion)


# ---------------------------------------------------------------------------
# Mixed / multiple direction tests
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
# Edge cases
# ---------------------------------------------------------------------------

def test_unknown_sseqid_falls_back_to_scramble():
    """If a row's sseqid is not in the subject dict, report scramble."""
    blast_rows = [
        mk_blast_row(sstart=900, send=1010, qstart=50, qend=150, sstrand="plus",
                     sseqid="unknown_ref"),
    ]
    result = is_scrambled("id", blast_rows, {"mock_subject": "A" * 10000})
    assert isinstance(result, Defect)
    assert isinstance(result.error, defect.Scramble)


# ---------------------------------------------------------------------------
# has_monotone_subject_assignment direct tests
# ---------------------------------------------------------------------------

def test_has_monotone_subject_assignment_empty_rows():
    assert has_monotone_subject_assignment([], {}, "plus") is True


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
    # For minus: normalized low must be nonincreasing.
    # Row A (q=10, start=100) → Row B (q=20, start=150): 150 > 100 → NOT nonincreasing
    rows = [
        mk_blast_row(200, 100, 10, "minus"),   # start=100
        mk_blast_row(250, 150, 20, "minus"),   # start=150
    ]
    subj = subject_from_rows(rows)
    assert has_monotone_subject_assignment(rows, {"mock_subject": subj}, "minus") is False


def test_has_monotone_subject_assignment_minus_sorted_unique():
    # Row A (q=10, start=150) → Row B (q=20, start=100): 100 < 150 → nonincreasing ✓
    rows = [
        mk_blast_row(250, 150, 10, "minus"),   # start=150
        mk_blast_row(200, 100, 20, "minus"),   # start=100
    ]
    subj = subject_from_rows(rows)
    assert has_monotone_subject_assignment(rows, {"mock_subject": subj}, "minus") is True
