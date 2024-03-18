from typing import Iterable
from math import floor

from cfeintact.reference_index import ReferenceIndex
from cfeintact.get_query_aminoacids_table import get_query_aminoacids_table
from cfeintact.get_biggest_protein import get_biggest_protein
from cfeintact.original_orf import OriginalORF
from cfeintact.mapped_orf import MappedORF
from cfeintact.aligned_sequence import AlignedSequence
from cfeintact.translate_to_aminoacids import translate_to_aminoacids

import cfeintact.detailed_aligner as detailed_aligner


def has_start_codon(orf):
    return orf.aminoacids.startswith("M")


def has_stop_codon(orf):
    return orf.aminoacids.endswith("*")


def find_closest(aminoacids, start, direction, target):
    assert isinstance(start, int) or start.is_integer()
    start = int(start)

    distance = 0
    n = len(aminoacids) - 1

    while start + distance >= 0 and start + distance <= n:
        if aminoacids[start + distance] == target:
            return start + distance
        distance += direction

    if target == '*':
        return n
    else:
        return 0


def find_candidate_positions(aligned_sequence: AlignedSequence, e: OriginalORF) -> Iterable[MappedORF]:
    q_start = ReferenceIndex(e.start).mapto(aligned_sequence)
    q_end = ReferenceIndex(e.end).mapto(aligned_sequence)
    expected_aminoacids = e.aminoacids
    expected_protein = expected_aminoacids.strip("*")
    visited_set = set()
    query_aminoacids_table = get_query_aminoacids_table(aligned_sequence.this)

    for frame in range(3):
        aminoacids = query_aminoacids_table[frame]
        for start_direction in [-1, +1]:
            for end_direction in [-1, +1]:
                if has_start_codon(e):
                    q_start_a = floor((q_start - frame) / 3)
                    closest_start_a = find_closest(aminoacids, q_start_a, start_direction, 'M')
                    closest_start = (closest_start_a * 3) + frame
                else:
                    closest_start = q_start
                    closest_start_a = floor(closest_start / 3)

                if has_stop_codon(e):
                    q_end_a = floor((q_end - frame) / 3)
                    closest_end_a = find_closest(aminoacids, q_end_a, end_direction, '*')
                    closest_end = (closest_end_a * 3) + 3 + frame - 1
                else:
                    closest_end = q_end
                    closest_end_a = floor(closest_end / 3) + 1

                if (closest_start, closest_end) in visited_set:
                    continue
                else:
                    visited_set.add((closest_start, closest_end))

                got_nucleotides = str(aligned_sequence.slice(closest_start, closest_end).this.seq)
                got_aminoacids = translate_to_aminoacids(got_nucleotides)
                got_protein = get_biggest_protein(has_start_codon(e), got_aminoacids)
                if has_start_codon(e) and has_stop_codon(e):
                    dist = detailed_aligner.align(got_protein, expected_protein).distance()
                else:
                    dist = detailed_aligner.align(got_aminoacids, expected_aminoacids).distance()

                orf = OriginalORF(
                    name=e.name,
                    start=closest_start,
                    end=closest_end,
                    nucleotides=got_nucleotides,
                    aminoacids=got_aminoacids,
                    protein=got_protein,
                    deletion_tolerence=e.deletion_tolerence,
                    is_small=e.is_small,
                )
                yield MappedORF(
                    reference=e,
                    query=orf,
                    orientation="forward",
                    distance=dist,
                )


def find_orf(aligned_sequence: AlignedSequence, e: OriginalORF) -> MappedORF:
    candidates = find_candidate_positions(aligned_sequence, e)
    return min(candidates, key=lambda x: x.distance)
