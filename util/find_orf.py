from util.reference_index import ReferenceIndex
from util.get_query_aminoacids_table import get_query_aminoacids_table
from util.get_biggest_protein import get_biggest_protein
from util.candidate_orf import CandidateORF

import util.detailed_aligner as detailed_aligner

def has_start_codon(orf):
    return orf.aminoacids.startswith("M")


def has_stop_codon(orf):
    return orf.aminoacids.endswith("*")


def find_closest(aminoacids, start, direction, target):
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


def find_candidate_positions(aligned_sequence, e):
    q_start = ReferenceIndex(e.start).mapto(aligned_sequence)
    q_end = ReferenceIndex(e.end).mapto(aligned_sequence)
    expected_aminoacids = e.aminoacids
    expected_protein = expected_aminoacids.strip("*")
    q_start_a = q_start // 3
    q_end_a = q_end // 3
    n = len(aligned_sequence.this.seq) - 1
    visited_set = set()
    query_aminoacids_table = get_query_aminoacids_table(aligned_sequence.this)
    if isinstance(query_aminoacids_table, Exception):
        raise query_aminoacids_table

    for frame in range(3):
        aminoacids = query_aminoacids_table[frame]
        for start_direction in [-1, +1]:
            for end_direction in [-1, +1]:
                closest_start_a = q_start_a if not has_start_codon(e) else find_closest(aminoacids, q_start_a, start_direction, 'M')
                closest_end_a = q_end_a if not has_stop_codon(e) else find_closest(aminoacids, q_end_a, end_direction, '*')
                got_aminoacids = aminoacids[closest_start_a:closest_end_a + 1]
                if got_aminoacids in visited_set:
                    continue
                else:
                    visited_set.add(got_aminoacids)

                closest_start = min(n, (closest_start_a * 3) + frame)
                closest_end = min(n + 1, (closest_end_a * 3) + 3 + frame)
                got_protein = get_biggest_protein(has_start_codon(e), got_aminoacids)
                dist = detailed_aligner.align(got_protein, expected_protein).distance()
                yield CandidateORF(e.name, closest_start, closest_end, e.start, e.end,
                                    "forward", dist, got_protein, got_aminoacids)


def find_orf(aligned_sequence, e):
    candidates = find_candidate_positions(aligned_sequence, e)
    return min(candidates, key=lambda x: x.distance)
