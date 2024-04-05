from typing import Iterable, Dict, Tuple
from math import floor

from cfeintact.get_query_aminoacids_table import get_query_aminoacids_table
from cfeintact.get_biggest_protein import get_biggest_protein
from cfeintact.original_orf import OriginalORF
from cfeintact.mapped_orf import MappedORF
from cfeintact.aligned_sequence import AlignedSequence
from cfeintact.translate_to_aminoacids import translate_to_aminoacids


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
    q_start = aligned_sequence.coordinate_mapping.ref_to_query.right_min(e.start) or 0
    q_end = aligned_sequence.coordinate_mapping.ref_to_query.left_max(e.end) or 0
    assert q_start is not None
    assert q_end is not None
    visited_set = set()
    query_aminoacids_table = get_query_aminoacids_table(aligned_sequence.this)

    for frame in range(3):
        aminoacids = query_aminoacids_table[frame]
        for start_direction in [-1, +1]:
            for end_direction in [-1, +1]:
                if e.has_start_codon:
                    q_start_a = floor((q_start - frame) / 3)
                    closest_start_a = find_closest(aminoacids, q_start_a, start_direction, 'M')
                    closest_start = (closest_start_a * 3) + frame
                else:
                    closest_start = q_start
                    closest_start_a = floor(closest_start / 3)

                if e.has_stop_codon:
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

                got_nucleotides = str(aligned_sequence.this.seq[closest_start:closest_end + 1])
                got_aminoacids = translate_to_aminoacids(got_nucleotides)
                got_protein = get_biggest_protein(e.has_start_codon, got_aminoacids)

                orf = OriginalORF(
                    name=e.name,
                    start=closest_start,
                    end=closest_end,
                    nucleotides=got_nucleotides,
                    aminoacids=got_aminoacids,
                    protein=got_protein,
                    max_deletions=e.max_deletions,
                    max_insertions=e.max_insertions,
                    max_distance=e.max_distance,
                    is_small=e.is_small,
                )

                yield MappedORF(
                    reference=e,
                    query=orf,
                    orientation="forward",
                )


FIND_ORF_CACHE: Dict[Tuple[str, str, str], MappedORF] = {}

def find_orf(aligned_sequence: AlignedSequence, e: OriginalORF) -> MappedORF:
    assert aligned_sequence.this.id is not None
    assert aligned_sequence.reference.id is not None
    key: Tuple[str, str, str] = (aligned_sequence.this.id, aligned_sequence.reference.id, e.name)

    candidates = find_candidate_positions(aligned_sequence, e)
    result = min(candidates, key=lambda x: x.distance)

    if key not in FIND_ORF_CACHE:
        FIND_ORF_CACHE[key] = result

    return FIND_ORF_CACHE[key]
