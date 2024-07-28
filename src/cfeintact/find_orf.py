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
    q_start = aligned_sequence.coordinate_mapping.ref_to_query.right_min(e.start)
    q_end = aligned_sequence.coordinate_mapping.ref_to_query.left_max(e.end)
    if q_start is None:
        q_start = 0
    if q_end is None:
        q_end = len(aligned_sequence.this.seq) - 1

    region_nucleotides = str(aligned_sequence.this.seq[q_start:q_end + 1])
    region_aminoacids = translate_to_aminoacids(region_nucleotides)

    visited_set = set()
    query_aminoacids_table = get_query_aminoacids_table(aligned_sequence.this)
    initial_frame = q_start % 3

    for frame_shift in range(3):
        q_start_in_aa = (q_start - initial_frame) / 3
        q_end_in_aa = (q_end - 2 - initial_frame) / 3
        q_start_in_aa = int(q_start_in_aa)
        q_end_in_aa = floor(q_end_in_aa)
        frame = (initial_frame + frame_shift) % 3
        aminoacids = query_aminoacids_table[frame]

        if e.has_start_codon:
            closest_start_in_aa = find_closest(aminoacids, q_start_in_aa, -1, 'M')
            closest_start = (closest_start_in_aa * 3) + frame
        else:
            closest_start = q_start
            closest_start_in_aa = floor(closest_start / 3)

        if e.has_stop_codon:
            closest_end_in_aa = find_closest(aminoacids, q_end_in_aa, +1, '*')
            closest_end = (closest_end_in_aa * 3) + 3 + frame - 1
        else:
            closest_end = q_end
            closest_end_in_aa = floor(closest_end / 3) + 1

        if (closest_start, closest_end) in visited_set:
            continue
        else:
            visited_set.add((closest_start, closest_end))

        protein_nucleotides = str(aligned_sequence.this.seq[closest_start:closest_end + 1])
        protein_aminoacids = translate_to_aminoacids(protein_nucleotides)
        found_protein = get_biggest_protein(e.has_start_codon, protein_aminoacids)
        if found_protein:
            protein, protein_aa_start, protein_aa_end = found_protein
        else:
            protein = ""

        orf = OriginalORF(
            name=e.name,
            start=q_start,
            end=q_end,
            nucleotides=region_nucleotides,
            aminoacids=region_aminoacids,
            protein=protein,
            max_deletions=e.max_deletions,
            max_insertions=e.max_insertions,
            max_distance=e.max_distance,
            is_small=e.is_small,
        )

        mapped = MappedORF(
            reference=e,
            query=orf,
            orientation="forward",
        )

        yield mapped


FIND_ORF_CACHE: Dict[Tuple[str, str, str], MappedORF] = {}

def find_orf(aligned_sequence: AlignedSequence, e: OriginalORF) -> MappedORF:
    assert aligned_sequence.this.id is not None
    assert aligned_sequence.reference.id is not None
    key: Tuple[str, str, str] = (aligned_sequence.this.id, aligned_sequence.reference.id, e.name)

    if key not in FIND_ORF_CACHE:
        candidates = find_candidate_positions(aligned_sequence, e)
        result = min(candidates, key=lambda x: x.distance)
        FIND_ORF_CACHE[key] = result

    return FIND_ORF_CACHE[key]
