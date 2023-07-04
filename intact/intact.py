import json
import os
import re
import subprocess
import sys
import uuid
import tempfile
import csv
from dataclasses import dataclass
from collections import Counter
import Bio
from Bio import AlignIO, Seq, SeqIO, SeqRecord, Align
from scipy.stats import fisher_exact
from jarowinkler import jaro_similarity

import util.constants as const
import util.subtypes as st
import util.wrappers as wrappers
import util.log as log


WRONGORFNUMBER_ERROR = "WrongORFNumber"
MISPLACEDORF_ERROR   = "MisplacedORF"
LONGDELETION_ERROR   = "LongDeletion"
DELETIONINORF_ERROR  = "DeletionInOrf"
INSERTIONINORF_ERROR  = "InsertionInOrf"
INTERNALSTOP_ERROR   = "InternalStopInOrf"
SCRAMBLE_ERROR       = "Scramble"
NONHIV_ERROR         = "NonHIV"
INTERNALINVERSION_ERROR = "InternalInversion"

FRAMESHIFTINORF_ERROR  = "FrameshiftInOrf"
MSDMUTATED_ERROR = "MajorSpliceDonorSiteMutated"
PSIDELETION_ERROR    = "PackagingSignalDeletion"
PSINOTFOUND_ERROR    = "PackagingSignalNotComplete"
RREDELETION_ERROR    = "RevResponseElementDeletion"
HYPERMUTATION_ERROR  = "APOBECHypermutationDetected"


@dataclass
class IntactnessError:
    sequence_name: str
    error: str
    message: str

@dataclass
class ORF:
    name: str
    orientation: str
    start: int
    end: int
    distance: str
    aminoacids: str
    nucleotides: str

@dataclass
class ExpectedORF:
    name: str
    start: int
    end: int
    deletion_tolerence: int

    @staticmethod
    def subtyped(pos_mapping, name, start, end, deletion_tolerence):
        vpr_defective_insertion_pos = 5772
        start = start if start < vpr_defective_insertion_pos else start - 1
        end = end if end < vpr_defective_insertion_pos else end - 1

        start_s = pos_mapping[start - 1]
        end_s = pos_mapping[end]
        return ExpectedORF(name, start_s, end_s, deletion_tolerence)

@dataclass
class ReceivedORF:
    start: int
    end: int
    deleted_count: int
    inserted_count: int

@dataclass
class CandidateORF:
    name: str
    start: int
    end: int
    orientation: str
    distance: float
    aminoseq: str
    expectedaminoseq: str

@dataclass
class BlastRow:
    qseqid: str
    qlen: int
    sseqid: str
    sgi: str
    slen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    length: int
    pident: float
    nident: float
    btop: int
    stitle: str
    sstrand: str

    @staticmethod
    def init(row):
        it = iter(row)
        return BlastRow(
            qseqid=next(it),
            qlen=int(next(it)),
            sseqid=next(it),
            sgi=next(it),
            slen=int(next(it)),
            qstart=int(next(it)),
            qend=int(next(it)),
            sstart=int(next(it)),
            send=int(next(it)),
            evalue=float(next(it)),
            bitscore=float(next(it)),
            length=int(next(it)),
            pident=float(next(it)),
            nident=float(next(it)),
            btop=next(it),
            stitle=next(it),
            sstrand=next(it),
        )


def iterate_values_from_tsv(file_path):
    with open(file_path, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            yield row


def iterate_blast_rows_from_tsv(file_path):
    previous_key = None
    values = []

    for row in iterate_values_from_tsv(file_path):
        key = row[0]
        typed = BlastRow.init(row)

        if key != previous_key and previous_key is not None:
            yield values
            values = []

        values.append(typed)
        previous_key = key

    if values:
        yield values


def blast_iterate(subtype, input_file):
    with tempfile.NamedTemporaryFile() as output_file:
        db_file = st.alignment_file(subtype)
        wrappers.blast(db_file, input_file, output_file.name)
        for seq in iterate_blast_rows_from_tsv(output_file.name):
            yield seq


def blast_iterate_inf(subtype, input_file):
    for seq in blast_iterate(subtype, input_file):
        yield seq

    while True:
        yield []


def is_sorted(lst):
    last = None
    for x in lst:
        if last is not None and x < last:
            return False
        last = x
    return True


def most_frequent_element(lst):
    counter = Counter(lst)
    most_common = counter.most_common(1)
    return most_common[0][0] if most_common else None


def remove_5_prime(blast_rows):
    # HIV 5' region can easily map to its 3' region because they are identical.
    # Such a maping would not constitute a scramble, so we ignore the 5' region for this check.
    return [x for x in blast_rows if x.sstart > 622 and x.send > 622]


def contains_internal_inversion(seqid, blast_rows):
    ignored_5_prime = remove_5_prime(blast_rows)

    if not ignored_5_prime:
        # No alignment.
        # It should be an error normally, yet not an internal inversion error.
        return None

    all_same = len(set(x.sstrand for x in ignored_5_prime)) == 1
    if not all_same:
        # Some parts of the sequence were aligned
        # in forward direction (plus)
        # and some in reverse (minus).
        # This indicates an internal inversion.
        return IntactnessError(seqid, INTERNALINVERSION_ERROR,
                               "Sequence contains an internal inversion.")
    else:
        return None


def is_scrambled(seqid, blast_rows):
    ignored_5_prime = remove_5_prime(blast_rows)

    if not ignored_5_prime:
        # No alignment.
        # It should be an error normally, yet not a scramble error.
        return None

    ignored_5_prime.sort(key=lambda x: x.qstart)
    direction = most_frequent_element(x.sstrand for x in ignored_5_prime)
    if direction == "plus" and is_sorted(x.sstart for x in ignored_5_prime):
        return None
    elif direction == "minus" and is_sorted(x.send for x in reversed(ignored_5_prime)):
        return None
    else:
        return IntactnessError(seqid, SCRAMBLE_ERROR,
                               f"Sequence is {direction}-scrambled.")


def is_nonhiv(seqid, seqlen, blast_rows):
    aligned_length = sum(abs(x.qend - x.qstart) + 1 for x in blast_rows)
    total_length = blast_rows[0].qlen if blast_rows else 1
    ratio = aligned_length / total_length

    if ratio < 0.8 and seqlen > total_length * 0.6:
        return IntactnessError(seqid, NONHIV_ERROR,
                               "Sequence contains unrecognized parts. "
                               "It is probably a Human/HIV Chimera sequence.")
    else:
        return None

def _getPositions(pattern, string):
    """
    Hidden function used to get interator of all intances of pattern in string

    Keyword Args:

        pattern -- regex pattern to search for.
        string -- input string to search in.

    Return:
        set of matched indices

    """
    return set(m.start() for m in re.finditer(pattern, string))
#/end _getPositions

def isHypermut(aln):
    """
    APOBEC3G/F hypermutation scan and test based on Rose and Korber, Bioinformatics (2000).
    Briefly, scans reference for APOBEC possible signatures and non-signatures and performs
    fisher test based on ratio of G->A in ref -> query at these signatures.

    Keyword Args:
        aln -- AlignIO object where aln[0] is the reference and aln[1] is the query.

    Internal Args:
        newAln -- all non-indel positions (i.e. removes all sites 
            where a '-' exists in either sequence.
        ref_str -- joined string from newAln[0]
        qry_str -- joined string from newAln[1]

        mutIndex -- APOBEC signature in qry_str
        ctlIndex -- non-APOBEC signature in qry_str
        gIndex -- all 'G' positions in ref
        aIndex -- all 'A' positions in qry
"""
    mutIndex = _getPositions('(?=.[AG][^C])', str(aln[1].seq).upper())
    ctlIndex = _getPositions('(?=.([CT].|[AG]C))', str(aln[1].seq).upper())
    gIndex = _getPositions('G', str(aln[0].seq).upper())
    aIndex = _getPositions('A', str(aln[1].seq).upper())

    mutPoss = gIndex.intersection(mutIndex)
    ctlPoss = gIndex.intersection(ctlIndex)

    nMut = len(mutPoss.intersection(aIndex))
    nCtl = len(ctlPoss.intersection(aIndex))

    mutTot = len(gIndex.intersection(mutIndex))
    ctlTot = len(gIndex.intersection(ctlIndex))

    _, pval = fisher_exact(
            [[nMut, mutTot - nMut],
             [nCtl, ctlTot - nCtl]],
            alternative = 'greater'
        )

    if pval < 0.05:
        return IntactnessError(
                aln[1].id, HYPERMUTATION_ERROR,
                "Query sequence shows evidence of APOBEC3F/G-mediated" +
                " hypermutation (p = " + str(pval) + ")."
        )
    #/end if
    return None
#/end isHypermut


def has_long_deletion(sequence, alignment):
    """
    Determines whether the sequence has a long deletion in it.
    Keyword Args:

        sequence -- the query sequence.
        alignment -- multiple sequence alignment object containing the
                     reference and query sequence.
    """
    # NOTE: This is the same check that HIVSeqInR uses.
    if len(sequence.seq) < 8000:
        return IntactnessError(sequence.id,
                               LONGDELETION_ERROR,
                               "Query sequence contains a long deletion.")
    return None
#/end has_long_deletion


def has_mutated_major_splice_donor_site(alignment, 
                                        splice_donor_start_pos, 
                                        splice_donor_end_pos,
                                        splice_donor_sequence):
    """
    Determines whether the major splice donor site is mutated.
    Keyword Args:
        
        alignment -- multiple sequence alignment object containing the 
                     reference and query sequence.
        splice_donor_start_pos -- first position of splice donor site
        splice_donor_end_pos -- last position of splice donor site
        splice_donor_sequence - sequence of splice donor site
    """

    sd_begin = [m.start() for m in re.finditer(r"[^-]",
                str(alignment[0].seq))][splice_donor_start_pos]
    
    sd_end = [m.start() for m in re.finditer(r"[^-]",
              str(alignment[0].seq))][splice_donor_end_pos]
    
    sd = alignment[1].seq[sd_begin:(sd_end + 1)]

    # splice donor site is missing from sequence
    if all([x == "-" for x in sd]):
        return IntactnessError(
                alignment[1].id, MSDMUTATED_ERROR,
                "Query sequence has a missing splice donor site, " 
                + "".join(sd.upper()) + "."
                )

    if sd.upper() != splice_donor_sequence.upper():

        return IntactnessError(
                alignment[1].id, MSDMUTATED_ERROR,
                "Query sequence has a mutated splice donor site, " 
                + "".join(sd.upper()) + "."
                )

    return None

    

    
    

def has_packaging_signal(alignment, psi_locus, psi_tolerance):
    """
    Determines presence and possible intactness of HIV 
    Packaging Signal Region.
    
    
    Keyword Args:
        
        alignment -- multiple sequence alignment object containing the 
                     reference and query sequence.
        psi_locus -- tuple containing start and end coordinates of PSI wrt
                     the reference being used.
        psi_tolerance -- number of deletions in query tolerated to be intact
        
        
    Internal Args:
        
        packaging_begin -- Aligned PSI start position.
        packaging_end -- Aligned PSI end position.
        query_start -- beginning position of query sequence.
        query_psi -- extracted PSI region from query
        query_psi_deletions -- number of deletions in query PSI region
        
        
    Return:
        
        PSINOTFOUND_ERROR -- IntactnessError denoting query does not encompass
                             the complete PSI region.
        PSIDELETION_ERROR -- IntactnessError denoting query likely contains
                             defective PSI.
        None -- Denotes intact PSI.    
    """
    packaging_begin = [m.start() for m in re.finditer(r"[^-]",
                       str(alignment[0].seq))][psi_locus[0]]
    query_start = [m.start() for m in re.finditer(r"[^-]",
                   str(alignment[1].seq))][0]
    packaging_end = [m.start() for m in re.finditer(r"[^-]",
                     str(alignment[0].seq))][psi_locus[1]]
    # if query_start > packaging_begin:
    #     return IntactnessError(
    #             alignment[1].id, PSINOTFOUND_ERROR,
    #             "Query Start at reference position " + str(query_start)
    #             + ". Does not encompass PSI at positions "
    #             + str(packaging_begin) + " to " + str(packaging_end) + "."
    #             )
    # #/end if
    query_psi = str(alignment[1].seq[packaging_begin:packaging_end])
    query_psi_deletions = len(re.findall(r"-", query_psi))
    if query_psi_deletions > psi_tolerance:
        return IntactnessError(
                alignment[1].id, PSIDELETION_ERROR,
                "Query Sequence exceeds maximum deletion tolerance in PSI. " + 
                "Contains " + str(query_psi_deletions) + " deletions with max "
                + "tolerance of " + str(psi_tolerance) + " deletions."
                )
    #/end if
    return None
#/end def has_packaging_signal

def has_rev_response_element(alignment, rre_locus, rre_tolerance):
    """
    Determines presence and possible intactness of HIV 
    Packaging Signal Region.
    
    
    Keyword Args:
        
        alignment -- multiple sequence alignment object containing the 
                     reference and query sequence.
        rre_locus -- tuple containing start and end coordinates of RRE wrt
                     the reference being used.
        RRE_tolerance -- number of deletions in query tolerated to be intact
        
        
    Internal Args:
        
        rre_begin -- Aligned RRE start position.
        rre_end -- Aligned RRE end position.
        query_rre -- extracted RRE region from query
        query_rre_deletions -- number of deletions in query RRE region
        
        
    Return:
        
        RREDELETION_ERROR -- IntactnessError denoting query likely contains
                             defective RRE.
        None -- Denotes intact RRE.    
    """
    rre_begin = [m.start() for m in re.finditer(r"[^-]",
                       str(alignment[0].seq))][rre_locus[0]]
    rre_end = [m.start() for m in re.finditer(r"[^-]",
                     str(alignment[0].seq))][rre_locus[1]]
    query_rre = str(alignment[1].seq[rre_begin:rre_end])
    query_rre_deletions = len(re.findall(r"-", query_rre))
    if query_rre_deletions > rre_tolerance:
        return IntactnessError(
                alignment[1].id, RREDELETION_ERROR,
                "Query Sequence exceeds maximum deletion tolerance in RRE. " + 
                "Contains " + str(query_rre_deletions) + " deletions with max "
                + "tolerance of " + str(rre_tolerance) + " deletions."
                )
    #/end if
    return None
#/end def has_rev_response_element


def alignment_score(alignment):
    """
    Simple score for an alignment (just to find out if it's forward or
    reverse - absolutely not a true genetic distance)
    """

    return sum([a==b for a, b in zip(alignment[0].seq, alignment[1].seq)])

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2

def has_reading_frames(
    alignment, sequence, reference, is_small,
    expected, error_bar, reverse = False
):
    """
    Check for presence of small reading frames
    """

    import util.coordinates as coords
    coordinates_mapping = coords.map_positions(alignment[0], alignment[1].seq)
    reverse_coordinates_mapping = coords.map_positions(alignment[1], alignment[0].seq)

    reference_aligned_mapping = coords.map_nonaligned_to_aligned_positions(reference, alignment[0].seq)
    query_aligned_mapping = coords.map_nonaligned_to_aligned_positions(sequence, alignment[1].seq)

    def translate(seq, frame = 0, to_stop = False):
        for_translation = seq[frame:]
        for_translation += 'N' * ({0: 0, 1: 2, 2: 1}[len(for_translation) % 3])
        return Seq.translate(for_translation, to_stop = to_stop)

    try:
        query_aminoacids_table = [translate(sequence.seq, i) for i in range(3)]
    except Bio.Data.CodonTable.TranslationError as e:
        log.error(e)
        return [], []

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

    def find_candidate_positions(e, q_start, q_end):
        expected_nucleotides = str(reference.seq[e.start:e.end])
        expected_aminoacids = translate(expected_nucleotides)
        q_start = coordinates_mapping[e.start]
        q_end = coordinates_mapping[e.end]
        got_nucleotides = sequence.seq[q_start:q_end]
        got_aminoacids = translate(got_nucleotides)
        q_start_a = q_start // 3
        q_end_a = q_end // 3
        has_start_codon = expected_aminoacids[0] == 'M'
        has_stop_codon = expected_aminoacids[-1] == '*'
        n = len(sequence.seq) - 1

        for frame in range(3):
            aminoacids = query_aminoacids_table[frame]
            for start_direction in [-1, +1]:
                for end_direction in [-1, +1]:
                    closest_start_a = q_start_a if not has_start_codon else find_closest(aminoacids, q_start_a, start_direction, 'M')
                    closest_end_a = q_end_a if not has_stop_codon else find_closest(aminoacids, q_end_a, end_direction, '*')
                    got_aminoacids = aminoacids[closest_start_a:closest_end_a + 1]
                    dist = 1 - jaro_similarity(got_aminoacids, expected_aminoacids)
                    closest_start = min(n, (closest_start_a * 3) + frame)
                    closest_end = min(n + 1, (closest_end_a * 3) + 3 + frame)
                    yield CandidateORF(e.name, closest_start, closest_end, "forward",
                                       dist, got_aminoacids, expected_aminoacids)

    def find_real_correspondence(e):
        q_start = coordinates_mapping[e.start]
        q_end = coordinates_mapping[e.end]
        candidates = list(find_candidate_positions(e, q_start, q_end))
        return min(candidates, key=lambda x: x.distance)

    def get_indel_impact(alignment):
        shift = 0
        impacted = 0
        counter = 0
        good = 0
        for (x, y) in zip(alignment[0], alignment[1]):
            if x == "-" and y == "-":
                continue
            if x == "-" and y != "-":
                shift += 1
            if x != "-" and y == "-":
                shift -= 1
            if shift % 3 != 0:
                counter += 1
                good = 0
            else:
                good += 1
                if good > 5:
                    impacted += counter if counter > 30 else 0
                    counter = 0

        impacted += counter if counter > 30 else 0
        return impacted

    errors = []
    matches = []
    for e in expected:
        best_match = find_real_correspondence(e)
        matches.append(best_match)

        got_protein = best_match.aminoseq.split("*")[0]
        exp_protein = best_match.expectedaminoseq.split("*")[0]

        got_nucleotides = sequence.seq[best_match.start:best_match.start + len(got_protein) * 3].upper()
        exp_nucleotides = reference.seq[e.start:e.end].upper()
        if got_nucleotides:
            orf_alignment = aligner.align(exp_nucleotides, got_nucleotides)[0]
            best_match.distance = -1 * (orf_alignment.score / len(exp_nucleotides))
        else:
            orf_alignment = (exp_nucleotides, "-" * len(exp_nucleotides))
            best_match.distance = len(exp_nucleotides) ** 2

        deletions = max(0, len(exp_protein) - len(got_protein)) * 3
        insertions = max(0, len(got_protein) - len(exp_protein)) * 3

        # Max deletion allowed in ORF exceeded
        if deletions > e.deletion_tolerence:

            if "*" in best_match.aminoseq[:(-e.deletion_tolerence // 3)]:
                errors.append(IntactnessError(
                    sequence.id, INTERNALSTOP_ERROR,
                    ("Smaller " if is_small else "")
                    + "ORF " + str(e.name) + " at " + str(e.start)
                    + "-" + str(e.end)
                    + " contains an internal stop codon"
                ))
            else:
                errors.append(IntactnessError(
                    sequence.id, DELETIONINORF_ERROR,
                    ("Smaller " if is_small else "")
                    + "ORF " + str(e.name) + " at " + str(e.start)
                    + "-" + str(e.end)
                    + " can have maximum deletions "
                    + str(e.deletion_tolerence) + ", got "
                    + str(deletions)
                ))

            continue

        # Max insertions allowed in ORF exceeded
        if insertions > 3 * e.deletion_tolerence:

            errors.append(IntactnessError(
                sequence.id, INSERTIONINORF_ERROR,
                ("Smaller " if is_small else "")
                + "ORF " + str(e.name) + " at " + str(e.start)
                + "-" + str(e.end)
                + " can have maximum insertions "
                + str(3 * e.deletion_tolerence) + ", got "
                + str(insertions)
            ))

            continue

        impacted_by_indels = get_indel_impact(orf_alignment)

        # Check for frameshift in ORF
        if impacted_by_indels >= e.deletion_tolerence + 0.10 * len(exp_nucleotides):

            errors.append(IntactnessError(
                sequence.id, FRAMESHIFTINORF_ERROR,
                ("Smaller " if is_small else "")
                + "ORF " + str(e.name) + " at " + str(e.start) 
                + "-" + str(e.end) 
                + " contains out of frame indels that impact " + str(impacted_by_indels)
                + " positions."
            ))

            continue


    return matches, errors


def iterate_sequences(input_file):
    with open(input_file, 'r') as in_handle:
        for sequence in SeqIO.parse(in_handle, "fasta"):
            yield sequence


def iterate_empty_lists():
    while True:
        yield []


def with_blast_rows(blast_it, sequence_it):
    blast_rows = next(blast_it)
    for sequence in sequence_it:
        if blast_rows and blast_rows[0].qseqid == sequence.id:
            yield (sequence, blast_rows)
            blast_rows = next(blast_it)
        else:
            yield (sequence, [])


def intact( working_dir,
            input_file,
            subtype,
            include_packaging_signal,
            include_rre,
            check_major_splice_donor_site,
            run_hypermut,
            check_long_deletion,
            check_nonhiv,
            check_scramble,
            check_internal_inversion,
            include_small_orfs,
            hxb2_forward_orfs = const.DEFAULT_FORWARD_ORFs,
            hxb2_reverse_orfs = const.DEFAULT_REVERSE_ORFS,
            hxb2_small_orfs = const.DEFAULT_SMALL_FORWARD_ORFS,
            hxb2_psi_locus = const.DEFAULT_PSI_LOCUS,
            hxb2_rre_locus = const.DEFAULT_RRE_LOCUS,
            hxb2_msd_site_locus = const.DEFAULT_MSD_SITE_LOCUS,
            min_orf_length = const.DEFAULT_ORF_LENGTH,
            error_bar = const.DEFAULT_ERROR_BAR):
    """
    Check if a set of consensus sequences in a FASTA file is intact.

    Args:
        input_folder: folder of files from NGS machine.

    Returns:
        Name of a file containing all consensus sequences.
    """

    pos_mapping = st.map_hxb2_positions_to_subtype(subtype)
    pos_subtype_mapping = {
        "forward": st.map_subtype_positions_to_hxb2("forward", subtype),
        "reverse": st.map_subtype_positions_to_hxb2("reverse", subtype),
    }

    # convert ORF positions to appropriate subtype
    forward_orfs, reverse_orfs, small_orfs = [
    [
        ExpectedORF.subtyped(pos_mapping, n, s, e, delta) \
        for (n, s, e, delta) in orfs
    ] \
    for orfs in [hxb2_forward_orfs, hxb2_reverse_orfs, hxb2_small_orfs]
    ]

    # convert PSI locus and RRE locus to appropriate subtype
    psi_locus = [pos_mapping[x] for x in hxb2_psi_locus]
    rre_locus = [pos_mapping[x] for x in hxb2_rre_locus]

    reference = st.subtype_sequence(subtype)
    blast_it = blast_iterate_inf(subtype, input_file) if check_internal_inversion or check_nonhiv or check_scramble else iterate_empty_lists()

    intact_file = os.path.join(working_dir, "intact.fasta")
    non_intact_file = os.path.join(working_dir, "nonintact.fasta")
    orf_file = os.path.join(working_dir, "orfs.json")
    error_file = os.path.join(working_dir, "errors.json")

    with open(intact_file, 'w') as intact_writer, \
         open(non_intact_file, 'w') as nonintact_writer, \
         open(orf_file, 'w') as orfs_writer, \
         open(error_file, 'w') as errors_writer:

        for (sequence, blast_rows) in with_blast_rows(blast_it, iterate_sequences(input_file)):
            sequence_errors = []

            reverse_sequence = SeqRecord.SeqRecord(Seq.reverse_complement(sequence.seq),
                                                   id = sequence.id + " [REVERSED]",
                                                   name = sequence.name
                                                   )

            alignment = wrappers.mafft([reference, sequence])
            reverse_alignment = wrappers.mafft([reference, reverse_sequence])

            forward_score = alignment_score(alignment)
            reverse_score = alignment_score(reverse_alignment)
            if alignment_score(reverse_alignment) > alignment_score(alignment):
                print("Reversing sequence " + sequence.id + "; forward score "
                      + str(forward_score) + "; reverse score " + str(reverse_score))
                alignment = reverse_alignment
                sequence = reverse_sequence

            sequence_orfs, orf_errors = has_reading_frames(
                alignment, sequence, reference, False,
                forward_orfs, error_bar)
            sequence_errors.extend(orf_errors)

            sequence_small_orfs, small_orf_errors = has_reading_frames(
                alignment, sequence, reference, True,
                small_orfs, error_bar, reverse = False)
            if include_small_orfs:
                sequence_errors.extend(small_orf_errors)

            hxb2_found_orfs = [ORF(
                o.name,
                o.orientation,
                pos_subtype_mapping[o.orientation][o.start],
                pos_subtype_mapping[o.orientation][o.end],
                o.distance,
                str(o.aminoseq),
                str(sequence[o.start:o.end].seq),
            ) for o in sorted(sequence_orfs + sequence_small_orfs, key=lambda o: o.start)]

            if include_packaging_signal:
                missing_psi_locus = has_packaging_signal(alignment,
                                                         psi_locus,
                                                         const.PSI_ERROR_TOLERANCE)
                if missing_psi_locus is not None:
                    sequence_errors.append(missing_psi_locus)

            if include_rre:
                missing_rre_locus = has_rev_response_element(alignment,
                                                             rre_locus,
                                                             const.RRE_ERROR_TOLERANCE
                                                             )
                if missing_rre_locus is not None:
                    sequence_errors.append(missing_rre_locus)

            if check_major_splice_donor_site:
                mutated_splice_donor_site = has_mutated_major_splice_donor_site(
                    alignment,
                    pos_mapping[hxb2_msd_site_locus],
                    pos_mapping[hxb2_msd_site_locus + 1],
                    const.DEFAULT_MSD_SEQUENCE)
                if mutated_splice_donor_site is not None:
                    sequence_errors.append(mutated_splice_donor_site)

            if run_hypermut is not None:
                hypermutated = isHypermut(alignment)

                if hypermutated is not None:
                    sequence_errors.append(hypermutated)

            if check_long_deletion is not None:
                long_deletion = has_long_deletion(sequence, alignment)
                if long_deletion:
                    sequence_errors.append(long_deletion)

            if check_nonhiv:
                error = is_nonhiv(sequence.id, len(sequence), blast_rows)
                if error:
                    sequence_errors.append(error)

            if check_scramble:
                error = is_scrambled(sequence.id, blast_rows)
                if error:
                    sequence_errors.append(error)

            if check_internal_inversion:
                error = contains_internal_inversion(sequence.id, blast_rows)
                if error:
                    sequence_errors.append(error)

            if len(sequence_errors) == 0:
                SeqIO.write([sequence], intact_writer, "fasta")
                intact_writer.flush()
            else:
                SeqIO.write([sequence], nonintact_writer, "fasta")
                nonintact_writer.flush()

            # add the small orf errors after the intactness check if not included
            if not include_small_orfs:
                sequence_errors.extend(small_orf_errors)

            orfs_writer.write(json.dumps({sequence.id: [x.__dict__ for x in hxb2_found_orfs]},
                                         indent=4))
            orfs_writer.flush()

            errors_writer.write(json.dumps({sequence.id: [x.__dict__ for x in sequence_errors]},
                                           indent=4))
            errors_writer.flush()

    return intact_file, non_intact_file, orf_file, error_file
#/end def intact
#/end intact.py
