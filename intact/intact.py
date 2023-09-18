import json
import os
import re
import csv
import dataclasses
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
import util.coordinates as coords
from util.blastrow import BlastRow


WRONGORFNUMBER_ERROR    = "WrongORFNumber"
MISPLACEDORF_ERROR      = "MisplacedORF"
LONGDELETION_ERROR      = "LongDeletion"
DELETIONINORF_ERROR     = "DeletionInOrf"
INSERTIONINORF_ERROR    = "InsertionInOrf"
INTERNALSTOP_ERROR      = "InternalStopInOrf"
SCRAMBLE_ERROR          = "Scramble"
NONHIV_ERROR            = "NonHIV"
INTERNALINVERSION_ERROR = "InternalInversion"
ALIGNMENT_FAILED        = "AlignmentFailed" # Happens when mafft process fails
INVALID_CODON           = "InvalidCodon" # Happens when there is an invalid codon in ORF

FRAMESHIFTINORF_ERROR   = "FrameshiftInOrf"
MSDMUTATED_ERROR        = "MajorSpliceDonorSiteMutated"
PSIDELETION_ERROR       = "PackagingSignalDeletion"
PSINOTFOUND_ERROR       = "PackagingSignalNotComplete"
RREDELETION_ERROR       = "RevResponseElementDeletion"
HYPERMUTATION_ERROR     = "APOBECHypermutationDetected"


@dataclass
class IntactnessError:
    sequence_name: str
    error: str
    message: str

@dataclass
class ExpectedORF:
    name: str
    start: int
    end: int
    deletion_tolerence: int
    nucleotides: str
    aminoacids: str
    protein: str

    @staticmethod
    def subtyped(reference, pos_mapping, name, start, end, deletion_tolerence):
        vpr_defective_insertion_pos = 5772
        start = start if start < vpr_defective_insertion_pos else start - 1
        end = end if end < vpr_defective_insertion_pos else end - 1

        start_s = pos_mapping[start - 1]
        end_s = pos_mapping[end]

        nucleotides = str(reference.seq[start_s:end_s])
        aminoacids = translate(nucleotides)
        has_start_codon = translate(reference.seq[(start - 1):end]).startswith("M")
        protein = get_biggest_protein(has_start_codon, aminoacids)

        return ExpectedORF(name=name,
                           start=start_s,
                           end=end_s,
                           deletion_tolerence=deletion_tolerence,
                           nucleotides=nucleotides,
                           aminoacids=aminoacids,
                           protein=protein,
                           )

@dataclass
class CandidateORF:
    name: str
    start: int
    end: int
    orientation: str
    distance: float
    protein: str
    aminoacids: str

@dataclass
class FoundORF:
    name: str
    orientation: str
    start: int
    end: int
    distance: str
    protein: str
    aminoacids: str
    nucleotides: str

@dataclass
class HolisticInfo:
    qlen: int = dataclasses.field(default=None)
    hypermutation_probablility: float = dataclasses.field(default=None)
    inferred_subtype: str  = dataclasses.field(default=None)
    blast_matched_qlen: int  = dataclasses.field(default=None) # number of query nucleotides matched to a known reference sequence
    blast_sseq_coverage: float  = dataclasses.field(default=None) # percentage of reference sequence covered by the query sequence
    blast_qseq_coverage: float  = dataclasses.field(default=None) # percentage of the query sequence covered by reference sequence
    blast_n_conseqs: int  = dataclasses.field(default=None) # number of blast conseqs in the resulting match


def iterate_values_from_csv(file_path):
    with open(file_path, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            yield row


def iterate_blast_rows_from_csv(file_path):
    previous_key = None
    values = []

    for row in iterate_values_from_csv(file_path):
        key = row['qseqid']
        typed = BlastRow.init(row)

        if key != previous_key and previous_key is not None:
            yield values
            values = []

        values.append(typed)
        previous_key = key

    if values:
        yield values


def blast_iterate(subtype, input_file, working_dir):
    with open(os.path.join(working_dir, 'blast.csv'), 'w') as output_file:
        db_file = st.alignment_file(subtype)
        wrappers.blast(db_file, input_file, output_file.name)
        for seq in iterate_blast_rows_from_csv(output_file.name):
            yield seq


def blast_iterate_inf(subtype, input_file, working_dir):
    for seq in blast_iterate(subtype, input_file, working_dir):
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


def is_nonhiv(holistic, seqid, seqlen, blast_rows):
    aligned_length = sum(abs(x.qend - x.qstart) + 1 for x in blast_rows)
    holistic.blast_matched_qlen = blast_rows[0].qlen if blast_rows else 1
    holistic.blast_qseq_coverage = aligned_length / holistic.blast_matched_qlen

    if holistic.blast_qseq_coverage < 0.8 and seqlen > holistic.blast_matched_qlen * 0.6:
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

def isHypermut(holistic, aln):
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

    holistic.hypermutation_probablility = 1 - pval

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


def get_biggest_protein(has_start_codon, aminoacids):
    def skip_to_startcodon(x):
        index = x.find("M")
        if index >= 0:
            return x[index:]
        else:
            return ""

    parts = aminoacids.split("*")
    subparts = [skip_to_startcodon(x) for x in parts] if has_start_codon else parts
    longest = max(subparts, key=len)
    return longest


aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2


def has_start_codon(orf):
    return orf.aminoacids[0] == "M"


def has_stop_codon(orf):
    return orf.aminoacids[-1] == "*"


def translate(seq, frame = 0, to_stop = False):
    for_translation = seq[frame:]
    for_translation += 'N' * ({0: 0, 1: 2, 2: 1}[len(for_translation) % 3])
    return Seq.translate(for_translation, to_stop = to_stop)


def has_reading_frames(
    alignment, sequence, reference, is_small,
    expected, error_bar, reverse = False
):
    """
    Check for presence of small reading frames
    """

    coordinates_mapping = coords.map_positions(alignment[0], alignment[1].seq)
    reverse_coordinates_mapping = coords.map_positions(alignment[1], alignment[0].seq)

    reference_aligned_mapping = coords.map_nonaligned_to_aligned_positions(reference, alignment[0].seq)
    query_aligned_mapping = coords.map_nonaligned_to_aligned_positions(sequence, alignment[1].seq)

    errors = []
    matches = []

    try:
        query_aminoacids_table = [translate(sequence.seq, i) for i in range(3)]
    except Bio.Data.CodonTable.TranslationError as e:
        log.error(e)
        err = IntactnessError(sequence.id, INVALID_CODON, e.args[0])
        errors.append(err)
        return matches, errors

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
        expected_nucleotides = e.nucleotides
        expected_aminoacids = e.aminoacids
        q_start = coordinates_mapping[e.start]
        q_end = coordinates_mapping[e.end - 1 if e.end == len(coordinates_mapping) else e.end]
        got_nucleotides = sequence.seq[q_start:q_end]
        got_aminoacids = translate(got_nucleotides)
        q_start_a = q_start // 3
        q_end_a = q_end // 3
        n = len(sequence.seq) - 1

        for frame in range(3):
            aminoacids = query_aminoacids_table[frame]
            for start_direction in [-1, +1]:
                for end_direction in [-1, +1]:
                    closest_start_a = q_start_a if not has_start_codon(e) else find_closest(aminoacids, q_start_a, start_direction, 'M')
                    closest_end_a = q_end_a if not has_stop_codon(e) else find_closest(aminoacids, q_end_a, end_direction, '*')
                    got_aminoacids = aminoacids[closest_start_a:closest_end_a + 1]
                    dist = 1 - jaro_similarity(got_aminoacids, expected_aminoacids)
                    closest_start = min(n, (closest_start_a * 3) + frame)
                    closest_end = min(n + 1, (closest_end_a * 3) + 3 + frame)
                    got_protein = get_biggest_protein(has_start_codon(e), got_aminoacids)
                    yield CandidateORF(e.name, closest_start, closest_end, "forward",
                                       dist, got_protein, got_aminoacids)

    def find_real_correspondence(e):
        q_start = coordinates_mapping[e.start]
        q_end = coordinates_mapping[e.end - 1 if e.end == len(coordinates_mapping) else e.end]
        candidates = list(find_candidate_positions(e, q_start, q_end))
        return min(candidates, key=lambda x: x.distance)

    def get_indel_impact(alignment):
        shift = 0
        impacted = 0
        counter = 0
        good = 0
        for (x, y) in zip(alignment[0], alignment[1]):
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

    for e in expected:
        best_match = find_real_correspondence(e)
        matches.append(best_match)

        got_protein = best_match.protein
        exp_protein = e.protein

        deletions = max(0, len(exp_protein) - len(got_protein)) * 3
        insertions = max(0, len(got_protein) - len(exp_protein)) * 3

        if got_protein and exp_protein:
            orf_alignment = aligner.align(exp_protein, got_protein)[0]
            best_match.distance = aligner.match_score - (orf_alignment.score / len(exp_protein))
        else:
            orf_alignment = (exp_protein, "-" * len(exp_protein))
            best_match.distance = aligner.match_score - ((aligner.open_gap_score + len(exp_protein) * aligner.extend_gap_score) / max(1, len(exp_protein)))

        # Max deletion allowed in ORF exceeded
        if deletions > e.deletion_tolerence:

            limit = e.deletion_tolerence // 3
            limited_aminoacids = best_match.aminoacids[limit:-limit]

            if "*" in limited_aminoacids:
                errors.append(IntactnessError(
                    sequence.id, INTERNALSTOP_ERROR,
                    ("Smaller " if is_small else "")
                    + "ORF " + str(e.name) + " at " + str(e.start)
                    + "-" + str(e.end)
                    + " contains an internal stop codon at " + str(e.start + (limit + limited_aminoacids.index("*")) * 3)
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

        got_nucleotides = sequence.seq[best_match.start:best_match.start + len(got_protein) * 3].upper()
        exp_nucleotides = reference.seq[e.start:e.end].upper()
        if got_nucleotides and exp_nucleotides:
            orf_alignment = aligner.align(exp_nucleotides, got_nucleotides)[0]
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


class OutputWriter:
    def __init__(self, working_dir, fmt):
        self.fmt = fmt
        self.working_dir = working_dir
        self.intact_path = os.path.join(working_dir, "intact.fasta")
        self.non_intact_path = os.path.join(working_dir, "nonintact.fasta")
        self.orf_path = os.path.join(working_dir, f"orfs.{fmt}")
        self.holistic_path = os.path.join(working_dir, f"holistic.{fmt}")
        self.error_path = os.path.join(working_dir, f"errors.{fmt}")

        if fmt not in ("json", "csv"):
            raise ValueError(f"Unrecognized output format {fmt}")


    def __enter__(self, *args):
        self.intact_file = open(self.intact_path, 'w')
        self.nonintact_file = open(self.non_intact_path, 'w')
        self.orfs_file = open(self.orf_path, 'w')
        self.holistic_file = open(self.holistic_path, 'w')
        self.errors_file = open(self.error_path, 'w')

        if self.fmt == "json":
            self.orfs = {}
            self.holistic = {}
            self.errors = {}
        elif self.fmt == "csv":
            self.orfs_writer = csv.writer(self.orfs_file)
            self.orfs_header = ['seqid'] + [field.name for field in dataclasses.fields(CandidateORF)]
            self.orfs_writer.writerow(self.orfs_header)
            self.holistic_writer = csv.writer(self.holistic_file)
            self.holistic_header = ['seqid'] + [field.name for field in dataclasses.fields(HolisticInfo)]
            self.holistic_writer.writerow(self.holistic_header)
            self.errors_writer = csv.writer(self.errors_file)
            self.errors_header = [field.name for field in dataclasses.fields(IntactnessError)]
            self.errors_writer.writerow(self.errors_header)

        return self


    def __exit__(self, *args):
        if self.fmt == "json":
            json.dump(self.orfs, self.orfs_file, indent=4)
            json.dump(self.holistic, self.holistic_file, indent=4)
            json.dump(self.errors, self.errors_file, indent=4)

        self.intact_file.close()
        self.nonintact_file.close()
        self.orfs_file.close()
        self.holistic_file.close()
        self.errors_file.close()

        log.info('Intact sequences written to ' + self.intact_path)
        log.info('Non-intact sequences written to ' + self.non_intact_path)
        log.info('ORFs for all sequences written to ' + self.orf_path)
        log.info('Holistic info for all sequences written to ' + self.holistic_path)
        log.info('Intactness error information written to ' + self.error_path)
        if os.path.exists(os.path.join(self.working_dir, 'blast.csv')):
            log.info('Blast output written to ' + os.path.join(self.working_dir, 'blast.csv'))


    def write(self, sequence, is_intact, orfs, errors, holistic):
        fasta_file = self.intact_file if is_intact else self.nonintact_file
        SeqIO.write([sequence], fasta_file, "fasta")
        fasta_file.flush()

        if self.fmt == "json":
            self.orfs[sequence.id] = orfs
            self.holistic[sequence.id] = holistic
            self.errors[sequence.id] = errors
        elif self.fmt == "csv":
            for orf in orfs:
                self.orfs_writer.writerow([(sequence.id if key == 'seqid' else orf[key]) for key in self.orfs_header])
            self.holistic_writer.writerow([(sequence.id if key == 'seqid' else holistic[key]) for key in self.holistic_header])
            for error in errors:
                self.errors_writer.writerow([error[key] for key in self.errors_header])

def strip_sequence_dashes(seq):
    return SeqRecord.SeqRecord(
        Seq.Seq(str(seq.seq).replace("-","").replace("\n", "")),
        id = seq.id, name = seq.name)

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
            output_csv,
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

    subtype_choices = {}
    with open(st.alignment_file(subtype), 'r') as in_handle:
        for sequence in SeqIO.parse(in_handle, "fasta"):
            subtype_choices[sequence.id] = strip_sequence_dashes(sequence)

    hxb2_reference = st.HXB2()

    with OutputWriter(working_dir, "csv" if output_csv else "json") as writer:

        blast_it = blast_iterate_inf(subtype, input_file, working_dir) if check_internal_inversion or check_nonhiv or check_scramble or 1 < len(subtype_choices) else iterate_empty_lists()
        for (sequence, blast_rows) in with_blast_rows(blast_it, iterate_sequences(input_file)):

            blast_rows_statistics = {}
            for blast_row in blast_rows:
                blast_rows_statistics[blast_row.sseqid] = blast_rows_statistics.get(blast_row.sseqid, 0) + abs(blast_row.qend - blast_row.qstart)

            if blast_rows_statistics:
                reference_name = max(blast_rows_statistics, key=blast_rows_statistics.get)
            else:
                reference_name = sorted(subtype_choices.keys())[0]

            reference = subtype_choices[reference_name]
            hxb2_alignment = wrappers.mafft([hxb2_reference, reference])
            pos_mapping = coords.map_positions(hxb2_alignment[0], hxb2_alignment[1])

            # convert ORF positions to appropriate subtype
            forward_orfs, reverse_orfs, small_orfs = [
            [
                ExpectedORF.subtyped(reference, pos_mapping, n, s, e, delta) \
                for (n, s, e, delta) in orfs
            ] \
            for orfs in [hxb2_forward_orfs, hxb2_reverse_orfs, hxb2_small_orfs]
            ]

            # convert PSI locus and RRE locus to appropriate subtype
            psi_locus = [pos_mapping[x] for x in hxb2_psi_locus]
            rre_locus = [pos_mapping[x] for x in hxb2_rre_locus]

            sequence_errors = []
            holistic = HolisticInfo()
            holistic.qlen = len(sequence)
            holistic.inferred_subtype = reference_name
            holistic.blast_n_conseqs = len(blast_rows)

            aligned_reference_length = sum(abs(x.send - x.sstart) + 1 for x in blast_rows)
            blast_matched_slen = blast_rows[0].slen if blast_rows else 1
            holistic.blast_sseq_coverage = aligned_reference_length / blast_matched_slen

            reverse_sequence = SeqRecord.SeqRecord(Seq.reverse_complement(sequence.seq),
                                                   id = sequence.id + " [REVERSED]",
                                                   name = sequence.name
                                                   )

            try:
                alignment = wrappers.mafft([reference, sequence])
                reverse_alignment = wrappers.mafft([reference, reverse_sequence])
            except wrappers.AlignmentFailure:
                err = IntactnessError(sequence.id, ALIGNMENT_FAILED, "Alignment failed for this sequence. It probably contains invalid symbols.")
                sequence_errors.append(err)
                errors = [x.__dict__ for x in sequence_errors]
                holistic = holistic.__dict__
                writer.write(sequence, is_intact=False, orfs=[], errors=errors, holistic=holistic)
                continue

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

            hxb2_found_orfs = [FoundORF(
                o.name,
                o.orientation,
                o.start,
                o.end,
                o.distance,
                str(o.protein),
                str(o.aminoacids),
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
                hypermutated = isHypermut(holistic, alignment)

                if hypermutated is not None:
                    sequence_errors.append(hypermutated)

            if check_long_deletion is not None:
                long_deletion = has_long_deletion(sequence, alignment)
                if long_deletion:
                    sequence_errors.append(long_deletion)

            if check_nonhiv:
                error = is_nonhiv(holistic, sequence.id, len(sequence), blast_rows)
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

            is_intact = len(sequence_errors) == 0

            # add the small orf errors after the intactness check if not included
            if not include_small_orfs:
                sequence_errors.extend(small_orf_errors)

            orfs = [x.__dict__ for x in hxb2_found_orfs]
            errors = [x.__dict__ for x in sequence_errors]
            holistic = holistic.__dict__
            writer.write(sequence, is_intact, orfs, errors, holistic=holistic)
#/end def intact
#/end intact.py
