import json
import os
import re
import csv
import dataclasses
from dataclasses import dataclass
from collections import Counter
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData
from scipy.stats import fisher_exact
from typing import Optional, Dict, List, Iterable, Union, Tuple

import cfeintact.constants as const
import cfeintact.subtypes as st
import cfeintact.wrappers as wrappers
import cfeintact.log as log
import cfeintact.defect as defect
from cfeintact.defect import Defect, ORFDefect
from cfeintact.aligned_sequence import AlignedSequence
from cfeintact.blastrow import BlastRow
from cfeintact.initialize_orf import initialize_orf
from cfeintact.original_orf import OriginalORF
from cfeintact.mapped_orf import MappedORF
from cfeintact.find_orf import find_orf


@dataclass(frozen=True)
class FoundORF:
    name: str
    start: int
    end: int
    orientation: str
    distance: float
    indel_impact: int
    protein: str
    aminoacids: str
    nucleotides: str
    subtype_start: int
    subtype_end: int
    subtype_aminoacids: str
    subtype_nucleotides: str


@dataclass
class HolisticInfo:
    intact: Optional[bool] = dataclasses.field(default=None)
    qlen: Optional[int] = dataclasses.field(default=None)
    hypermutation_probablility: Optional[float] = dataclasses.field(default=None)
    inferred_subtype: Optional[str] = dataclasses.field(default=None)
    # number of query nucleotides matched to a known reference sequence
    blast_matched_qlen: Optional[int] = dataclasses.field(default=None)
    # percentage of reference sequence covered by the query sequence
    blast_sseq_coverage: Optional[float] = dataclasses.field(default=None)
    # percentage of the query sequence covered by reference sequence
    blast_qseq_coverage: Optional[float] = dataclasses.field(default=None)
    # percentage of the query sequence covered by reference sequence
    blast_sseq_orfs_coverage: Optional[float] = dataclasses.field(default=None)
    # start position of the region used for orfs coverage
    orfs_start: Optional[int] = dataclasses.field(default=None)
    # end position of the region used for orfs coverage
    orfs_end: Optional[int] = dataclasses.field(default=None)
    # number of blast conseqs in the resulting match
    blast_n_conseqs: Optional[int] = dataclasses.field(default=None)


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
        with st.alignment_file(subtype) as db_file:
            wrappers.blast(str(db_file), input_file, output_file.name)

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


def contains_internal_inversion(qseqid: str, blast_rows: List[BlastRow]) -> Optional[Defect]:
    ignored_5_prime = remove_5_prime(blast_rows)

    if not ignored_5_prime:
        # No alignment.
        # It should be an error normally, yet not an internal inversion error.
        return None

    all_same = len(set(x.sstrand for x in ignored_5_prime)) == 1
    if all_same:
        return None

    # Some parts of the sequence were aligned
    # in forward direction (plus)
    # and some in reverse (minus).
    # This indicates an internal inversion.
    return Defect(qseqid, defect.InternalInversion())


def is_scrambled(qseqid, blast_rows):
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
        return Defect(qseqid, defect.Scramble(direction))


def is_nonhiv(holistic, qseqid, seqlen, blast_rows):
    aligned_length = sum(abs(x.qend - x.qstart) + 1 for x in blast_rows)
    holistic.blast_matched_qlen = blast_rows[0].qlen if blast_rows else 1
    holistic.blast_qseq_coverage = aligned_length / holistic.blast_matched_qlen

    if holistic.blast_qseq_coverage < 0.8 and seqlen > holistic.blast_matched_qlen * 0.6:
        return Defect(qseqid, defect.NonHIV())
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
# /end _getPositions


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
        alternative='greater'
    )

    holistic.hypermutation_probablility = 1 - pval

    if pval < 0.05:
        return Defect(
            aln[1].id,
            defect.APOBECHypermutationDetected(pval),
        )
    # /end if
    return None
# /end isHypermut


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
        return Defect(sequence.id, defect.LongDeletion())

    return None
# /end has_long_deletion


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
    context = alignment[1].seq[(sd_begin - 10):(sd_end + 1 + 10)]

    if sd.upper() != splice_donor_sequence.upper():
        return Defect(
            alignment[1].id,
            defect.MajorSpliceDonorSiteMutated(''.join(sd.upper()), ''.join(context.upper())),
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

        PSINOTFOUND_ERROR -- Defect denoting query does not encompass
                             the complete PSI region.
        PSIDELETION_ERROR -- Defect denoting query likely contains
                             defective PSI.
        None -- Denotes intact PSI.
    """
    packaging_begin = [m.start() for m in re.finditer(r"[^-]",
                       str(alignment[0].seq))][psi_locus[0]]
    packaging_end = [m.start() for m in re.finditer(r"[^-]",
                     str(alignment[0].seq))][psi_locus[1]]
    # query_options = [m and m.start() for m in re.finditer(r"[^-]",
    #                  str(alignment[1].seq))]
    # query_start = query_options[0] if query_options else ''
    # if query_start > packaging_begin:
    #     return Defect(
    #             alignment[1].id, PSINOTFOUND_ERROR,
    #             "Query Start at reference position " + str(query_start)
    #             + ". Does not encompass PSI at positions "
    #             + str(packaging_begin) + " to " + str(packaging_end) + "."
    #             )
    # #/end if
    query_psi = str(alignment[1].seq[packaging_begin:packaging_end])
    query_psi_deletions = len(re.findall(r"-", query_psi))
    if query_psi_deletions > psi_tolerance:
        return Defect(
            alignment[1].id,
            defect.PackagingSignalDeletion(query_psi_deletions, psi_tolerance),
        )
    # /end if
    return None
# /end def has_packaging_signal


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

        RREDELETION_ERROR -- Defect denoting query likely contains
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
        return Defect(
            alignment[1].id,
            defect.RevResponseElementDeletion(query_rre_deletions, rre_tolerance),
        )
    # /end if
    return None
# /end def has_rev_response_element


def check_reading_frame_shift(reference: SeqRecord,
                              sequence: SeqRecord,
                              best_match: MappedORF) \
        -> Optional[defect.FrameshiftInOrf]:

    q = best_match.query
    e = best_match.reference

    impacted_by_indels = best_match.indel_impact
    if impacted_by_indels >= max(e.max_deletions, e.max_insertions) + 3:
        return defect.FrameshiftInOrf(e=e, q=q, impacted_positions=impacted_by_indels)

    return None


def check_reading_frame_deletions(e: OriginalORF, q: OriginalORF) \
        -> Optional[Union[defect.InternalStopInOrf, defect.DeletionInOrf]]:

    got_protein = q.protein
    exp_protein = e.protein
    deletions = max(0, len(exp_protein) - len(got_protein)) * 3

    # Max deletion allowed in ORF exceeded
    if deletions > e.max_deletions:

        limit = e.max_deletions // 3
        limited_aminoacids = q.aminoacids[limit:-limit]

        if "*" in limited_aminoacids:
            position = q.start + (limit + limited_aminoacids.index('*')) * 3
            return defect.InternalStopInOrf(e=e, q=q, position=position)
        else:
            return defect.DeletionInOrf(e=e, q=q, deletions=deletions)

    return None


def check_reading_frame_insertions(check_distance: bool, best_match: MappedORF, e: OriginalORF, q: OriginalORF) \
        -> Optional[defect.InsertionInOrf]:

    if not check_distance:
        return None

    if best_match.distance <= e.max_distance:
        return None

    got_protein = q.protein
    exp_protein = e.protein
    insertions = max(0, len(got_protein) - len(exp_protein)) * 3

    # Max insertions allowed in ORF exceeded
    if insertions > e.max_insertions:
        return defect.InsertionInOrf(e=e, q=q, insertions=insertions)

    return None


def check_reading_frame_distance(check_distance: bool, best_match: MappedORF, e: OriginalORF, q: OriginalORF) \
        -> Optional[defect.SequenceDivergence]:

    if not check_distance:
        return None

    if best_match.distance > e.max_distance:
        return defect.SequenceDivergence(q=q, e=e, distance=best_match.distance)
    else:
        return None


def check_reading_frame(check_distance: bool,
                        aligned_sequence: AlignedSequence,
                        expected: Iterable[OriginalORF],
                        reverse: bool = False) \
        -> Tuple[List[MappedORF], List[defect.Defect]]:

    sequence = aligned_sequence.this
    reference = aligned_sequence.reference
    errors = []
    matches = []

    def add_error(defect):
        if defect:
            assert sequence.id is not None
            errors.append(Defect(sequence.id, defect))

    for e in expected:
        best_match = find_orf(aligned_sequence, e)
        matches.append(best_match)
        q = best_match.query

        add_error(check_reading_frame_shift(reference, sequence=sequence, best_match=best_match))
        add_error(check_reading_frame_deletions(e=e, q=q))
        add_error(check_reading_frame_insertions(check_distance, best_match, e=e, q=q))
        add_error(check_reading_frame_distance(check_distance, best_match, e=e, q=q))

    return matches, errors


def iterate_sequences(input_file):
    with open(input_file, 'r') as in_handle:
        for sequence in SeqIO.parse(in_handle, "fasta"):
            yield sequence.upper()


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
        self.subtypes_path = os.path.join(working_dir, "subtypes.fasta")
        self.orf_path = os.path.join(working_dir, f"orfs.{fmt}")
        self.holistic_path = os.path.join(working_dir, f"holistic.{fmt}")
        self.error_path = os.path.join(working_dir, f"errors.{fmt}")

        self.subtypes = set()

        if fmt not in ("json", "csv"):
            raise ValueError(f"Unrecognized output format {fmt}")

    def __enter__(self, *args):
        self.intact_file = open(self.intact_path, 'w')
        self.nonintact_file = open(self.non_intact_path, 'w')
        self.subtypes_file = open(self.subtypes_path, 'w')
        self.orfs_file = open(self.orf_path, 'w')
        self.holistic_file = open(self.holistic_path, 'w')
        self.errors_file = open(self.error_path, 'w')

        if self.fmt == "json":
            self.orfs = {}
            self.holistic = {}
            self.errors = {}
        elif self.fmt == "csv":
            self.orfs_writer = csv.writer(self.orfs_file)
            self.orfs_header = ['qseqid'] + [field.name for field in dataclasses.fields(FoundORF)]
            self.orfs_writer.writerow(self.orfs_header)
            self.holistic_writer = csv.writer(self.holistic_file)
            self.holistic_header = ['qseqid'] + [field.name for field in dataclasses.fields(HolisticInfo)]
            self.holistic_writer.writerow(self.holistic_header)
            self.errors_writer = csv.DictWriter(
                self.errors_file, fieldnames=["qseqid", "error", "message", "orf"])
            self.errors_writer.writeheader()

        return self

    def __exit__(self, *args):
        if self.fmt == "json":
            json.dump(self.orfs, self.orfs_file, indent=4)
            json.dump(self.holistic, self.holistic_file, indent=4)
            json.dump(self.errors, self.errors_file, indent=4)

        self.intact_file.close()
        self.nonintact_file.close()
        self.subtypes_file.close()
        self.orfs_file.close()
        self.holistic_file.close()
        self.errors_file.close()

        log.info('Intact sequences written to ' + self.intact_path)
        log.info('Non-intact sequences written to ' + self.non_intact_path)
        log.info('Subtype sequences written to ' + self.subtypes_path)
        log.info('ORFs for all sequences written to ' + self.orf_path)
        log.info('Holistic info for all sequences written to ' + self.holistic_path)
        log.info('Intactness error information written to ' + self.error_path)
        if os.path.exists(os.path.join(self.working_dir, 'blast.csv')):
            log.info('Blast output written to ' + os.path.join(self.working_dir, 'blast.csv'))

    def write(self, sequence, subtype, is_intact, orfs, defects, holistic):
        fasta_file = self.intact_file if is_intact else self.nonintact_file
        SeqIO.write([sequence], fasta_file, "fasta")
        fasta_file.flush()

        if subtype.id not in self.subtypes:
            self.subtypes.add(subtype.id)
            SeqIO.write([subtype], self.subtypes_file, "fasta")
            self.subtypes_file.flush()

        errors_dicts = [{
            "qseqid": d.qseqid,
            "error": d.error.__class__.__name__,
            "message": str(d.error),
            "orf": d.error.q.name if isinstance(d.error, ORFDefect) else None,
        } for d in defects]

        if self.fmt == "json":
            self.orfs[sequence.id] = orfs
            self.holistic[sequence.id] = holistic
            self.errors[sequence.id] = errors_dicts
        elif self.fmt == "csv":
            for orf in orfs:
                row = [(sequence.id if key == 'qseqid' else orf[key]) for key in self.orfs_header]
                self.orfs_writer.writerow(row)
            self.holistic_writer.writerow(
                [(sequence.id if key == 'qseqid' else holistic[key]) for key in self.holistic_header])
            for error in errors_dicts:
                self.errors_writer.writerow(error)


def read_hxb2_orfs(aligned_subtype: AlignedSequence,
                   orfs: const.ORFsDefinition,
                   is_small: bool) -> Iterable[OriginalORF]:
    for (name, start, end, max_deletions, max_insertions, max_distance) in orfs:
        # Decrement is needed because original coordinates are 1-based.
        start = start - 1
        end = end - 1

        # Offset by the vpr bug in the original HXB2
        vpr_defective_insertion_pos = 5771
        start = start if start < vpr_defective_insertion_pos else start - 1
        end = end if end < vpr_defective_insertion_pos else end - 1

        yield initialize_orf(aligned_subtype, name=name,
                             start=start, end=end,
                             max_deletions=max_deletions,
                             max_insertions=max_insertions,
                             max_distance=max_distance,
                             is_small=is_small)


VALID_DNA_CHARACTERS = IUPACData.ambiguous_dna_letters.upper()


def find_invalid_subsequences(sequence):
    invalid_subsequences = []
    current_subsequence = []
    start_position = None

    for idx, symbol in enumerate(sequence.seq):
        if symbol not in VALID_DNA_CHARACTERS:
            if start_position is None:
                start_position = idx
            current_subsequence.append(symbol)
        else:
            if current_subsequence:
                end_position = idx - 1
                invalid_subsequences.append(
                    {
                        'sequence': ''.join(current_subsequence),
                        'start': start_position,
                        'end': end_position
                    }
                )
                current_subsequence = []
                start_position = None

    if current_subsequence:
        end_position = len(sequence.seq) - 1
        invalid_subsequences.append(
            {
                'sequence': ''.join(current_subsequence),
                'start': start_position,
                'end': end_position
            }
        )

    return invalid_subsequences


def intact(working_dir: str,
           input_file: str,
           subtype: str,
           check_packaging_signal: bool,
           check_rre: bool,
           check_major_splice_donor_site: bool,
           check_hypermut: bool,
           check_long_deletion: bool,
           check_nonhiv: bool,
           check_scramble: bool,
           check_internal_inversion: bool,
           check_unknown_nucleotides: bool,
           check_small_orfs: bool,
           check_distance: bool,
           output_csv: bool,
           hxb2_forward_orfs: const.ORFsDefinition = const.DEFAULT_FORWARD_ORFs,
           hxb2_reverse_orfs: const.ORFsDefinition = const.DEFAULT_REVERSE_ORFS,
           hxb2_small_orfs: const.ORFsDefinition = const.DEFAULT_SMALL_FORWARD_ORFS,
           hxb2_psi_locus: Tuple[int, int] = const.DEFAULT_PSI_LOCUS,
           hxb2_rre_locus: Tuple[int, int] = const.DEFAULT_RRE_LOCUS,
           hxb2_msd_site_locus: int = const.DEFAULT_MSD_SITE_LOCUS,
           ) -> None:
    """
    Check if a set of consensus sequences in a FASTA file is cfeintact.

    Args:
        input_folder: folder of files from NGS machine.

    Returns:
        Name of a file containing all consensus sequences.
    """

    subtype_choices = {}
    with st.alignment_file(subtype) as path:
        with open(path, 'r') as in_handle:
            for sequence in SeqIO.parse(in_handle, "fasta"):
                subtype_choices[sequence.id] = sequence

    def analyse_single_sequence(writer, sequence, blast_rows):
        sequence_errors = []
        holistic = HolisticInfo()

        invalid_subsequences = find_invalid_subsequences(sequence)
        if invalid_subsequences:

            if check_unknown_nucleotides:
                error_details = ', '.join(
                    f"{subseq['sequence']} (start: {subseq['start']}, end: {subseq['end']})"
                    for subseq in invalid_subsequences
                )
                err = Defect(sequence.id, defect.UnknownNucleotide(error_details))
                sequence_errors.append(err)

            seq = Seq.Seq(''.join(x for x in sequence.seq if x in VALID_DNA_CHARACTERS))
            sequence = SeqRecord(seq, id=sequence.id, name=sequence.name)

        holistic.qlen = len(sequence)
        holistic.blast_n_conseqs = len(blast_rows)

        aligned_reference_length = sum(
            abs(x.send - x.sstart) + 1 for x in blast_rows)
        blast_matched_slen = blast_rows[0].slen if blast_rows else 1
        holistic.blast_sseq_coverage = aligned_reference_length / blast_matched_slen

        blast_rows_statistics: Dict[str, int] = {}
        for blast_row in blast_rows:
            blast_rows_statistics[blast_row.sseqid] = blast_rows_statistics.get(
                blast_row.sseqid, 0) + abs(blast_row.qend - blast_row.qstart)

        if blast_rows_statistics:
            reference_name = max(blast_rows_statistics, key=lambda key: blast_rows_statistics.get(key, 0))
        else:
            reference_name = sorted(subtype_choices.keys())[0]

        blast_orientation_statistics = {"plus": 0, "minus": 0}
        for blast_row in blast_rows:
            if blast_row.qseqid == reference_name:
                blast_orientation_statistics[blast_row.sstrand] += abs(
                    blast_row.qend - blast_row.qstart)

        holistic.inferred_subtype = reference_name
        reference = subtype_choices[reference_name]
        aligned_subtype = AlignedSequence(this=reference, reference=st.HXB2())

        forward_aligned_sequence = AlignedSequence(
            this=sequence, reference=aligned_subtype.this)
        reverse_aligned_sequence = forward_aligned_sequence.reverse()

        if blast_orientation_statistics["minus"] < blast_orientation_statistics["plus"] \
           or reverse_aligned_sequence.alignment_score() <= forward_aligned_sequence.alignment_score():
            aligned_sequence = forward_aligned_sequence
        else:
            log.info("Reversing sequence " + sequence.id)
            aligned_sequence = reverse_aligned_sequence
            sequence = aligned_sequence.this

        # convert ORF positions to appropriate subtype
        forward_orfs, reverse_orfs = \
            [list(read_hxb2_orfs(aligned_subtype, orfs, is_small=False))
             for orfs in [hxb2_forward_orfs, hxb2_reverse_orfs]]
        small_orfs = list(read_hxb2_orfs(aligned_subtype, hxb2_small_orfs, is_small=True))

        holistic.orfs_start = min(forward_orfs, key=lambda e: e.start).start
        holistic.orfs_end = max(forward_orfs, key=lambda e: e.end).end

        def clamp(p):
            return max(min(p, holistic.orfs_end), holistic.orfs_start)
        aligned_reference_orfs_length = sum(abs(clamp(x.send + 1) - clamp(x.sstart)) for x in blast_rows)
        if holistic.orfs_end is not None and holistic.orfs_start is not None:
            blast_matched_orfs_slen = holistic.orfs_end - holistic.orfs_start
        else:
            blast_matched_orfs_slen = 0
        holistic.blast_sseq_orfs_coverage = aligned_reference_orfs_length / blast_matched_orfs_slen

        # convert PSI locus to appropriate subtype
        hxb2_psi_locus_start, hxb2_psi_locus_end = hxb2_psi_locus
        psi_locus_start = aligned_subtype.coordinate_mapping.ref_to_query.right_min(hxb2_psi_locus_start) or 0
        psi_locus_end = aligned_subtype.coordinate_mapping.ref_to_query.left_max(hxb2_psi_locus_end) or 0
        psi_locus = (psi_locus_start, psi_locus_end)

        # convert RRE locus to appropriate subtype
        hxb2_rre_locus_start, hxb2_rre_locus_end = hxb2_rre_locus
        rre_locus_start = aligned_subtype.coordinate_mapping.ref_to_query.right_min(hxb2_rre_locus_start) or 0
        rre_locus_end = aligned_subtype.coordinate_mapping.ref_to_query.left_max(hxb2_rre_locus_end) or 0
        rre_locus = (rre_locus_start, rre_locus_end)

        # convert MSD site to appropriate subtype
        hxb2_msd_site_locus_start, hxb2_msd_site_locus_end = (hxb2_msd_site_locus, hxb2_msd_site_locus + 1)
        msd_site_locus_start = aligned_subtype.coordinate_mapping.ref_to_query.right_min(hxb2_msd_site_locus_start) or 0
        msd_site_locus_end = aligned_subtype.coordinate_mapping.ref_to_query.left_max(hxb2_msd_site_locus_end) or 0

        alignment = aligned_sequence.alignment

        sequence_orfs, orf_errors = check_reading_frame(check_distance, aligned_sequence, forward_orfs)
        sequence_errors.extend(orf_errors)

        sequence_small_orfs, small_orf_errors = check_reading_frame(check_distance, aligned_sequence, small_orfs)
        if check_small_orfs:
            sequence_errors.extend(small_orf_errors)

        hxb2_found_orfs = [FoundORF(
            o.query.name,
            o.query.start,
            o.query.end,
            o.orientation,
            o.distance,
            o.indel_impact,
            str(o.query.protein),
            str(o.query.aminoacids),
            str(o.query.nucleotides),
            o.reference.start,
            o.reference.end,
            str(o.reference.aminoacids),
            str(o.reference.nucleotides),
        ) for o in sorted(sequence_orfs + sequence_small_orfs, key=lambda o: o.query.start)]

        if check_packaging_signal:
            missing_psi_locus = has_packaging_signal(alignment,
                                                     psi_locus,
                                                     const.PSI_ERROR_TOLERANCE)
            if missing_psi_locus is not None:
                sequence_errors.append(missing_psi_locus)

        if check_rre:
            missing_rre_locus = has_rev_response_element(alignment,
                                                         rre_locus,
                                                         const.RRE_ERROR_TOLERANCE
                                                         )
            if missing_rre_locus is not None:
                sequence_errors.append(missing_rre_locus)

        if check_major_splice_donor_site:
            mutated_splice_donor_site = has_mutated_major_splice_donor_site(
                alignment,
                msd_site_locus_start,
                msd_site_locus_end,
                const.DEFAULT_MSD_SEQUENCE)
            if mutated_splice_donor_site is not None:
                sequence_errors.append(mutated_splice_donor_site)

        if check_hypermut is not None:
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

        holistic.intact = len(sequence_errors) == 0

        # add the small orf errors after the intactness check if not included
        if not check_small_orfs:
            sequence_errors.extend(small_orf_errors)

        orfs = [x.__dict__ for x in hxb2_found_orfs]
        subtype = aligned_sequence.reference
        writer.write(sequence, subtype, holistic.intact,
                     orfs, sequence_errors, holistic.__dict__)

    with OutputWriter(working_dir, "csv" if output_csv else "json") as writer:

        should_run_blast = check_internal_inversion or check_nonhiv or check_scramble or 1 < len(subtype_choices)
        blast_it = blast_iterate_inf(subtype, input_file, working_dir) if should_run_blast else iterate_empty_lists()
        for (sequence, blast_rows) in with_blast_rows(blast_it, iterate_sequences(input_file)):
            analyse_single_sequence(writer, sequence, blast_rows)

# /end def intact
# /end cfeintact.py
