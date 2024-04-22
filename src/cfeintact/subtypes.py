import os
from Bio import Seq, SeqIO, SeqRecord
from pathlib import Path

import importlib.resources as resources


def get_reference_dir() -> Path:
    ret = resources.path('cfeintact', 'subtype_alignments')
    if isinstance(ret, str):
        return Path(ret)
    elif isinstance(ret, Path):
        return ret
    else:
        return ret.__enter__()


REFERENCE_DIR = get_reference_dir()


# This variable holds a mapping between subtype names and their sequences.
# Prevents re-reading the same sequence from disk multiple time.
SEQUENCE_CACHE = {}


def HXB2():
    """
    Return the sequence of HXB2, the standard HIV reference sequence
    """
    return subtype_sequence("HXB2")


def subtypes():
    """
List all currently available HIV subtypes
"""
    return [f.replace(".fasta", "") for f in os.listdir(REFERENCE_DIR) if f.endswith(".fasta")]


def alignment_file(subtype):
    """
Return an alignment file associated with an HIV subtype.

Args:
    subtype: folder in which to put temporary files.
"""
    return os.path.join(REFERENCE_DIR, subtype + ".fasta")


def subtype_sequence(subtype):
    """
    Return an example sequence associated with an HIV subtype.

    Args:
        subtype: folder in which to put temporary files.
    """

    if subtype not in SEQUENCE_CACHE:
        alignment = list(SeqIO.parse(alignment_file(subtype), "fasta"))
        SEQUENCE_CACHE[subtype] = SeqRecord.SeqRecord(
            Seq.Seq(str(alignment[0].seq).replace("-", "").replace("\n", "")),
            id=alignment[0].id,
            name=alignment[0].name
        )

    return SEQUENCE_CACHE[subtype]


def convert_from_aligned_to_reference(position, alignment):
    hxb2_pos = 0
    subtype_pos = 0
    for i in range(len(alignment[0])):
        if subtype_pos == position:
            return hxb2_pos
        if alignment[0][i] != "-":
            subtype_pos += 1
        if alignment[1][i] != "-":
            hxb2_pos += 1
