import os
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Iterator, List, Sequence, Optional, Dict
import contextlib

import importlib.resources as resources


@contextlib.contextmanager
def reference_dir() -> Iterator[Path]:
    ret = resources.path('cfeintact', 'subtype_alignments')
    if isinstance(ret, str):
        yield Path(ret)
    elif isinstance(ret, Path):
        yield ret
    else:
        with ret as path:
            yield path


# This variable holds a mapping between subtype names and their sequences.
# Prevents re-reading the same sequence from disk multiple time.
SEQUENCE_CACHE: Dict[str, SeqRecord] = {}


def HXB2():
    """
    Return the sequence of HXB2, the standard HIV reference sequence
    """
    return subtype_sequence("HXB2")


def subtypes():
    """
List all currently available HIV subtypes
"""
    with reference_dir() as REFERENCE_DIR:
        return [f.replace(".fasta", "") for f in os.listdir(REFERENCE_DIR) if f.endswith(".fasta")]


@contextlib.contextmanager
def alignment_file(subtype: str) -> Iterator[Path]:
    with reference_dir() as REFERENCE_DIR:
        yield Path(os.path.join(REFERENCE_DIR, subtype + ".fasta"))


def alignment_sequence(subtype: str) -> List[SeqRecord]:
    """
Return an alignment file associated with an HIV subtype.

Args:
    subtype: folder in which to put temporary files.
"""
    with alignment_file(subtype) as path:
        return list(SeqIO.parse(path, "fasta"))


def subtype_sequence(subtype: str) -> SeqRecord:
    """
    Return an example sequence associated with an HIV subtype.

    Args:
        subtype: folder in which to put temporary files.
    """

    if subtype not in SEQUENCE_CACHE:
        # alignment = list(SeqIO.parse(alignment_file(subtype), "fasta"))
        alignment = alignment_sequence(subtype)
        SEQUENCE_CACHE[subtype] = SeqRecord(
            Seq.Seq(str(alignment[0].seq).replace("-", "").replace("\n", "")),
            id=alignment[0].id,
            name=alignment[0].name
        )

    return SEQUENCE_CACHE[subtype]


def convert_from_aligned_to_reference(position: int, alignment: List[Sequence[str]]) -> Optional[int]:
    hxb2_pos = 0
    subtype_pos = 0
    for i in range(len(alignment[0])):
        if subtype_pos == position:
            return hxb2_pos
        if alignment[0][i] != "-":
            subtype_pos += 1
        if alignment[1][i] != "-":
            hxb2_pos += 1

    return None
